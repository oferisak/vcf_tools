# Get chromosomes and their lengths from the first VCF file
get_chromosomes <- function(vcf_file) {
    cmd <- paste("bcftools view -h", shQuote(vcf_file), "| grep '^##contig'")
    message("\n=== Running commands to get chromosomes and lengths")
    message(cmd)
    contig_lines <- system(cmd, intern = TRUE)
    if (length(contig_lines) == 0) {
        # Fallback: query unique chromosomes in data
        cmd2 <- paste("bcftools query -f '%CHROM\\n'", shQuote(vcf_file), "| sort -u")
        message(cmd2)
        chroms <- system(cmd2, intern = TRUE)
        return(data.frame(chrom = chroms, length = NA, stringsAsFactors = FALSE))
    }
    chrom_info <- data.frame(chrom = character(), length = numeric(), stringsAsFactors = FALSE)
    for (line in contig_lines) {
        chrom_match <- regexpr("ID=([^,>]+)", line, perl = TRUE)
        length_match <- regexpr("length=([0-9]+)", line, perl = TRUE)
        if (chrom_match > 0) {
            chrom <- substr(
                line,
                attr(chrom_match, "capture.start")[1],
                attr(chrom_match, "capture.start")[1] + attr(chrom_match, "capture.length")[1] - 1
            )
            if (length_match > 0) {
                length_val <- as.numeric(substr(
                    line,
                    attr(length_match, "capture.start")[1],
                    attr(length_match, "capture.start")[1] + attr(length_match, "capture.length")[1] - 1
                ))
            } else {
                length_val <- NA
            }
            chrom_info <- rbind(chrom_info, data.frame(chrom = chrom, length = length_val, stringsAsFactors = FALSE))
        }
    }
    return(chrom_info)
}

# Generate regions based on chunk size
generate_regions <- function(chrom_info, chunk_size) {
    regions <- character()
    for (i in seq_len(nrow(chrom_info))) {
        chrom <- chrom_info$chrom[i]
        if (grepl("_|\\-", chrom)) {
            next
        }
        chrom_length <- chrom_info$length[i]
        if (is.na(chrom_length)) {
            regions <- c(regions, chrom)
        } else {
            start_pos <- 1
            while (start_pos <= chrom_length) {
                end_pos <- min(start_pos + chunk_size - 1, chrom_length)
                regions <- c(regions, paste0(chrom, ":", start_pos, "-", end_pos))
                start_pos <- end_pos + 1
            }
        }
    }
    return(regions)
}

# Check required tools
check_tools <- function() {
    tools <- c("bcftools", "tabix")
    for (tool in tools) {
        if (system(paste("which", tool), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
            stop("Required tool not found: ", tool)
        }
    }
}

# Function to chunk a single VCF file for a specific region
chunk_vcf <- function(vcf_file, region, chunk_dir) {
    base_name <- basename(tools::file_path_sans_ext(tools::file_path_sans_ext(vcf_file)))
    region_name <- gsub("[:-]", "_", region)
    output_file <- file.path(chunk_dir, paste0(base_name, "_", region_name, ".vcf.gz"))
    cmd <- paste(
        "bcftools view -r", shQuote(region), "-f PASS", shQuote(vcf_file),
        "-Oz -o", shQuote(output_file)
    )
    message(cmd)
    res <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (res == 0 && file.exists(output_file) && file.size(output_file) > 0) {
        system(paste("tabix -p vcf", shQuote(output_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)
        return(output_file)
    } else {
        stop(glue("Error: failed to chunk region {region} in file {vcf_file}"))
        return(NULL)
    }
}

# Function to merge chunks for a region
merge_region_chunks <- function(chunk_files, region, merged_dir) {
    if (length(chunk_files) == 0) {
        return(NULL)
    }

    region_name <- gsub("[:-]", "_", region)
    merged_file <- file.path(merged_dir, paste0("merged_", region_name, ".vcf.gz"))

    if (length(chunk_files) > 1) {
        file_list <- glue("{merged_dir}/{region_name}_chunks_file_list.txt")
        # file_list <- tempfile(fileext = ".txt")
        writeLines(chunk_files, file_list)

        temp_merged <- tempfile(fileext = ".vcf.gz")
        cmd1 <- paste("bcftools merge -l", shQuote(file_list), "-Oz -o", shQuote(temp_merged))
        message(cmd1)
        res1 <- system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)

        if (res1 == 0 && file.exists(temp_merged)) {
            system(paste("tabix -p vcf", shQuote(temp_merged)), ignore.stdout = TRUE, ignore.stderr = TRUE)

            format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))
            cmd2 <- paste(
                "bcftools +setGT", shQuote(temp_merged), "-- -t . -n 0",
                "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
                "-Oz -o", shQuote(merged_file)
            )
            message(cmd2)
            res2 <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)

            unlink(temp_merged)
            unlink(paste0(temp_merged, ".tbi"))

            if (res2 != 0 || !file.exists(merged_file) || file.size(merged_file) == 0) {
                message(glue("Error: Normalizing GT for region {region_name} failed."))
                return(NULL)
            }
        } else {
            message(glue("Error: Merging chunks for region {region_name} failed."))
            return(NULL)
        }
    } else {
        single_chunk <- chunk_files[[1]]
        format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))
        cmd_single <- paste(
            "bcftools +setGT", shQuote(single_chunk), "-- -t . -n 0",
            "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
            "-Oz -o", shQuote(merged_file)
        )
        message(cmd_single)
        # Run the command and check for success
        res_single <- system(cmd_single, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (res_single != 0 || !file.exists(merged_file) || file.size(merged_file) == 0) {
            return(NULL)
        }
    }

    system(paste("tabix -p vcf", shQuote(merged_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Post-processing with fill-tags
    temp_fill <- glue("{merged_file}.tmp.vcf.gz")
    # temp_fill <- tempfile(fileext = ".vcf.gz")
    cmd3 <- paste(
        "bcftools annotate -x INFO", shQuote(merged_file),
        "| bcftools +fill-tags -Oz -o", shQuote(temp_fill),
        "-- -t AC,AN,AF,AC_Hom,AC_Het,AC_Hemi,NS"
    )
    message(cmd3)
    res3 <- system(cmd3, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (res3 != 0 || !file.exists(temp_fill) || file.size(temp_fill) == 0) {
        cat("Warning: Post-processing (fill-tags) failed for", merged_file, "\n")
    } else {
        # Overwrite merged_file with the new one
        file.rename(temp_fill, merged_file)
    }

    return(merged_file)
}

# Function to clean up chunk files
cleanup_chunks <- function(chunk_files) {
    for (cf in chunk_files) {
        if (!is.null(cf) && file.exists(cf)) {
            unlink(cf)
            unlink(paste0(cf, ".tbi"))
        }
    }
}

# ===============================
# REVISED parse_region_order():
#   1) numeric contigs (chr1..chr22) → 1..22
#   2) chrX → 23
#   3) chrY → 24
#   4) chrM (or chrMT) → 25
#   5) everything else → 26
# Then sort by (contig_order, start_pos).
# ===============================
parse_region_order <- function(filepath) {
    fname <- basename(filepath)
    parts <- strsplit(gsub("^merged_|\\.vcf\\.gz$", "", fname), "_")[[1]]
    chrom <- parts[1]
    pos <- parts[2]
    # determine contig_order
    if (grepl("^chr([0-9]+)$", chrom)) {
        contig_order <- as.numeric(sub("^chr", "", chrom))
    } else if (chrom == "chrX") {
        contig_order <- 23
    } else if (chrom == "chrY") {
        contig_order <- 24
    } else if (chrom %in% c("chrM", "chrMT")) {
        contig_order <- 25
    } else {
        contig_order <- 26
    }
    # extract start_pos
    start_pos <- suppressWarnings(as.numeric(strsplit(pos, "-")[[1]][1]))
    if (is.na(start_pos)) start_pos <- Inf
    return(c(contig_order, start_pos))
}

merge_vcf_files <- function(vcf_files_to_merge, output_dir = ".") {
    temp_merged <- tempfile(fileext = ".vcf.gz")
    temp_withgt <- tempfile(fileext = ".vcf.gz")
    cmd1 <- paste("bcftools merge -l", shQuote(vcf_files_to_merge), "-Oz -o", shQuote(temp_merged))
    message("\n=== Running commands to merge VCF files and normalize GT:")
    message(cmd1)
    res1 <- system(cmd1, ignore.stdout = TRUE, ignore.stderr = FALSE)

    if (res1 == 0 && file.exists(temp_merged)) {
        system(paste("tabix -p vcf", shQuote(temp_merged)), ignore.stdout = TRUE, ignore.stderr = FALSE)

        cmd2 <- paste(
            "bcftools +setGT", shQuote(temp_merged), "-- -t . -n 0",
            "| bcftools view",
            "-Oz -o", shQuote(temp_withgt)
        )
        message(cmd2)
        res2 <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = FALSE)

        unlink(temp_merged)
        unlink(paste0(temp_merged, ".tbi"))

        if (res2 != 0 || !file.exists(temp_withgt) || file.size(temp_withgt) == 0) {
            stop("Error: Adding GT to merged VCF file failed.")
            return(NULL)
        }
    } else {
        stop("Error: Merging final batch VCF files failed.")
        return(NULL)
    }

    system(paste("tabix -p vcf", shQuote(temp_withgt)), ignore.stdout = TRUE, ignore.stderr = FALSE)

    # Post-processing with fill-tags
    final_merged <- glue("{output_dir}/final_merged.vcf.gz")
    # temp_fill <- tempfile(fileext = ".vcf.gz")
    cmd3 <- paste(
        "bcftools annotate -x INFO", shQuote(temp_withgt),
        "| bcftools +fill-tags -Oz -o", shQuote(final_merged),
        "-- -t AC,AN,AF,AC_Hom,AC_Het,AC_Hemi,NS"
    )
    message(cmd3)
    res3 <- system(cmd3, ignore.stdout = TRUE, ignore.stderr = FALSE)
    if (res3 != 0 || !file.exists(final_merged) || file.size(final_merged) == 0) {
        cat("Warning: Post-processing (fill-tags) failed for", final_merged, "\n")
    }
}
