#!/usr/bin/env Rscript

# VCF Processing Script (Parallelized, with unique temp files)
# Processes VCF files by chunking, merging, and optionally joining
# Requires: bcftools, tabix

library(parallel)

process_vcf_files <- function(vcf_folder, region_chunk_size, join_chunks = FALSE,
                              format_fields = c("GT", "AD"), n_cores = 1) {
    # Validate inputs
    if (!dir.exists(vcf_folder)) {
        stop("VCF folder does not exist: ", vcf_folder)
    }
    if (!is.numeric(region_chunk_size) || region_chunk_size <= 0) {
        stop("region_chunk_size must be a positive numeric value")
    }
    if (!is.logical(join_chunks)) {
        stop("join_chunks must be TRUE or FALSE")
    }
    if (!is.character(format_fields) || length(format_fields) == 0) {
        stop("format_fields must be a character vector with at least one field")
    }
    if (!is.numeric(n_cores) || n_cores < 1) {
        stop("n_cores must be a positive integer")
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
    check_tools()

    # Get VCF files (compressed or not)
    vcf_files_raw <- list.files(vcf_folder, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
    if (length(vcf_files_raw) == 0) {
        stop("No VCF files found in: ", vcf_folder)
    }

    cat("Found", length(vcf_files_raw), "VCF file(s)\n")
    cat("FORMAT fields to include:", paste(format_fields, collapse = ", "), "\n")
    cat("Filtering for PASS variants only\n")
    cat("Using", n_cores, "core(s) for parallel steps\n\n")

    # Create output directories
    chunk_dir <- file.path(vcf_folder, "chunks")
    merged_dir <- file.path(vcf_folder, "merged_chunks")
    dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

    # Index VCF files if not already indexed; compress if needed
    cat("Indexing/compressing VCF files...\n")
    vcf_files <- character()
    for (orig_vcf in vcf_files_raw) {
        vcf_path <- orig_vcf
        # If uncompressed, compress now
        if (!grepl("\\.gz$", vcf_path)) {
            gz_path <- paste0(vcf_path, ".gz")
            cmd_bgzip <- paste("bgzip -c", shQuote(vcf_path), ">", shQuote(gz_path))
            system(cmd_bgzip, ignore.stdout = TRUE, ignore.stderr = TRUE)
            vcf_path <- gz_path
        }
        # Index if not already
        if (!file.exists(paste0(vcf_path, ".tbi")) && !file.exists(paste0(vcf_path, ".csi"))) {
            system(paste("tabix -p vcf", shQuote(vcf_path)), ignore.stdout = TRUE, ignore.stderr = TRUE)
        }
        vcf_files <- c(vcf_files, vcf_path)
    }

    # Get chromosomes and their lengths from the first VCF file
    get_chromosomes <- function(vcf_file) {
        cmd <- paste("bcftools view -h", shQuote(vcf_file), "| grep '^##contig'")
        contig_lines <- system(cmd, intern = TRUE)
        if (length(contig_lines) == 0) {
            # Fallback: query unique chromosomes in data
            cmd2 <- paste("bcftools query -f '%CHROM\\n'", shQuote(vcf_file), "| sort -u")
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

    cat("Getting chromosome information...\n")
    chrom_info <- get_chromosomes(vcf_files[1])
    regions <- generate_regions(chrom_info, region_chunk_size)
    cat("Generated", length(regions), "region(s) for processing\n\n")

    #
    # Step 1: Chunk each VCF file (in parallel per region)
    #
    cat("Step 1: Chunking VCF files (parallel per region)...\n")
    chunk_vcf <- function(vcf_file, region) {
        base_name <- basename(tools::file_path_sans_ext(tools::file_path_sans_ext(vcf_file)))
        region_name <- gsub("[:-]", "_", region)
        output_file <- file.path(chunk_dir, paste0(base_name, "_", region_name, ".vcf.gz"))
        cmd <- paste(
            "bcftools view -r", shQuote(region), "-f PASS", shQuote(vcf_file),
            "-Oz -o", shQuote(output_file)
        )
        res <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (res == 0 && file.exists(output_file) && file.size(output_file) > 0) {
            system(paste("tabix -p vcf", shQuote(output_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)
            return(output_file)
        }
        return(NULL)
    }

    chunk_results <- list()
    for (region in regions) {
        cat("  Chunking region:", region, "...\n")
        region_chunks_raw <- mclapply(
            vcf_files,
            function(v) chunk_vcf(v, region),
            mc.cores = n_cores
        )
        region_chunks <- Filter(Negate(is.null), region_chunks_raw)
        if (length(region_chunks) > 0) {
            region_chunks <- unlist(region_chunks, use.names = FALSE)
            chunk_results[[region]] <- region_chunks
        }
    }
    cat("Completed chunking.\n\n")

    #
    # Step 2: Merge chunks for each region (parallel over regions)
    #
    cat("Step 2: Merging chunks by region (parallel)...\n")
    merge_region <- function(region) {
        chunk_files <- chunk_results[[region]]
        if (length(chunk_files) == 0) {
            return(NULL)
        }

        region_name <- gsub("[:-]", "_", region)
        merged_file <- file.path(merged_dir, paste0("merged_", region_name, ".vcf.gz"))

        if (length(chunk_files) > 1) {
            file_list <- tempfile(fileext = ".txt")
            writeLines(chunk_files, file_list)

            temp_merged <- tempfile(fileext = ".vcf.gz")
            cmd1 <- paste("bcftools merge -l", shQuote(file_list), "-Oz -o", shQuote(temp_merged))
            res1 <- system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)

            if (res1 == 0 && file.exists(temp_merged)) {
                system(paste("tabix -p vcf", shQuote(temp_merged)), ignore.stdout = TRUE, ignore.stderr = TRUE)

                format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))
                cmd2 <- paste(
                    "bcftools +setGT", shQuote(temp_merged), "-- -t . -n 0",
                    "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
                    "-Oz -o", shQuote(merged_file)
                )
                res2 <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)

                unlink(temp_merged)
                unlink(paste0(temp_merged, ".tbi"))

                if (res2 != 0 || !file.exists(merged_file) || file.size(merged_file) == 0) {
                    unlink(file_list)
                    return(NULL)
                }
            } else {
                unlink(file_list)
                return(NULL)
            }

            unlink(file_list)
        } else {
            single_chunk <- chunk_files[[1]]
            format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))
            cmd_single <- paste(
                "bcftools +setGT", shQuote(single_chunk), "-- -t . -n 0",
                "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
                "-Oz -o", shQuote(merged_file)
            )
            res_single <- system(cmd_single, ignore.stdout = TRUE, ignore.stderr = TRUE)
            if (res_single != 0 || !file.exists(merged_file) || file.size(merged_file) == 0) {
                return(NULL)
            }
        }

        system(paste("tabix -p vcf", shQuote(merged_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)

        # === BEGIN FIXED: use a unique tempfile for fill-tags ===
        temp_fill <- tempfile(fileext = ".vcf.gz")
        cmd3 <- paste(
            "bcftools annotate -x INFO", shQuote(merged_file),
            "| bcftools +fill-tags -Oz -o", shQuote(temp_fill),
            "-- -t AC,AN,AF,AC_Hom,AC_Het,AC_Hemi,NS"
        )
        res3 <- system(cmd3, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (res3 != 0 || !file.exists(temp_fill) || file.size(temp_fill) == 0) {
            cat("Warning: Post-processing (fill-tags) failed for", merged_file, "\n")
        } else {
            # Overwrite merged_file with the new one
            file.rename(temp_fill, merged_file)
        }
        # ============================================================

        for (cf in chunk_files) {
            unlink(cf)
            unlink(paste0(cf, ".tbi"))
        }

        return(merged_file)
    }

    merged_results <- mclapply(names(chunk_results), merge_region, mc.cores = n_cores)
    merged_files <- unlist(merged_results, use.names = FALSE)
    cat("Completed merging. Created", length(merged_files), "merged region file(s).\n\n")

    # Clean up chunk directory if empty
    if (length(list.files(chunk_dir, recursive = TRUE)) == 0) {
        unlink(chunk_dir, recursive = TRUE)
        cat("Cleaned up chunk directory\n\n")
    }

    #
    # Step 3: Join all merged chunks if requested (with improved sorting)
    #
    if (join_chunks && length(merged_files) > 1) {
        cat("Step 3: Joining all merged chunks...\n")

        # Print out merged_files and their sizes
        cat("Merged files to concatenate (", length(merged_files), "):\n", sep = "")
        for (f in merged_files) {
            size_info <- if (file.exists(f)) file.info(f)$size else NA
            cat("  ", f, " (bytes:", size_info, ")\n", sep = "")
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

        order_mat <- t(sapply(merged_files, parse_region_order))
        ord_idx <- order(order_mat[, 1], order_mat[, 2])
        sorted_files <- merged_files[ord_idx]

        cat("Sorted files order:\n")
        print(sorted_files)

        # Write sorted file list to disk
        file_list <- tempfile(fileext = ".txt")
        cat("Writing file list to:", file_list, "\n")
        writeLines(sorted_files, file_list)

        # Show first 10 lines of file_list
        cat("First 10 lines of file_list:\n")
        flines <- readLines(file_list, n = 10)
        print(flines)
        if (length(flines) < length(sorted_files)) {
            cat("... (", length(sorted_files) - length(flines), "more lines) ...\n", sep = "")
        }

        final_output <- file.path(vcf_folder, "final_merged_all_samples_all_regions.vcf.gz")
        cat("Final output path:", final_output, "\n")

        # Build bcftools concat command
        cmd_concat <- paste("bcftools concat -f", shQuote(file_list), "-Oz -o", shQuote(final_output))
        cat("Running command:\n  ", cmd_concat, "\n", sep = "")

        # Use system2 to capture stdout/stderr
        concat_output <- system2("bcftools",
            args = c("concat", "-f", file_list, "-Oz", "-o", final_output),
            stdout = TRUE,
            stderr = TRUE
        )
        exit_status <- attr(concat_output, "status")
        if (is.null(exit_status)) exit_status <- 0

        cat("=== bcftools concat stdout+stderr begin ===\n")
        cat(paste(concat_output, collapse = "\n"), "\n")
        cat("=== bcftools concat stdout+stderr end ===\n")
        cat("Exit status:", exit_status, "\n")

        # Check if final_output exists and print its size
        if (file.exists(final_output)) {
            final_size <- file.info(final_output)$size
            cat("Final merged file exists; size (bytes):", final_size, "\n\n")
        } else {
            cat("Error: Final merged file does not exist after concat.\n\n")
        }

        unlink(file_list)
    } else if (join_chunks) {
        cat("Only one merged region file exists; skipping concatenation.\n\n")
    }

    cat("Processing complete!\n")
    cat("Merged region files directory:", merged_dir, "\n")
    if (join_chunks) {
        final_file <- file.path(vcf_folder, "final_merged_all_samples_all_regions.vcf.gz")
        if (file.exists(final_file)) {
            cat("Final merged file:", final_file, "\n")
        }
    }

    return(invisible(TRUE))
}

# Command line interface
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        cat("Usage: Rscript merge_vcfs.R <vcf_folder> <region_chunk_size> [join_chunks] [format_fields] [n_cores]\n")
        cat("  vcf_folder: Path to folder containing VCF files\n")
        cat("  region_chunk_size: Numeric value for chunk size (bp)\n")
        cat("  join_chunks: TRUE/FALSE (optional, default: FALSE)\n")
        cat("  format_fields: Comma-separated FORMAT fields (optional, default: GT,AD)\n")
        cat("  n_cores: Number of cores for parallel processing (optional, default: 1)\n")
        quit(status = 1)
    }

    vcf_folder <- args[1]
    region_chunk_size <- as.numeric(args[2])
    join_chunks <- if (length(args) >= 3) as.logical(args[3]) else FALSE
    format_fields <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else c("GT", "AD")
    n_cores <- if (length(args) >= 5) as.integer(args[5]) else 1

    process_vcf_files(vcf_folder, region_chunk_size, join_chunks, format_fields, n_cores)
}
