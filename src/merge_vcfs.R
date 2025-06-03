#!/usr/bin/env Rscript

# VCF Processing Script
# Processes VCF files by chunking, merging, and optionally joining
# Requires: bcftools, tabix

library(parallel)

process_vcf_files <- function(vcf_folder, region_chunk_size, join_chunks = FALSE, format_fields = c("GT", "AD")) {
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

    # Get VCF files
    vcf_files <- list.files(vcf_folder, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)

    if (length(vcf_files) == 0) {
        stop("No VCF files found in: ", vcf_folder)
    }

    cat("Found", length(vcf_files), "VCF files\n")
    cat("FORMAT fields to include:", paste(format_fields, collapse = ", "), "\n")
    cat("Filtering for PASS variants only\n")

    # Create output directories
    chunk_dir <- file.path(vcf_folder, "chunks")
    merged_dir <- file.path(vcf_folder, "merged_chunks")

    dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

    # Index VCF files if not already indexed
    cat("Indexing VCF files...\n")
    for (vcf_file in vcf_files) {
        if (!file.exists(paste0(vcf_file, ".tbi")) && !file.exists(paste0(vcf_file, ".csi"))) {
            if (!grepl("\\.gz$", vcf_file)) {
                # Compress and index
                system(paste("bgzip -c", shQuote(vcf_file), ">", shQuote(paste0(vcf_file, ".gz"))))
                vcf_file <- paste0(vcf_file, ".gz")
            }
            system(paste("tabix -p vcf", shQuote(vcf_file)))
        }
    }

    # Update vcf_files list to include compressed versions
    vcf_files <- list.files(vcf_folder, pattern = "\\.vcf\\.gz$", full.names = TRUE)

    # Get chromosomes and their lengths from the first VCF file
    get_chromosomes <- function(vcf_file) {
        cmd <- paste("bcftools view -h", shQuote(vcf_file), "| grep '^##contig'")
        contig_lines <- system(cmd, intern = TRUE)

        if (length(contig_lines) == 0) {
            # Fallback: get chromosomes from data
            cmd <- paste("bcftools query -f '%CHROM\\n'", shQuote(vcf_file), "| sort -u")
            chroms <- system(cmd, intern = TRUE)
            return(data.frame(chrom = chroms, length = NA, stringsAsFactors = FALSE))
        }

        # Parse contig lines
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

        for (i in 1:nrow(chrom_info)) {
            chrom <- chrom_info$chrom[i]
            chrom_length <- chrom_info$length[i]

            if (is.na(chrom_length)) {
                # If length is unknown, create single region for the chromosome
                regions <- c(regions, chrom)
            } else {
                # Create chunks
                start_pos <- 1
                chunk_num <- 1

                while (start_pos <= chrom_length) {
                    end_pos <- min(start_pos + chunk_size - 1, chrom_length)
                    region <- paste0(chrom, ":", start_pos, "-", end_pos)
                    regions <- c(regions, region)
                    start_pos <- end_pos + 1
                    chunk_num <- chunk_num + 1
                }
            }
        }

        return(regions)
    }

    cat("Getting chromosome information...\n")
    chrom_info <- get_chromosomes(vcf_files[1])
    regions <- generate_regions(chrom_info, region_chunk_size)

    cat("Generated", length(regions), "regions for processing\n")

    # Step 1: Chunk each VCF file
    cat("Step 1: Chunking VCF files...\n")

    chunk_vcf <- function(vcf_file, region) {
        base_name <- basename(tools::file_path_sans_ext(tools::file_path_sans_ext(vcf_file)))
        region_name <- gsub("[:-]", "_", region)
        output_file <- file.path(chunk_dir, paste0(base_name, "_", region_name, ".vcf.gz"))

        cmd <- paste(
            "bcftools view -r", shQuote(region), "-f PASS", shQuote(vcf_file),
            "-Oz -o", shQuote(output_file)
        )

        result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

        if (result == 0 && file.exists(output_file) && file.size(output_file) > 0) {
            # Index the chunk
            system(paste("tabix -p vcf", shQuote(output_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)
            return(output_file)
        } else {
            return(NULL)
        }
    }

    # Process all combinations of VCF files and regions
    total_chunks <- length(vcf_files) * length(regions)
    cat("Processing", total_chunks, "chunks...\n")

    chunk_results <- list()
    pb <- txtProgressBar(min = 0, max = total_chunks, style = 3)
    counter <- 0

    for (region in regions) {
        region_chunks <- character()

        for (vcf_file in vcf_files) {
            chunk_file <- chunk_vcf(vcf_file, region)
            if (!is.null(chunk_file)) {
                region_chunks <- c(region_chunks, chunk_file)
            }
            counter <- counter + 1
            setTxtProgressBar(pb, counter)
        }

        if (length(region_chunks) > 0) {
            chunk_results[[region]] <- region_chunks
        }
    }

    close(pb)

    # Step 2: Merge chunks for each region
    cat("\nStep 2: Merging chunks by region...\n")

    merged_files <- character()
    pb <- txtProgressBar(min = 0, max = length(chunk_results), style = 3)

    for (i in seq_along(chunk_results)) {
        region <- names(chunk_results)[i]
        chunk_files <- chunk_results[[i]]

        if (length(chunk_files) > 1) {
            region_name <- gsub("[:-]", "_", region)
            merged_file <- file.path(merged_dir, paste0("merged_", region_name, ".vcf.gz"))

            # Create file list for bcftools merge
            file_list <- tempfile(fileext = ".txt")
            writeLines(chunk_files, file_list)

            # First merge the files normally
            temp_merged <- tempfile(fileext = ".vcf.gz")
            cmd1 <- paste("bcftools merge -l", shQuote(file_list), "-Oz -o", shQuote(temp_merged))
            result1 <- system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)

            if (result1 == 0 && file.exists(temp_merged)) {
                # Index the temporary file
                system(paste("tabix -p vcf", shQuote(temp_merged)), ignore.stdout = TRUE, ignore.stderr = TRUE)

                # Set missing genotypes to reference and filter FORMAT fields (only PASS variants)
                format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))

                cmd2 <- paste(
                    "bcftools +setGT", shQuote(temp_merged), "-- -t . -n 0",
                    "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
                    "-Oz -o", shQuote(merged_file)
                )

                result <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)

                result <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)

                # Clean up temp file
                unlink(temp_merged)
                unlink(paste0(temp_merged, ".tbi"))
            } else {
                result <- 1
            }

            if (result == 0 && file.exists(merged_file) && file.size(merged_file) > 0) {
                system(paste("tabix -p vcf", shQuote(merged_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)
                merged_files <- c(merged_files, merged_file)

                # Clean up chunk files for this region
                for (chunk_file in chunk_files) {
                    unlink(chunk_file)
                    unlink(paste0(chunk_file, ".tbi"))
                }
            }

            unlink(file_list)
        } else if (length(chunk_files) == 1) {
            # Single file, process format fields and missing genotypes
            region_name <- gsub("[:-]", "_", region)
            merged_file <- file.path(merged_dir, paste0("merged_", region_name, ".vcf.gz"))

            format_keep <- paste0("FORMAT/", paste(c(format_fields, "FT"), collapse = ",FORMAT/"))

            # Process single file with missing genotype handling and format filtering (only PASS variants)
            format_keep <- paste0("FORMAT/", paste(format_fields, collapse = ",FORMAT/"))

            cmd <- paste(
                "bcftools +setGT", shQuote(chunk_files[1]), "-- -t . -n 0",
                "| bcftools annotate --remove", shQuote(paste0("^", format_keep)),
                "-Oz -o", shQuote(merged_file)
            )

            result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

            if (result == 0 && file.exists(merged_file) && file.size(merged_file) > 0) {
                system(paste("tabix -p vcf", shQuote(merged_file)), ignore.stdout = TRUE, ignore.stderr = TRUE)
                merged_files <- c(merged_files, merged_file)

                # Clean up the single chunk file
                unlink(chunk_files[1])
                unlink(paste0(chunk_files[1], ".tbi"))
            }

            # No temp files to clean up for single file case
        }

        setTxtProgressBar(pb, i)
    }

    close(pb)

    cat("\nCreated", length(merged_files), "merged region files\n")

    # Clean up empty chunk directory if all files were removed
    if (length(list.files(chunk_dir)) == 0) {
        unlink(chunk_dir, recursive = TRUE)
        cat("Cleaned up chunk files\n")
    }

    # Step 3: Join all chunks if requested
    if (join_chunks && length(merged_files) > 1) {
        cat("Step 3: Joining all merged chunks...\n")

        final_output <- file.path(vcf_folder, "final_merged_all_samples_all_regions.vcf.gz")

        # Sort merged files by genomic position
        sorted_files <- sort(merged_files)

        # Create file list
        file_list <- tempfile(fileext = ".txt")
        writeLines(sorted_files, file_list)

        cmd <- paste("bcftools concat -f", shQuote(file_list), "-Oz -o", shQuote(final_output))
        result <- system(cmd)

        if (result == 0 && file.exists(final_output)) {
            system(paste("tabix -p vcf", shQuote(final_output)))
            cat("Final merged file created:", final_output, "\n")
        } else {
            cat("Error creating final merged file\n")
        }

        unlink(file_list)
    } else if (join_chunks) {
        cat("Step 3: Only one merged file exists, no joining needed\n")
    }

    # Summary
    cat("\nProcessing complete!\n")
    cat("Chunk files:", file.path(chunk_dir, "*"), "\n")
    cat("Merged region files:", file.path(merged_dir, "*"), "\n")

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
        cat("Usage: Rscript vcf_processing.R <vcf_folder> <region_chunk_size> [join_chunks] [format_fields]\n")
        cat("  vcf_folder: Path to folder containing VCF files\n")
        cat("  region_chunk_size: Numeric value for chunk size (bp)\n")
        cat("  join_chunks: TRUE/FALSE (optional, default: FALSE)\n")
        cat("  format_fields: Comma-separated FORMAT fields (optional, default: GT,AD)\n")
        quit(status = 1)
    }

    vcf_folder <- args[1]
    region_chunk_size <- as.numeric(args[2])
    join_chunks <- if (length(args) >= 3) as.logical(args[3]) else FALSE
    format_fields <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("GT", "AD")

    process_vcf_files(vcf_folder, region_chunk_size, join_chunks, format_fields)
}