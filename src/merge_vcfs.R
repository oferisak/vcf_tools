#!/usr/bin/env Rscript

# !/usr/bin/env Rscript


library(parallel)
library(glue)
source("../lib/merge_vcf_helper_functions.R")

process_vcf_files <- function(vcf_input, output_folder, region_chunk_size, join_chunks = FALSE,
                              format_fields = c("GT", "AD"), n_cores = 1, prefix = NULL) {
    # Validate inputs
    if (dir.exists(vcf_input)) {
        vcf_files_raw <- list.files(vcf_input, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
    } else if (file.exists(vcf_input)) {
        vcf_files_raw <- readLines(vcf_input)
    } else {
        stop("Input must be a directory or a text file listing VCF paths: ", vcf_input)
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
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }


    if (length(vcf_files_raw) == 0) {
        stop("No VCF files found in: ", vcf_input)
    }

    cat("Found", length(vcf_files_raw), "VCF file(s)\n")
    cat("FORMAT fields to include:", paste(format_fields, collapse = ", "), "\n")
    cat("Filtering for PASS variants only\n")
    cat("Using", n_cores, "core(s) for parallel steps\n\n")

    # Create output directories
    chunk_dir <- file.path(output_folder, "chunks")
    merged_dir <- file.path(output_folder, "merged_chunks")
    dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
    log_file_name <- ifelse(is.null(prefix), "commands.log", paste0(prefix, "_commands.log"))
    log_file_connection <- file(file.path(output_folder, log_file_name), open = "a")
    sink(log_file_connection, append = TRUE, type = "message")

    # Index VCF files if not already indexed; compress if needed
    cat("Indexing/compressing VCF files...\n")
    vcf_files <- character()
    for (orig_vcf in vcf_files_raw) {
        vcf_path <- orig_vcf
        # If uncompressed, compress now
        if (!grepl("\\.gz$", vcf_path)) {
            gz_path <- paste0(vcf_path, ".gz")
            cmd_bgzip <- paste("bgzip -c", shQuote(vcf_path), ">", shQuote(gz_path))
            message(cmd_bgzip)
            system(cmd_bgzip, ignore.stdout = TRUE, ignore.stderr = TRUE)
            vcf_path <- gz_path
        }
        # Index if not already
        if (!file.exists(paste0(vcf_path, ".tbi")) && !file.exists(paste0(vcf_path, ".csi"))) {
            system(paste("tabix -p vcf", shQuote(vcf_path)), ignore.stdout = TRUE, ignore.stderr = TRUE)
        }
        vcf_files <- c(vcf_files, vcf_path)
    }


    cat("Getting chromosome information...\n")
    chrom_info <- get_chromosomes(vcf_files[1])
    regions <- generate_regions(chrom_info, region_chunk_size)
    cat("Generated", length(regions), "region(s) for processing\n\n")


    #
    # Main processing loop: Process regions in batches
    #
    cat("Processing regions in batches of", n_cores, "...\n")

    # Split regions into batches of n_cores
    region_batches <- split(regions, ceiling(seq_along(regions) / n_cores))
    merged_files <- character(0)

    for (batch_idx in seq_along(region_batches)) {
        current_regions <- region_batches[[batch_idx]]
        cat(
            "Processing batch", batch_idx, "of", length(region_batches),
            "(", length(current_regions), "regions)...\n"
        )

        # Step 1: Chunk VCF files for current batch of regions (parallel)
        cat("  Chunking regions in parallel...\n")
        message("=== Chunking regions ===")
        batch_chunk_results <- list()

        for (region in current_regions) {
            cat("    Chunking region:", region, "...\n")
            region_chunks_raw <- mclapply(
                vcf_files,
                function(v) chunk_vcf(v, region, chunk_dir),
                mc.cores = n_cores
            )
            region_chunks <- Filter(Negate(is.null), region_chunks_raw)
            if (length(region_chunks) > 0) {
                region_chunks <- unlist(region_chunks, use.names = FALSE)
                batch_chunk_results[[region]] <- region_chunks
            }
        }

        cat("  Merging chunks for batch regions in parallel...\n")

        # Wrapper function for parallel merging that includes cleanup
        merge_and_cleanup <- function(region) {
            chunk_files <- batch_chunk_results[[region]]
            cat("    Merging region:", region, "...\n")

            merged_file <- merge_region_chunks(chunk_files, region, merged_dir)
            if (!is.null(merged_file)) {
                cat("    Successfully merged region:", region, "\n")
            } else {
                cat("    Warning: Failed to merge region:", region, "\n")
            }

            # Clean up chunk files immediately after merging
            cleanup_chunks(chunk_files)

            return(merged_file)
        }
        message("=== Merging regions ===")
        # Run merging in parallel across regions in the current batch
        batch_merged_results <- mclapply(
            names(batch_chunk_results),
            merge_and_cleanup,
            mc.cores = n_cores
        )

        # Collect successful merges
        batch_merged_files <- Filter(Negate(is.null), batch_merged_results)
        if (length(batch_merged_files) > 0) {
            merged_files <- c(merged_files, unlist(batch_merged_files, use.names = FALSE))
        }

        cat("  Completed batch", batch_idx, "\n\n")
    }

    cat("Processing complete. Created", length(merged_files), "merged region file(s).\n")

    # Clean up chunk directory if empty
    if (dir.exists(chunk_dir) && length(list.files(chunk_dir, recursive = TRUE)) == 0) {
        unlink(chunk_dir, recursive = TRUE)
        cat("Cleaned up chunk directory\n")
    }
    #
    # Step 3: Join all merged chunks if requested (with improved sorting)
    #
    if (join_chunks && length(merged_files) > 1) {
        cat("Step 3: Joining all merged chunks...\n")
        message("=== Joining merged chunks ===")
        # Print out merged_files and their sizes
        cat("Merged files to concatenate (", length(merged_files), "):\n", sep = "")
        for (f in merged_files) {
            size_info <- if (file.exists(f)) file.info(f)$size else NA
            cat("  ", f, " (bytes:", size_info, ")\n", sep = "")
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
        final_output_name <- ifelse(is.null(prefix), "final_merged_all_samples_all_regions.vcf.gz",
            paste0(prefix, "_final_merged_all_samples_all_regions.vcf.gz")
        )
        final_output <- file.path(output_folder, final_output_name)
        cat("Final output path:", final_output, "\n")

        # Build bcftools concat command
        message("concatenating files with bcftools concat...")
        cmd_concat <- paste("bcftools concat -f", shQuote(file_list), "-Oz -o", shQuote(final_output))
        message(cmd_concat)
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
            cat("Final merged file exists; size (bytes):", final_size, ", Indexing..", "\n\n")
            system(paste("tabix -p vcf", shQuote(final_output)), ignore.stdout = TRUE, ignore.stderr = TRUE)
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
        final_file <- file.path(output_folder, "final_merged_all_samples_all_regions.vcf.gz")
        if (file.exists(final_file)) {
            cat("Final merged file:", final_file, "\n")
        }
    }
    # Close the log file
    sink()
    return(invisible(TRUE))
}

# Command line interface
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        cat("Usage: Rscript merge_vcfs.R <vcf_input> <region_chunk_size> [join_chunks] [format_fields] [n_cores]\n")
        cat("  vcf_input: Path to folder / text file containing VCF files\n")
        cat("  output_folder: Path to folder containing VCF files\n")
        cat("  region_chunk_size: Numeric value for chunk size (bp)\n")
        cat("  join_chunks: TRUE/FALSE (optional, default: FALSE)\n")
        cat("  format_fields: Comma-separated FORMAT fields (optional, default: GT,AD)\n")
        cat("  n_cores: Number of cores for parallel processing (optional, default: 1)\n")
        quit(status = 1)
    }

    vcf_input <- args[1]
    output_folder <- args[2]
    region_chunk_size <- as.numeric(args[3])
    join_chunks <- if (length(args) >= 4) as.logical(args[4]) else FALSE
    format_fields <- if (length(args) >= 5 && nzchar(args[5])) strsplit(args[5], ",")[[1]] else c("GT", "AD")
    n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1
    prefix <- if (length(args) >= 7) args[7] else NULL

    process_vcf_files(vcf_input, output_folder, region_chunk_size, join_chunks, format_fields, n_cores, prefix)
}
