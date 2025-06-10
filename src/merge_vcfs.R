#!/usr/bin/env Rscript

# !/usr/bin/env Rscript

library(parallel)
library(glue)
source("../lib/merge_vcf_helper_functions.R")

process_vcf_files <- function(vcf_input, output_folder, region_chunk_size, join_chunks = FALSE,
                              format_fields = c("GT", "AD"), n_cores = 1, vcf_batch_size = NULL) {
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
    if (!is.null(vcf_batch_size) && (!is.numeric(vcf_batch_size) || vcf_batch_size < 1)) {
        stop("vcf_batch_size must be NULL or a positive integer")
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
    cat("Using", n_cores, "core(s) for parallel steps\n")

    # VCF batching logic
    if (is.null(vcf_batch_size) || vcf_batch_size >= length(vcf_files_raw)) {
        cat("Processing all VCFs together (no batching)\n\n")
        vcf_batch_size <- length(vcf_files_raw)
        vcf_batches <- list(vcf_files_raw)
    } else {
        cat("VCF batch size:", vcf_batch_size, "\n")
        vcf_batches <- split(vcf_files_raw, ceiling(seq_along(vcf_files_raw) / vcf_batch_size))
        cat("Created", length(vcf_batches), "VCF batch(es)\n\n")
    }

    # Store paths to batch final VCFs for ultimate merging
    batch_final_vcfs <- character(0)

    # Process each VCF batch
    for (batch_idx in seq_along(vcf_batches)) {
        current_vcf_batch <- vcf_batches[[batch_idx]]
        cat(
            "=== Processing VCF Batch", batch_idx, "of", length(vcf_batches),
            "(", length(current_vcf_batch), "VCFs) ===\n"
        )

        # Create batch-specific output directory
        batch_output_folder <- file.path(output_folder, paste0("batch_", batch_idx))
        if (!dir.exists(batch_output_folder)) {
            dir.create(batch_output_folder, recursive = TRUE)
        }

        # Process current batch of VCFs
        batch_final_vcf <- process_single_vcf_batch(
            current_vcf_batch, batch_output_folder, region_chunk_size,
            join_chunks = TRUE, format_fields, n_cores, batch_idx
        )

        if (!is.null(batch_final_vcf) && file.exists(batch_final_vcf)) {
            batch_final_vcfs <- c(batch_final_vcfs, batch_final_vcf)
            cat("Batch", batch_idx, "completed successfully. Final VCF:", batch_final_vcf, "\n\n")
        } else {
            cat("Warning: Batch", batch_idx, "failed to produce final VCF\n\n")
        }
    }

    # Final step: Merge all batch final VCFs if we have multiple batches
    if (length(vcf_batches) > 1 && length(batch_final_vcfs) > 1) {
        cat("=== Merging all batch final VCFs into ultimate final VCF ===\n")
        ultimate_final_vcf <- merge_batch_final_vcfs(batch_final_vcfs, output_folder)

        if (!is.null(ultimate_final_vcf) && file.exists(ultimate_final_vcf)) {
            cat("Ultimate final VCF created:", ultimate_final_vcf, "\n")
        } else {
            cat("Warning: Failed to create ultimate final VCF\n")
        }
    } else if (length(batch_final_vcfs) == 1) {
        # Copy single batch final VCF to main output folder
        ultimate_final_vcf <- file.path(output_folder, "final_merged_all_samples_all_regions.vcf.gz")
        file.copy(batch_final_vcfs[1], ultimate_final_vcf, overwrite = TRUE)
        cat("Single batch processed. Final VCF:", ultimate_final_vcf, "\n")
    }

    cat("Processing complete!\n")
    cat("Output directory:", output_folder, "\n")
    if (exists("ultimate_final_vcf") && file.exists(ultimate_final_vcf)) {
        cat("Ultimate final merged file:", ultimate_final_vcf, "\n")
    }

    return(invisible(TRUE))
}

# Process a single batch of VCF files (extracted from original logic)
process_single_vcf_batch <- function(vcf_files_raw, output_folder, region_chunk_size,
                                     join_chunks = TRUE, format_fields, n_cores, batch_idx) {
    # Create output directories
    chunk_dir <- file.path(output_folder, "chunks")
    merged_dir <- file.path(output_folder, "merged_chunks")
    dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

    # Index VCF files if not already indexed; compress if needed
    cat("  Indexing/compressing VCF files for batch", batch_idx, "...\n")
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

    cat("  Getting chromosome information for batch", batch_idx, "...\n")
    chrom_info <- get_chromosomes(vcf_files[1])
    regions <- generate_regions(chrom_info, region_chunk_size)
    cat("  Generated", length(regions), "region(s) for batch", batch_idx, "\n")

    # Main processing loop: Process regions in batches
    cat("  Processing regions in batches of", n_cores, "for batch", batch_idx, "...\n")

    # Split regions into batches of n_cores
    region_batches <- split(regions, ceiling(seq_along(regions) / n_cores))
    merged_files <- character(0)

    for (region_batch_idx in seq_along(region_batches)) {
        current_regions <- region_batches[[region_batch_idx]]
        cat(
            "    Processing region batch", region_batch_idx, "of", length(region_batches),
            "(", length(current_regions), "regions) for VCF batch", batch_idx, "...\n"
        )

        # Step 1: Chunk VCF files for current batch of regions (parallel)
        cat("      Chunking regions in parallel...\n")
        batch_chunk_results <- list()

        for (region in current_regions) {
            cat("        Chunking region:", region, "...\n")
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

        cat("      Merging chunks for batch regions in parallel...\n")

        # Wrapper function for parallel merging that includes cleanup
        merge_and_cleanup <- function(region) {
            chunk_files <- batch_chunk_results[[region]]
            cat("        Merging region:", region, "...\n")

            merged_file <- merge_region_chunks(chunk_files, region, merged_dir)
            if (!is.null(merged_file)) {
                cat("        Successfully merged region:", region, "\n")
            } else {
                cat("        Warning: Failed to merge region:", region, "\n")
            }

            # Clean up chunk files immediately after merging
            cleanup_chunks(chunk_files)

            return(merged_file)
        }

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

        cat("      Completed region batch", region_batch_idx, "for VCF batch", batch_idx, "\n")
    }

    cat("  Region processing complete for batch", batch_idx, ". Created", length(merged_files), "merged region file(s).\n")

    # Clean up chunk directory if empty
    if (dir.exists(chunk_dir) && length(list.files(chunk_dir, recursive = TRUE)) == 0) {
        unlink(chunk_dir, recursive = TRUE)
        cat("  Cleaned up chunk directory for batch", batch_idx, "\n")
    }

    # Join all merged chunks for this batch
    if (join_chunks && length(merged_files) > 1) {
        cat("  Joining all merged chunks for batch", batch_idx, "...\n")

        # Print out merged_files and their sizes
        cat("  Merged files to concatenate (", length(merged_files), "):\n", sep = "")
        for (f in merged_files) {
            size_info <- if (file.exists(f)) file.info(f)$size else NA
            cat("    ", f, " (bytes:", size_info, ")\n", sep = "")
        }

        order_mat <- t(sapply(merged_files, parse_region_order))
        ord_idx <- order(order_mat[, 1], order_mat[, 2])
        sorted_files <- merged_files[ord_idx]

        cat("  Sorted files order:\n")
        print(sorted_files)

        # Write sorted file list to disk
        file_list <- tempfile(fileext = ".txt")
        cat("  Writing file list to:", file_list, "\n")
        writeLines(sorted_files, file_list)

        final_output <- file.path(output_folder, paste0("final_merged_batch_", batch_idx, "_all_samples_all_regions.vcf.gz"))
        cat("  Final output path for batch", batch_idx, ":", final_output, "\n")

        # Build bcftools concat command
        cmd_concat <- paste("bcftools concat -f", shQuote(file_list), "-Oz -o", shQuote(final_output))
        cat("  Running command:\n    ", cmd_concat, "\n", sep = "")

        # Use system2 to capture stdout/stderr
        concat_output <- system2("bcftools",
            args = c("concat", "-f", file_list, "-Oz", "-o", final_output),
            stdout = TRUE,
            stderr = TRUE
        )
        exit_status <- attr(concat_output, "status")
        if (is.null(exit_status)) exit_status <- 0

        cat("  === bcftools concat stdout+stderr begin ===\n")
        cat("  ", paste(concat_output, collapse = "\n  "), "\n", sep = "")
        cat("  === bcftools concat stdout+stderr end ===\n")
        cat("  Exit status:", exit_status, "\n")

        # Check if final_output exists and print its size
        if (file.exists(final_output)) {
            final_size <- file.info(final_output)$size
            cat("  Final merged file for batch", batch_idx, "exists; size (bytes):", final_size, "\n")
            # Index the final output
            system(paste("tabix -p vcf", shQuote(final_output)), ignore.stdout = TRUE, ignore.stderr = TRUE)
            cat("  Indexed final output for batch", batch_idx, "\n")
            unlink(file_list)
            return(final_output)
        } else {
            cat("  Error: Final merged file for batch", batch_idx, "does not exist after concat.\n")
            unlink(file_list)
            return(NULL)
        }
    } else if (join_chunks && length(merged_files) == 1) {
        cat("  Only one merged region file exists for batch", batch_idx, "; copying as final output.\n")
        final_output <- file.path(output_folder, paste0("final_merged_batch_", batch_idx, "_all_samples_all_regions.vcf.gz"))
        file.copy(merged_files[1], final_output, overwrite = TRUE)
        # Index the final output file
        system(paste("tabix -p vcf", shQuote(final_output)), ignore.stdout = TRUE, ignore.stderr = TRUE)
        cat("  Indexed final output for batch", batch_idx, "\n")
        return(final_output)
    } else {
        cat("  join_chunks is FALSE, skipping concatenation for batch", batch_idx, ".\n")
        return(NULL)
    }
}

# Merge all batch final VCFs into ultimate final VCF
merge_batch_final_vcfs <- function(batch_final_vcfs, output_folder) {
    if (length(batch_final_vcfs) <= 1) {
        return(batch_final_vcfs[1])
    }

    cat("Merging", length(batch_final_vcfs), "batch final VCFs...\n")

    # Print batch final VCFs and their sizes
    for (i in seq_along(batch_final_vcfs)) {
        f <- batch_final_vcfs[i]
        size_info <- if (file.exists(f)) file.info(f)$size else NA
        cat("  Batch", i, "final VCF:", f, "(bytes:", size_info, ")\n")
    }

    # Create file list for merging
    file_list <- tempfile(fileext = ".txt")
    writeLines(batch_final_vcfs, file_list)

    ultimate_final_output <- file.path(output_folder, "final_merged_all_samples_all_regions.vcf.gz")
    cat("Ultimate final output path:", ultimate_final_output, "\n")

    # First merge all batch VCFs using bcftools merge
    temp_merged <- tempfile(fileext = ".vcf.gz")
    cmd_merge <- paste("bcftools merge -l", shQuote(file_list), "-Oz -o", shQuote(temp_merged))
    cat("Running merge command:\n  ", cmd_merge, "\n", sep = "")

    merge_output <- system2("bcftools",
        args = c("merge", "-l", file_list, "-Oz", "-o", temp_merged),
        stdout = TRUE,
        stderr = TRUE
    )
    merge_exit_status <- attr(merge_output, "status")
    if (is.null(merge_exit_status)) merge_exit_status <- 0

    cat("=== bcftools merge stdout+stderr begin ===\n")
    cat(paste(merge_output, collapse = "\n"), "\n")
    cat("=== bcftools merge stdout+stderr end ===\n")
    cat("Merge exit status:", merge_exit_status, "\n")

    if (merge_exit_status != 0 || !file.exists(temp_merged)) {
        cat("Error: Failed to merge batch final VCFs\n")
        unlink(file_list)
        return(NULL)
    }

    # Index the temporary merged file
    system(paste("tabix -p vcf", shQuote(temp_merged)), ignore.stdout = TRUE, ignore.stderr = TRUE)


    temp_setgt <- tempfile(fileext = ".vcf.gz")
    # set genotype
    cmd_setgt <- paste(
        "bcftools +setGT", shQuote(temp_merged), "-- -t . -n 0",
        " | bcftools view - -Oz -o", shQuote(temp_setgt)
    )
    message(glue("Running setGT command:\n  {cmd_setgt}\n"))
    res_setgt <- system(cmd_setgt, ignore.stdout = TRUE, ignore.stderr = FALSE)

    if (res_setgt != 0 || !file.exists(temp_setgt) || file.size(temp_setgt) == 0) {
        stop("Error: Failed to set genotype to reference in final merged VCF")
    }
    message("Indexing temporary setGT VCF...")
    system(paste("tabix -p vcf", shQuote(temp_setgt)), ignore.stdout = TRUE, ignore.stderr = TRUE)
    # Final step: fill-tag
    cmd_filltags <- paste(
        "bcftools +fill-tags", shQuote(temp_setgt),
        "-Oz -o", shQuote(ultimate_final_output),
        "-- -t AC,AN,AF,AC_Hom,AC_Het,AC_Hemi,NS"
    )

    cat("Running fill-tags command:\n  ", cmd_filltags, "\n", sep = "")
    filltags_result <- system(cmd_filltags, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Clean up temporary files
    unlink(temp_merged)
    unlink(paste0(temp_merged, ".tbi"))
    unlink(temp_setgt)
    unlink(paste0(temp_setgt, ".tbi"))
    unlink(temp_annotated)
    unlink(paste0(temp_annotated, ".tbi"))
    unlink(file_list)

    if (filltags_result == 0 && file.exists(ultimate_final_output) && file.size(ultimate_final_output) > 0) {
        # Index final output
        system(paste("tabix -p vcf", shQuote(ultimate_final_output)), ignore.stdout = TRUE, ignore.stderr = TRUE)

        final_size <- file.info(ultimate_final_output)$size
        cat("Ultimate final merged file created successfully; size (bytes):", final_size, "\n")
        return(ultimate_final_output)
    } else {
        cat("Error: Failed to create ultimate final merged file\n")
        return(NULL)
    }
}

# Command line interface
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 3) {
        cat("Usage: Rscript merge_vcfs.R <vcf_input> <output_folder> <region_chunk_size> [join_chunks] [format_fields] [n_cores] [vcf_batch_size]\n")
        cat("  vcf_input: Path to folder / text file containing VCF files\n")
        cat("  output_folder: Path to output folder\n")
        cat("  region_chunk_size: Numeric value for chunk size (bp)\n")
        cat("  join_chunks: TRUE/FALSE (optional, default: FALSE)\n")
        cat("  format_fields: Comma-separated FORMAT fields (optional, default: GT,AD)\n")
        cat("  n_cores: Number of cores for parallel processing (optional, default: 1)\n")
        cat("  vcf_batch_size: Number of VCFs to process per batch (optional, default: process all together)\n")
        quit(status = 1)
    }

    vcf_input <- args[1]
    output_folder <- args[2]
    region_chunk_size <- as.numeric(args[3])
    join_chunks <- if (length(args) >= 4) as.logical(args[4]) else FALSE
    format_fields <- if (length(args) >= 5 && nzchar(args[5])) strsplit(args[5], ",")[[1]] else c("GT", "AD")
    n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1
    vcf_batch_size <- if (length(args) >= 7) as.integer(args[7]) else NULL

    process_vcf_files(vcf_input, output_folder, region_chunk_size, join_chunks, format_fields, n_cores, vcf_batch_size)
}
