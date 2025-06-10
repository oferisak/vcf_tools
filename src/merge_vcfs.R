#!/usr/bin/env Rscript

# !/usr/bin/env Rscript

################################################################################
#                         VCF Processing Script                               #
################################################################################
#
# DESCRIPTION:
#   This script processes multiple VCF (Variant Call Format) files by:
#   1. Chunking each VCF file into genomic regions
#   2. Merging chunks across samples for each region
#   3. Optionally concatenating all regions into a single final file
#
#   The script is designed for large-scale genomic data processing with
#   parallelization support and handles both compressed (.vcf.gz) and
#   uncompressed (.vcf) files.
#
# FEATURES:
#   - Parallel processing support (configurable number of cores)
#   - Automatic VCF compression and indexing
#   - Filters for PASS variants only
#   - Customizable FORMAT field retention
#   - Proper genomic coordinate sorting
#   - Automatic cleanup of intermediate files
#
# REQUIREMENTS:
#   - R packages: parallel
#   - External tools: bcftools, tabix, bgzip (from htslib/samtools suite)
#   - Operating System: Linux/Unix/macOS with bash shell
#
# INSTALLATION OF DEPENDENCIES:
#   # Install R (if not already installed)
#   # On Ubuntu/Debian: sudo apt-get install r-base
#   # On CentOS/RHEL: sudo yum install R
#   # On macOS: brew install r
#
#   # Install bcftools and related tools
#   # On Ubuntu/Debian: sudo apt-get install bcftools tabix
#   # On CentOS/RHEL: sudo yum install bcftools htslib
#   # On macOS: brew install bcftools htslib
#
#   # Or compile from source: http://www.htslib.org/download/
#
# USAGE:
#   Rscript vcf_processor.R <vcf_input> <chunk_size> [join_chunks] [format_fields] [n_cores]
#
# PARAMETERS:
#   vcf_input      : Path to directory containing VCF files (.vcf or .vcf.gz)
#   output_folder   : Path to directory where output files will be saved
#   chunk_size      : Size of genomic regions in base pairs (e.g., 1000000 for 1Mb)
#   join_chunks     : TRUE/FALSE - whether to concatenate all regions (default: FALSE)
#   format_fields   : Comma-separated FORMAT fields to keep (default: "GT,AD")
#   n_cores         : Number of CPU cores for parallel processing (default: 1)
#
# EXAMPLES:
#   # Basic usage -  Custom FORMAT fields and parallel processing
#   Rscript vcf_processor.R /path/to/vcf/files /path/to/output/folder 500000 TRUE "GT,AD,DP,GQ" 4

# OUTPUT:
#   - chunks/           : Temporary directory with individual region chunks (deleted if empty)
#   - merged_chunks/    : Directory with merged VCF files per region
#   - final_merged_all_samples_all_regions.vcf.gz : Final concatenated file (if join_chunks=TRUE)
#
# WORKFLOW:
#   1. Input validation and tool checking
#   2. VCF file discovery, compression, and indexing
#   3. Chromosome/contig information extraction
#   4. Region generation based on chunk size
#   5. Parallel chunking of each VCF by region (filters PASS variants only)
#   6. Parallel merging of chunks across samples for each region
#   7. FORMAT field filtering and INFO field recalculation
#   8. Optional concatenation of all regions with proper genomic sorting
#
# PERFORMANCE NOTES:
#   - Processing time scales with: number of samples, genome size, variant density
#   - Memory usage is generally low due to streaming processing
#   - Disk space needed: ~2-3x input file size for intermediate files
#   - Recommended chunk sizes: 100kb-10Mb depending on variant density
#
# TROUBLESHOOTING:
#   - "Required tool not found": Install bcftools/tabix/bgzip
#   - "No VCF files found": Check file extensions (.vcf or .vcf.gz)
#   - "Permission denied": Ensure write access to output directory
#   - Out of disk space: Reduce chunk_size or free up space
#   - Memory issues: Reduce n_cores parameter

library(parallel)
library(glue)

process_vcf_files <- function(vcf_input, output_folder, region_chunk_size, join_chunks = FALSE,
                              format_fields = c("GT", "AD"), n_cores = 1) {
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
            if (grepl("_|\-", chrom)) {
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

    cat("Getting chromosome information...\n")
    chrom_info <- get_chromosomes(vcf_files[1])
    regions <- generate_regions(chrom_info, region_chunk_size)
    cat("Generated", length(regions), "region(s) for processing\n\n")

    # Function to chunk a single VCF file for a specific region
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

    # Function to merge chunks for a region
    merge_region_chunks <- function(chunk_files, region) {
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

        # Post-processing with fill-tags
        temp_fill <- glue("{merged_file}.tmp.vcf.gz")
        # temp_fill <- tempfile(fileext = ".vcf.gz")
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
        batch_chunk_results <- list()

        for (region in current_regions) {
            cat("    Chunking region:", region, "...\n")
            region_chunks_raw <- mclapply(
                vcf_files,
                function(v) chunk_vcf(v, region),
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

            merged_file <- merge_region_chunks(chunk_files, region)
            if (!is.null(merged_file)) {
                cat("    Successfully merged region:", region, "\n")
            } else {
                cat("    Warning: Failed to merge region:", region, "\n")
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

        final_output <- file.path(output_folder, "final_merged_all_samples_all_regions.vcf.gz")
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
        final_file <- file.path(output_folder, "final_merged_all_samples_all_regions.vcf.gz")
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

    process_vcf_files(vcf_input, output_folder, region_chunk_size, join_chunks, format_fields, n_cores)
}
