#!/usr/bin/env Rscript

# Required library
if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "http://cran.us.r-project.org")
}
library(optparse)
library(glue)
source("../lib/merge_vcf_helper_functions.R")

# Define command line options
option_list <- list(
    make_option(c("-b", "--batch_size"), type = "integer", help = "Batch size (positive integer)", metavar = "number"),
    make_option(c("-i", "--input_path"), type = "character", help = "Input path: directory or text file with file paths", metavar = "path"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for batch files", metavar = "path"),
    make_option(c("-c", "--chunk_size"), type = "integer", default = 30000000, help = "Chunk size for VCF files (default: 1000000)", metavar = "number"),
    make_option(c("-j", "--join_chunks"), type = "logical", default = TRUE, help = "Join chunks after processing (default: TRUE)", metavar = "TRUE/FALSE"),
    make_option(c("-f", "--fields"), type = "character", default = "GT,AD", help = "Fields to retain in merged VCF (default: GT,AD)", metavar = "fields"),
    make_option(c("-n", "--n_cores"), type = "integer", default = 10, help = "Number of cores to use for merging (default: 1)", metavar = "number"),
    make_option(c("-M", "--only_merge_final", type = "logical", default = FALSE, metavar = "TRUE/FALSE"))
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$batch_size) || opt$batch_size <= 0) {
    stop("Error: You must specify a positive integer for batch size with -b or --batch_size.")
}
if (is.null(opt$input_path)) {
    stop("Error: You must specify an input path with -i or --input_path.")
}
if (is.null(opt$output_dir)) {
    stop("Error: You must specify an output directory with -o or --output_dir.")
}

batch_size <- opt$batch_size
input_path <- opt$input_path
output_dir <- opt$output_dir
chunk_size <- opt$chunk_size
join_chunks <- as.logical(opt$join_chunks)
fields <- opt$fields
n_cores <- opt$n_cores
only_merge_final <- as.logical(opt$only_merge_final)


if (!only_merge_final) {
    # Prepare output directory
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        dir.create(glue("{output_dir}/vcf_batch_files"), recursive = TRUE)
    }


    # Gather file paths
    if (dir.exists(input_path)) {
        files <- list.files(input_path, full.names = TRUE, pattern = "\\.vcf\\.gz$")
    } else if (file.exists(input_path)) {
        files <- readLines(input_path)
    } else {
        stop("Error: input_path must be an existing directory or a file listing paths.")
    }

    n_files <- length(files)
    n_batches <- ceiling(n_files / batch_size)

    if (n_files == 0) {
        message("No files found in input; nothing to do.")
        quit(status = 0)
    }

    # Split into batches and write out
    for (i in seq_len(n_batches)) {
        start_idx <- (i - 1) * batch_size + 1
        end_idx <- min(i * batch_size, n_files)
        batch <- files[start_idx:end_idx]
        message(sprintf("=== Processing batch %d: %d files", i, length(batch)))
        out_file <- file.path(glue("{output_dir}/vcf_batch_files"), sprintf("batch_%03d.txt", i))
        writeLines(batch, out_file)
        message(sprintf("Wrote %d paths to %s", length(batch), out_file))
    }

    # merge each batch
    batch_files <- list.files(glue("{output_dir}/vcf_batch_files"), full.names = TRUE)

    for (i in 1:length(batch_files)) {
        batch_file <- batch_files[i]
        message("++++++++++++++++++++++++++++++++++++++++++++++")
        message(sprintf("Running on batch file: %s", batch_file))
        message("++++++++++++++++++++++++++++++++++++++++++++++")
        merge_vcfs_command <- glue("./merge_vcfs.R {batch_file} {output_dir} {chunk_size} {join_chunks} {fields} {n_cores} batch_{i}")
        system(merge_vcfs_command)
        write("===== Full batch analysis command (run this to retry running the batch) =====", file = glue("{output_dir}/batch_{i}_commands.log"), append = TRUE)
        write(merge_vcfs_command, file = glue("{output_dir}/batch_{i}_commands.log"), append = TRUE)
    }
}

# merge final batch files
final_batch_files <- list.files(glue("{output_dir}"), full.names = TRUE, pattern = "batch_\\d+_final_merged.+\\.vcf\\.gz$")
final_batch_files_list <- glue("{output_dir}/final_batch_files.txt")
writeLines(final_batch_files, file.path(final_batch_files_list))
merge_vcf_files(final_batch_files_list, output_dir)
