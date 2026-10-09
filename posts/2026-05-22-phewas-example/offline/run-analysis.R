#!/usr/bin/env Rscript
# Run from posts/2026-05-22-phewas-example. Optional input/output directories.
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "."
output_dir <- if (length(args) > 1L) args[2] else "chd-followup-results"
source("offline/R/offline-pipeline.R")
disable_remote_access()
source("offline/R/chd-inputs.R")
result <- run_offline_analysis(reconstruct_chd_inputs(input_dir), output_dir)
print(result$summary$families)
