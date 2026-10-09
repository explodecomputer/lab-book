#!/usr/bin/env Rscript
source("offline/R/offline-pipeline.R")
source_offline_modules()
source("offline/R/chd-inputs.R")
testthat::test_dir("offline/tests", reporter = "summary", stop_on_failure = TRUE)
