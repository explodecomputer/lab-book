.plot_fixture <- function() {
  sizes <- c(7, 10, 1, 2, 3)
  x <- data.frame(mediator_id = rep(c("ieu-b-5138", "ukb-d-30780_raw", "ukb-b-10831", "ukb-b-5779", "ukb-b-2303"), sizes),
    mediator = rep(c("SBP", "LDL", "Smoking", "Alcohol", "BMI"), sizes),
    target_snp = paste0("rs", seq_len(sum(sizes))), mediator_scale = "Saved GWAS units; scale unavailable",
    outcome_name = "CHD", outcome_units = "Saved GWAS units; scale unavailable",
    effect_allele = "A", a = 0.2, a_se = 0.02, t = 0.1, t_se = 0.01,
    comparator_b = 0.5, comparator_se = 0.05, comparator_lo = 0.402,
    comparator_hi = 0.598, comparator_nsnp = 100, comparator_status = "available",
    target_b = 0.6, target_lo = 0.45, target_hi = 0.75,
    q = 0.1, q_pval = 0.75, window_bp = 1e6, n_exclusion_loci = 41)
  x$a[1] <- x$target_b[1] <- x$target_lo[1] <- x$target_hi[1] <- x$q[1] <- x$q_pval[1] <- NA_real_
  annotate_comparisons(x)
}

testthat::test_that("plot exports preserve complete page inventories and unknown units", {
  x <- .plot_fixture()
  directory <- tempfile("offline-plots-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  comparisons <- export_comparison_plot(x, file.path(directory, "comparisons.pdf"))
  decomposition <- export_decomposition_plot(annotate_decomposition(compute_decomposition(x)), file.path(directory, "decomposition.pdf"))
  pooled <- x[!duplicated(x$mediator_id), ]
  pooled$b <- pooled$comparator_b; pooled$lo <- pooled$comparator_lo; pooled$hi <- pooled$comparator_hi
  pooled$nsnp <- pooled$comparator_nsnp; pooled$method <- "IVW (multiplicative random effects)"
  pooled$pval <- 0.01; pooled$pval_bh <- 0.02
  pooled$n_excluded_target <- 3; pooled$n_excluded_region <- 4
  pooled_pngs <- export_pooled_plot(pooled, file.path(directory, "pooled.pdf"))
  testthat::expect_length(comparisons, 5)
  testthat::expect_length(decomposition, 4)
  testthat::expect_length(pooled_pngs, 1)
  files <- c(comparisons, decomposition, pooled_pngs,
    file.path(directory, c("comparisons.pdf", "decomposition.pdf", "pooled.pdf")))
  testthat::expect_true(all(file.exists(files)))
  testthat::expect_true(all(file.info(files)$size > 1000))
  testthat::expect_equal(basename(comparisons), c("comparisons.png", sprintf("comparisons-page-%02d.png", 2:5)))
  testthat::expect_match(.plot_units(x, "outcome_units", "outcome"), "scale unknown")
  testthat::expect_false(grepl("GBC|log odds", .plot_units(x, "outcome_units", "outcome")))
})

testthat::test_that("comparison plots reject multiple comparators within a mediator", {
  x <- .plot_fixture()
  x$comparator_b[2] <- 0.8
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  testthat::expect_error(export_comparison_plot(x, path), "one shared comparator")
})
