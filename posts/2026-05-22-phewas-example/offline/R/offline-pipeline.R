# Repo-local, outcome-neutral runner. Statistical modules are derived from
# MRCIEU/gbc-mr@3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
source_offline_modules <- function(module_dir = "offline/R") {
  for (module in c("core", "pooled", "decomposition", "inference", "plots"))
    source(file.path(module_dir, paste0("offline-", module, ".R")), local = globalenv())
}

disable_remote_access <- function() {
  Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
             ALL_PROXY = "http://127.0.0.1:9", NO_PROXY = "")
  block <- quote(stop("Remote access disabled: this report uses local inputs only"))
  for (package in c("TwoSampleMR", "ieugwasr")) {
    if (requireNamespace(package, quietly = TRUE)) {
      functions <- if (package == "TwoSampleMR") c("extract_instruments", "extract_outcome_data") else
        c("api_query", "tophits", "phewas", "gwasinfo", "gwasinfo_files", "variants_rsid", "ld_matrix", "ld_clump")
      for (fun in functions) if (exists(fun, envir = asNamespace(package), inherits = FALSE))
        invisible(capture.output(trace(fun, tracer = block, where = asNamespace(package), print = FALSE)))
    }
  }
  for (fun in c("download.file", "download.packages", "install.packages"))
    invisible(capture.output(trace(fun, tracer = block, where = asNamespace("utils"), print = FALSE)))
  invisible(TRUE)
}

run_offline_analysis <- function(inputs, output_dir = "chd-followup-results",
                                 figures_dir = file.path(output_dir, "figures"),
                                 module_dir = "offline/R") {
  source_offline_modules(module_dir)
  disable_remote_access()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  comparisons <- annotate_comparisons(compute_comparisons(inputs))
  pooled <- compute_pooled(inputs, comparisons = comparisons)
  decomposition <- annotate_decomposition(compute_decomposition(comparisons))
  summary <- summarise_inference(comparisons, pooled, decomposition)
  sensitivity <- do.call(rbind, lapply(c(FALSE, TRUE), function(mhc) {
    x <- annotate_comparisons(compute_comparisons(inputs, window_bp = 2e6, mhc = mhc))
    x$analysis <- if (mhc) "+/-2 Mb plus extended MHC chr6 25-34 Mb" else "+/-2 Mb"
    x$q_family <- if (mhc) "sensitivity_pair_heterogeneity_2mb_mhc" else "sensitivity_pair_heterogeneity_2mb"
    x
  }))
  covariance <- decomposition_covariance_sensitivity(comparisons)
  covariance <- do.call(rbind, lapply(sort(unique(covariance$rho)), function(rho) {
    x <- annotate_decomposition(covariance[covariance$rho == rho, ])
    suffix <- paste0("_hypothetical_rho_", rho)
    x$i_family <- paste0(x$i_family, suffix)
    x$d_family <- paste0(x$d_family, suffix)
    x
  }))
  add_family <- function(name, p) {
    data.frame(family = name, inventory = length(p), tested = sum(is.finite(p)),
      unavailable = sum(!is.finite(p)), nominal_below_005 = sum(p < 0.05, na.rm = TRUE),
      bh_below_005 = sum(bh_available(p) < 0.05, na.rm = TRUE))
  }
  for (name in unique(sensitivity$q_family))
    summary$families <- rbind(summary$families, add_family(name, sensitivity$q_pval[sensitivity$q_family == name]))
  for (component in c("i", "d")) {
    family_column <- paste0(component, "_family")
    for (name in unique(covariance[[family_column]]))
      summary$families <- rbind(summary$families,
        add_family(name, covariance[[paste0(component, "_pval")]][covariance[[family_column]] == name]))
  }
  leads <- inputs$lead_inventory
  leads$retained_pairs <- vapply(leads$target_snp, function(snp) sum(comparisons$target_snp == snp), integer(1))
  leads$outcome_available <- leads$target_snp %in% inputs$outgwasf$SNP
  leads$followup_status <- ifelse(leads$retained_pairs > 0, "Retained exact follow-up pair(s)",
    "No retained saved mediator association; still excluded from every comparator")
  tables <- list(comparisons = comparisons, pooled_mr = pooled, decomposition = decomposition,
    window_sensitivity = sensitivity, covariance_sensitivity = covariance,
    inference_families = summary$families, per_snp_summary = summary$per_snp,
    lead_inventory = leads, exact_associations = inputs$int_chd,
    instrument_associations = inputs$int_inst, outcome_associations = inputs$outgwasf,
    provenance = inputs$provenance, extraction_counts = inputs$extraction_counts)
  source_version <- "3bbac89e254acfbfa3385aae9b4520e4e8b5932c"
  tables$implementation_provenance <- data.frame(source = "MRCIEU/gbc-mr", commit = source_version,
    adaptation = "Repo-local outcome metadata, complete lead inventory, local VCF reconstruction and readable plots")
  packages <- c("TwoSampleMR", "ieugwasr", "gwasvcf", "GenomicRanges", "Rsamtools", "VariantAnnotation", "knitr", "rmarkdown", "testthat")
  tables$package_versions <- data.frame(package = c("R", packages),
    version = c(as.character(getRversion()), vapply(packages, function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "not installed", character(1))))
  stopifnot(nrow(comparisons) == nrow(canonical_pairs(inputs)),
    nrow(decomposition) == nrow(comparisons),
    all(abs(decomposition$t - decomposition$i - decomposition$d) < 1e-10, na.rm = TRUE))
  pi <- match(comparisons$mediator_id, pooled$mediator_id)
  stopifnot(identical(comparisons$comparator_b, pooled$b[pi]),
    identical(comparisons$comparator_se, pooled$se[pi]),
    identical(comparisons$comparator_nsnp, pooled$nsnp[pi]))
  for (name in names(tables)) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    write.csv(tables[[name]], path, row.names = FALSE, na = "")
    stopifnot(nrow(read.csv(path)) == nrow(tables[[name]]))
  }
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
    file.path(output_dir, "session-info.txt"))
  plots <- list(comparisons = export_comparison_plot(comparisons, file.path(figures_dir, "comparisons.pdf")),
    decomposition = export_decomposition_plot(decomposition, file.path(figures_dir, "decomposition.pdf")),
    pooled = export_pooled_plot(pooled, file.path(figures_dir, "pooled.pdf")))
  list(inputs = inputs, comparisons = comparisons, pooled = pooled, decomposition = decomposition,
    summary = summary, tables = tables, plots = plots)
}
