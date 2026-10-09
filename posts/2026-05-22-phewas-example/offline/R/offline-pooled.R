# Adapted from MRCIEU/gbc-mr 3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
# Pooled mediator MR uses the exact same shared fits as the comparison panels.
# Pass comparisons from compute_comparisons to reuse fits without recomputing.
compute_pooled <- function(inputs, window_bp = 1e6, mhc = FALSE,
                           comparisons = NULL, comparators = NULL) {
  pairs <- canonical_pairs(inputs)
  if (!is.null(comparisons)) {
    if (is.null(attr(comparisons, "comparators")))
      stop("comparisons must retain its shared comparator attribute")
    if (any(comparisons$window_bp != window_bp) || any(comparisons$mhc_exclusion != mhc))
      stop("Shared comparator settings do not match pooled settings")
    comparators <- attr(comparisons, "comparators")
  }
  if (is.null(comparators)) comparators <- compute_shared_comparators(inputs, window_bp, mhc)
  mediators <- unique(pairs[, c("mediator_id", "mediator", "mediator_scale",
                               "outcome_id", "outcome_name", "outcome_units", "mr_effect_scale")])
  if (anyDuplicated(mediators$mediator_id)) stop("Conflicting mediator labels or scales")
  rows <- lapply(seq_len(nrow(mediators)), function(i) {
    med <- mediators[i, , drop = FALSE]
    shared <- comparators[[med$mediator_id]]
    if (is.null(shared) || shared$window_bp != window_bp || shared$mhc != mhc)
      stop("Missing or incompatible shared comparator for ", med$mediator_id)
    count <- shared$counts
    hc <- shared$harmonisation_counts
    selected <- shared$selection
    targets <- unique(pairs$target_snp[pairs$mediator_id == med$mediator_id])
    med$n_targets <- length(targets)
    med$target_snps <- paste(targets, collapse = ";")
    med$window_bp <- window_bp
    med$mhc_exclusion <- mhc
    med$exclusion_family <- "pooled_mediator_mr"
    med$n_saved_before_exclusions <- count$n_saved
    med$n_retained_before_outcome_matching <- count$n_retained
    matched <- if (!is.null(selected)) selected$instruments$SNP %in% inputs$outgwasf$SNP else logical()
    med$n_outcome_matched_before_harmonisation <- if (!is.null(selected)) sum(matched) else NA_integer_
    med$n_harmonised <- hc$n_usable
    med$n_removed_during_harmonisation <- if (!is.null(selected)) sum(matched) - hc$n_usable else NA_integer_
    med$n_harmonisation_missing_outcome <- hc$n_missing_outcome
    med$n_harmonisation_palindromic <- hc$n_palindromic
    med$n_harmonisation_incompatible <- hc$n_incompatible
    for (reason in c("region", "target", "mhc", "coordinate_missing", "coordinate_mismatch")) {
      med[[paste0("n_excluded_", reason, "_after_outcome_matching")]] <- if (!is.null(selected))
        sum(selected$excluded_instruments$SNP %in% inputs$outgwasf$SNP &
          selected$excluded_instruments$exclusion_reason == reason) else NA_integer_
    }
    med$n_excluded_after_outcome_matching <- if (!is.null(selected))
      sum(selected$excluded_instruments$SNP %in% inputs$outgwasf$SNP) else NA_integer_
    # The existing selector count is the single authoritative count here.
    count <- count[, setdiff(names(count), names(med)), drop = FALSE]
    cbind(med, count, hc, shared$estimate)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$family <- "pooled_mediator_mr"
  result$pval_bh <- NA_real_
  available <- result$status == "available" & is.finite(result$pval) & result$pval >= 0 & result$pval <= 1
  result$pval_bh[available] <- stats::p.adjust(result$pval[available], method = "BH")
  result$family_inventory <- nrow(result)
  result$family_tested <- sum(available)
  result$family_unavailable <- sum(!available)
  result
}

# Diagnostic-only comparison with the pooled MR values printed in the original
# stored HTML report. The historical exposure label sometimes contains its
# OpenGWAS ID; a plain display-label fallback is explicitly marked when used.
compare_pooled_legacy_diagnostic <- function(pooled, legacy) {
  stopifnot(all(c("mediator_id", "b", "se", "pval") %in% names(pooled)),
            all(c("exposure", "b", "se", "pval") %in% names(legacy)))
  legacy_id <- if ("mediator_id" %in% names(legacy)) as.character(legacy$mediator_id) else rep(NA_character_, nrow(legacy))
  legacy_id[!nzchar(legacy_id)] <- NA_character_
  mediator_name <- if ("mediator" %in% names(pooled)) as.character(pooled$mediator) else rep(NA_character_, nrow(pooled))
  id_match <- match(as.character(pooled$mediator_id), legacy_id)
  label_match <- match(mediator_name, as.character(legacy$exposure))
  index <- id_match
  method <- rep("unmatched", nrow(pooled))
  method[!is.na(id_match)] <- "mediator_id"
  fallback <- is.na(index) & !is.na(label_match)
  index[fallback] <- label_match[fallback]
  method[fallback] <- "display_label_only_no_legacy_id"
  data.frame(
    mediator_id = pooled$mediator_id,
    mediator = mediator_name,
    legacy_exposure = ifelse(is.na(index), NA_character_, as.character(legacy$exposure[index])),
    matching_method = method,
    corrected_b = pooled$b,
    corrected_se = pooled$se,
    corrected_pval = pooled$pval,
    corrected_nsnp = pooled$nsnp,
    legacy_b = ifelse(is.na(index), NA_real_, as.numeric(legacy$b[index])),
    legacy_se = ifelse(is.na(index), NA_real_, as.numeric(legacy$se[index])),
    legacy_pval = ifelse(is.na(index), NA_real_, as.numeric(legacy$pval[index])),
    legacy_nsnp = ifelse(is.na(index), NA_integer_, as.integer(legacy$nsnp[index])),
    delta_b = ifelse(is.na(index), NA_real_, pooled$b - as.numeric(legacy$b[index])),
    stringsAsFactors = FALSE
  )
}
