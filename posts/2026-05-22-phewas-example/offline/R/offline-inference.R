# Each BH family is separate; unavailable tests do not enter its denominator.
bh_available <- function(p) {
  result <- rep(NA_real_, length(p))
  available <- is.finite(p) & p >= 0 & p <= 1
  result[available] <- p.adjust(p[available], method = "BH")
  result
}

assert_unique_pairs <- function(x) {
  if (anyDuplicated(x[c("mediator_id", "target_snp")]))
    stop("Analytical summaries require one row per unique mediator/target pair")
  invisible(x)
}

annotate_comparisons <- function(x) {
  assert_unique_pairs(x)
  x$q_family <- "unique_pair_heterogeneity"
  x$q_padj <- bh_available(x$q_pval)
  x$q_family_inventory <- nrow(x)
  x$q_family_tested <- sum(is.finite(x$q_padj))
  x$q_family_unavailable <- sum(!is.finite(x$q_padj))
  x$interpretation <- ifelse(!is.finite(x$q_padj), "Unavailable comparison",
    ifelse(x$q_padj < 0.05, "Incompatible with common-effect model (BH < 0.05)",
      "No detected heterogeneity; mediation unresolved"))
  x$precision <- ifelse(!is.finite(x$comparator_lo) | !is.finite(x$comparator_hi),
    "Comparator unavailable", ifelse(x$comparator_lo <= 0 & x$comparator_hi >= 0,
      "Comparator CI includes zero", "Comparator CI excludes zero"))
  x
}

annotate_decomposition <- function(x) {
  assert_unique_pairs(x)
  x$i_family <- "unique_pair_indirect_component"
  x$d_family <- "unique_pair_residual_component"
  x$i_padj <- bh_available(x$i_pval)
  x$d_padj <- bh_available(x$d_pval)
  for (component in c("i", "d")) {
    available <- is.finite(x[[paste0(component, "_padj")]])
    x[[paste0(component, "_family_inventory")]] <- nrow(x)
    x[[paste0(component, "_family_tested")]] <- sum(available)
    x[[paste0(component, "_family_unavailable")]] <- sum(!available)
  }
  x
}

summarise_inference <- function(comparisons, pooled, decomposition) {
  assert_unique_pairs(comparisons)
  assert_unique_pairs(decomposition)
  if (anyDuplicated(pooled$mediator_id)) stop("Pooled inference requires unique mediator IDs")
  family <- function(name, p) data.frame(family = name, inventory = length(p),
    tested = sum(is.finite(bh_available(p))), unavailable = sum(!is.finite(bh_available(p))),
    nominal_below_005 = sum(is.finite(bh_available(p)) & p < 0.05, na.rm = TRUE),
    bh_below_005 = sum(bh_available(p) < 0.05, na.rm = TRUE))
  families <- rbind(family("unique_pair_heterogeneity", comparisons$q_pval),
    family("pooled_mediator_mr", pooled$pval),
    family("unique_pair_indirect_component", decomposition$i_pval),
    family("unique_pair_residual_component", decomposition$d_pval))
  per_snp <- do.call(rbind, lapply(unique(comparisons$target_snp), function(snp) {
    x <- comparisons[comparisons$target_snp == snp, ]
    data.frame(target_snp = snp, pairs = nrow(x), tested = sum(is.finite(x$q_pval)),
      unavailable = sum(!is.finite(x$q_pval)),
      incompatible_bh = sum(x$q_padj < 0.05, na.rm = TRUE),
      no_detected_heterogeneity = sum(x$q_padj >= 0.05, na.rm = TRUE))
  }))
  list(families = families, per_snp = per_snp)
}
