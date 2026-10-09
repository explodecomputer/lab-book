# Adapted from MRCIEU/gbc-mr 3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
# Outcome-neutral offline analysis. SNP and GWAS IDs are keys; labels never identify groups.
# Positions are saved coordinates, not a verified LD reference or build conversion.
load_saved_inputs <- function(path = "data/gbc_phewas_results-g2.Rdata") {
  saved <- new.env(parent = emptyenv())
  load(path, envir = saved)
  required <- c("int_chd", "int_inst", "outgwasf", "phewas_result")
  stopifnot(all(required %in% ls(saved)))
  wanted <- unique(c(saved$int_inst$SNP, saved$int_chd$SNP,
                     saved$phewas_result$rsid))
  # Retain exact targets from the full outcome, never from reduced outdat.
  outcome <- as.data.frame(saved$outgwasf[saved$outgwasf$SNP %in% wanted, ])
  list(int_chd = as.data.frame(saved$int_chd),
       int_inst = as.data.frame(saved$int_inst), outgwasf = outcome,
       phewas_result = as.data.frame(saved$phewas_result),
       inventory_counts = data.frame(n_saved_pair_rows = nrow(saved$int_chd),
                                     n_saved_instrument_rows = nrow(saved$int_inst),
                                     n_full_outcome_rows = nrow(saved$outgwasf),
                                     n_retained_outcome_rows = nrow(outcome),
                                     n_original_phewas_rows = nrow(saved$phewas_result)))
}

offline_deduplicate <- function(dat, role) {
  dat <- as.data.frame(dat)
  if (!nrow(dat)) return(list(data = dat, removed = 0L, conflicts = character()))
  key <- paste(dat[[paste0("id.", role)]], dat$SNP, sep = "\r")
  important <- intersect(c(paste0(c("beta", "se", "effect_allele", "other_allele",
                                    "chr", "pos", "mr_keep"), ".", role),
                            "chr", "pos", paste0(c("proxy", "target_snp", "build"), ".", role), "build"),
                         names(dat))
  conflict_keys <- unique(key[duplicated(key)])
  conflict_keys <- conflict_keys[vapply(conflict_keys, function(k) {
    nrow(unique(dat[key == k, important, drop = FALSE])) > 1L
  }, logical(1))]
  keep <- !duplicated(key) & !key %in% conflict_keys
  list(data = dat[keep, , drop = FALSE],
       removed = sum(duplicated(key)),
       conflicts = unique(dat$SNP[key %in% conflict_keys]))
}

offline_chr <- function(x) {
  x <- toupper(sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE))
  x[x %in% c("", "NA", "UNKNOWN")] <- NA_character_
  x[x == "23" & !is.na(x)] <- "X"
  x[x == "24" & !is.na(x)] <- "Y"
  x[x %in% c("M", "26")] <- "MT"
  x
}

offline_valid_coords <- function(chr, pos) {
  !is.na(chr) & chr %in% c(as.character(1:22), "X", "Y", "MT") &
    is.finite(pos) & pos > 0 & pos == floor(pos)
}

offline_coord <- function(dat, role) {
  chr_col <- if (paste0("chr.", role) %in% names(dat)) paste0("chr.", role) else "chr"
  pos_col <- if (paste0("pos.", role) %in% names(dat)) paste0("pos.", role) else "pos"
  chr <- if (chr_col %in% names(dat)) as.character(dat[[chr_col]]) else rep(NA_character_, nrow(dat))
  chr <- offline_chr(chr)
  pos <- if (pos_col %in% names(dat)) suppressWarnings(as.numeric(dat[[pos_col]])) else rep(NA_real_, nrow(dat))
  data.frame(chr = chr, pos = pos, stringsAsFactors = FALSE)
}

offline_outcome_metadata <- function(inputs) {
  metadata <- inputs$outcome_metadata
  if (is.null(metadata)) {
    # Compatibility for the pinned reference's saved GBC bundle.
    out <- inputs$outgwasf
    metadata <- list(name = if ("outcome" %in% names(out)) unique(out$outcome)[1] else "Outcome",
      id = unique(out$id.outcome)[1], units = "Saved GWAS units; scale unavailable", build = NA_character_)
  }
  for (field in c("name", "id", "units")) {
    if (length(metadata[[field]]) != 1L || is.na(metadata[[field]]) || !nzchar(metadata[[field]]))
      stop("outcome_metadata requires a nonempty ", field)
  }
  ids <- unique(as.character(inputs$outgwasf$id.outcome))
  ids <- ids[!is.na(ids)]
  if (length(ids) && (length(ids) != 1L || ids != metadata$id))
    stop("Outcome associations do not match outcome_metadata$id")
  metadata
}

offline_build <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("hg19", "b37", "grch37", "37")] <- "GRCh37"
  x[x %in% c("hg38", "b38", "grch38", "38")] <- "GRCh38"
  x[x %in% c("", "unknown", "na")] <- NA_character_
  x
}

offline_data_build <- function(dat, role) {
  field <- intersect(c(paste0("build.", role), "build"), names(dat))
  if (!length(field)) return(rep(NA_character_, nrow(dat)))
  offline_build(dat[[field[1]]])
}

canonical_pairs <- function(inputs) {
  target <- inputs$int_chd
  pairs <- unique(data.frame(mediator_id = as.character(target$id.outcome),
                             target_snp = as.character(target$SNP)))
  if (anyNA(pairs) || any(!nzchar(pairs$mediator_id) | !nzchar(pairs$target_snp)))
    stop("Pair inventory requires exact mediator IDs and SNPs")
  ix <- match(paste(pairs$mediator_id, pairs$target_snp, sep = "\r"),
              paste(target$id.outcome, target$SNP, sep = "\r"))
  pairs$pair_id <- paste(pairs$mediator_id, pairs$target_snp, sep = "::")
  pairs$mediator <- if ("originalname.outcome" %in% names(target)) {
    as.character(target$originalname.outcome[ix])
  } else if ("outcome" %in% names(target)) as.character(target$outcome[ix]) else pairs$mediator_id
  pairs$mediator_scale <- if ("units.outcome" %in% names(target)) as.character(target$units.outcome[ix]) else NA_character_
  missing_scale <- is.na(pairs$mediator_scale) | !nzchar(pairs$mediator_scale)
  pairs$mediator_scale[missing_scale] <- "Saved GWAS units; specific transformation/unit unavailable"
  loci <- known_lead_loci(inputs)
  li <- match(pairs$target_snp, loci$target_snp)
  pairs$target_chr <- loci$target_chr[li]
  pairs$target_pos <- loci$target_pos[li]
  if (!is.null(inputs$lead_inventory)) {
    # Pair labels/independent total effects use the declared lead coordinates.
    # A conflicting mediator record may invalidate harmonisation/comparators
    # without erasing an independently valid outcome association.
    inventory <- inputs$lead_inventory
    for (i in seq_len(nrow(pairs))) {
      rows <- inventory[inventory$target_snp == pairs$target_snp[i],
                        c("target_chr", "target_pos"), drop = FALSE]
      rows$target_chr <- offline_chr(rows$target_chr)
      rows$target_pos <- suppressWarnings(as.numeric(rows$target_pos))
      rows <- unique(rows)
      if (nrow(rows) == 1L && offline_valid_coords(rows$target_chr, rows$target_pos)) {
        pairs$target_chr[i] <- rows$target_chr
        pairs$target_pos[i] <- rows$target_pos
      }
    }
  }
  meta <- offline_outcome_metadata(inputs)
  pairs$outcome_id <- meta$id
  pairs$outcome_name <- meta$name
  pairs$outcome_units <- meta$units
  pairs$effect_scale <- paste(meta$units, "per target effect allele")
  pairs$mr_effect_scale <- paste(meta$units, "per", pairs$mediator_scale)
  pairs[order(pairs$mediator_id, pairs$target_snp), , drop = FALSE]
}

# A complete inventory is independent of outcome and retained pair coverage.
# Only the legacy pinned GBC bundle falls back to its original PheWAS inventory.
known_lead_loci <- function(inputs) {
  explicit <- !is.null(inputs$lead_inventory)
  if (explicit) {
    inventory <- as.data.frame(inputs$lead_inventory)
    if (!all(c("target_snp", "target_chr", "target_pos") %in% names(inventory)))
      stop("lead_inventory requires target_snp, target_chr and target_pos")
    snps <- sort(unique(as.character(inventory$target_snp)))
    if (anyNA(snps) || any(!nzchar(snps))) stop("Lead SNP identities must be nonempty")
  } else {
    snps <- sort(unique(c(as.character(inputs$phewas_result$rsid), as.character(inputs$int_chd$SNP))))
    snps <- snps[!is.na(snps) & nzchar(snps)]
  }
  loci <- data.frame(target_snp = snps, target_chr = NA_character_, target_pos = NA_real_)
  issues <- character()
  for (i in seq_along(snps)) {
    saved <- inputs$outgwasf[inputs$outgwasf$SNP == snps[i], , drop = FALSE]
    exact <- inputs$int_chd[inputs$int_chd$SNP == snps[i], , drop = FALSE]
    coords <- rbind(offline_coord(saved, "outcome"), offline_coord(exact, "outcome"))
    if (explicit) {
      lead <- inventory[inventory$target_snp == snps[i], , drop = FALSE]
      lc <- data.frame(chr = offline_chr(lead$target_chr),
                       pos = suppressWarnings(as.numeric(lead$target_pos)))
      if (any(!offline_valid_coords(lc$chr, lc$pos))) {
        issues <- c(issues, paste(snps[i], "missing lead coordinates"))
        next
      }
      coords <- rbind(lc, coords)
    }
    coords <- unique(coords[offline_valid_coords(coords$chr, coords$pos), ])
    if (nrow(coords) == 1L) {
      loci$target_chr[i] <- coords$chr
      loci$target_pos[i] <- coords$pos
    } else issues <- c(issues, paste(snps[i], "missing or conflicting coordinates"))
  }
  if (explicit) {
    metadata <- offline_outcome_metadata(inputs)
    declared <- offline_build(inputs$coordinate_build)
    expected <- offline_build(metadata$build)
    lb <- if ("build" %in% names(inventory)) offline_build(inventory$build) else declared
    if (length(expected) != 1L || is.na(expected) || !expected %in% c("GRCh37", "GRCh38") ||
        !length(lb) || anyNA(lb)) issues <- c(issues, "Missing or unsupported coordinate/build declaration")
    if (!length(declared) || anyNA(declared)) {
      for (entry in c("int_inst", "int_chd", "outgwasf")) {
        role <- if (entry == "int_inst") "exposure" else "outcome"
        if (anyNA(offline_data_build(inputs[[entry]], role)))
          issues <- c(issues, paste("Missing association coordinate build:", entry))
      }
    }
    builds <- c(lb, declared, offline_data_build(inputs$int_inst, "exposure"),
                offline_data_build(inputs$int_chd, "outcome"), offline_data_build(inputs$outgwasf, "outcome"))
    builds <- unique(builds[!is.na(builds)])
    if (length(expected) == 1L && (any(!builds %in% expected)))
      issues <- c(issues, "Conflicting coordinate builds; no build conversion is supported")
  }
  attr(loci, "issues") <- unique(issues)
  loci
}

# Compatibility alias; all new callers use the outcome-neutral name.
known_gbc_loci <- known_lead_loci

select_instruments <- function(inputs, mediator_id, target_snps,
                               window_bp = 1e6, mhc = FALSE) {
  stopifnot(length(mediator_id) == 1L, is.finite(window_bp), window_bp >= 0)
  raw <- inputs$int_inst[inputs$int_inst$id.exposure == mediator_id, , drop = FALSE]
  dd <- offline_deduplicate(raw, "exposure")
  inst <- dd$data
  coords <- offline_coord(inst, "exposure")
  loci <- known_lead_loci(inputs)
  if (!all(target_snps %in% loci$target_snp) ||
      !nrow(loci) || length(attr(loci, "issues")) || any(is.na(loci$target_chr) | loci$target_chr == "" | !is.finite(loci$target_pos))) {
    stop("Known lead coordinates unavailable or conflicting: cannot establish geographic exclusions; ",
         paste(attr(loci, "issues"), collapse = "; "))
  }
  reason <- rep(NA_character_, nrow(inst))
  mark <- function(which, value) {
    reason[which & is.na(reason)] <<- value
  }
  mark(inst$SNP %in% loci$target_snp, "target")
  mark(!offline_valid_coords(coords$chr, coords$pos), "coordinate_missing")
  for (i in seq_len(nrow(loci))) {
    regional <- coords$chr == loci$target_chr[i] &
      abs(coords$pos - loci$target_pos[i]) <= window_bp
    regional[is.na(regional)] <- FALSE
    mark(regional, "region")
  }
  if (mhc && any(loci$target_chr == "6")) {
    in_mhc <- coords$chr == "6" & coords$pos >= 25e6 & coords$pos <= 34e6
    in_mhc[is.na(in_mhc)] <- FALSE
    mark(in_mhc, "mhc")
  }
  out <- offline_deduplicate(inputs$outgwasf, "outcome")$data
  j <- match(inst$SNP, out$SNP)
  oc <- offline_coord(out, "outcome")
  mismatch <- !is.na(j) & !is.na(oc$chr[j]) & is.finite(oc$pos[j]) &
    (coords$chr != oc$chr[j] | coords$pos != oc$pos[j])
  mismatch[is.na(mismatch)] <- FALSE
  mark(mismatch, "coordinate_mismatch")
  counts <- data.frame(n_saved = nrow(raw), n_duplicates_removed = dd$removed,
                       n_duplicate_conflicts = length(dd$conflicts),
                       n_excluded_target = sum(reason == "target", na.rm = TRUE),
                       n_excluded_region = sum(reason == "region", na.rm = TRUE),
                       n_excluded_mhc = sum(reason == "mhc", na.rm = TRUE),
                       n_excluded_coordinate_missing = sum(reason == "coordinate_missing", na.rm = TRUE),
                       n_excluded_coordinate_mismatch = sum(reason == "coordinate_mismatch", na.rm = TRUE),
                       n_retained = sum(is.na(reason)))
  excluded <- inst[!is.na(reason), , drop = FALSE]
  excluded$exclusion_reason <- reason[!is.na(reason)]
  kept <- inst[is.na(reason), , drop = FALSE]
  counts$retained_SNPs <- paste(kept$SNP, collapse = ";")
  counts$retained_SNPs_before_outcome_matching <- counts$retained_SNPs
  counts$n_retained_before_outcome_matching <- nrow(kept)
  counts$n_exclusion_loci <- nrow(loci)
  counts$excluded_locus_snps <- paste(loci$target_snp, collapse = ";")
  list(instruments = kept, excluded_instruments = excluded, counts = counts,
       duplicate_conflict_snps = dd$conflicts)
}

offline_as_exposure <- function(dat) {
  names(dat) <- sub("\\.outcome$", ".exposure", names(dat))
  names(dat)[names(dat) == "outcome"] <- "exposure"
  if ("chr" %in% names(dat)) names(dat)[names(dat) == "chr"] <- "chr.exposure"
  if ("pos" %in% names(dat)) names(dat)[names(dat) == "pos"] <- "pos.exposure"
  dat
}

offline_is_proxy <- function(dat, role) {
  proxy <- rep(FALSE, nrow(dat))
  field <- paste0("proxy.", role)
  if (field %in% names(dat)) {
    flag <- tolower(as.character(dat[[field]]))
    proxy <- proxy | (!is.na(flag) & flag %in% c("true", "t", "1", "yes"))
  }
  target <- paste0("target_snp.", role)
  if (target %in% names(dat))
    proxy <- proxy | (!is.na(dat[[target]]) & dat[[target]] != dat$SNP)
  proxy
}

harmonise_local <- function(exposure, outcome, expected_build = NULL, declared_build = NULL) {
  ex <- offline_deduplicate(exposure, "exposure")
  ou <- offline_deduplicate(outcome, "outcome")
  e <- ex$data
  o <- ou$data
  outcome_ids <- unique(o$id.outcome[!is.na(o$id.outcome)])
  if (length(outcome_ids) > 1L) stop("Harmonisation requires a single outcome GWAS ID")
  match_ix <- match(e$SNP, o$SNP)
  counts <- data.frame(n_missing_outcome = sum(is.na(match_ix)),
                       n_duplicate_conflicts_exposure = length(ex$conflicts),
                       n_duplicate_conflicts_outcome = length(ou$conflicts),
                       n_coordinate_mismatch = 0L, n_coordinate_missing = 0L, n_build_mismatch = 0L, n_palindromic = 0L,
                       n_incompatible = 0L, n_invalid_association = 0L,
                       n_proxy = 0L, n_usable = 0L)
  e <- e[!is.na(match_ix), , drop = FALSE]
  if (!nrow(e)) {
    result <- data.frame()
    attr(result, "counts") <- counts
    attr(result, "reason") <- if (length(ex$conflicts)) "Conflicting duplicate target/exposure associations" else if (length(ou$conflicts)) "Conflicting duplicate outcome associations" else "No matching outcome associations"
    return(result)
  }
  o <- o[o$SNP %in% e$SNP, , drop = FALSE]
  ec <- offline_coord(e, "exposure")
  oc <- offline_coord(o, "outcome")[match(e$SNP, o$SNP), ]
  mismatch <- !is.na(ec$chr) & !is.na(oc$chr) & is.finite(ec$pos) & is.finite(oc$pos) &
    (ec$chr != oc$chr | ec$pos != oc$pos)
  mismatch[is.na(mismatch)] <- FALSE
  missing_coord <- !offline_valid_coords(ec$chr, ec$pos) | !offline_valid_coords(oc$chr, oc$pos)
  eb <- offline_data_build(e, "exposure")
  ob <- offline_data_build(o, "outcome")[match(e$SNP, o$SNP)]
  build_mismatch <- !is.na(eb) & !is.na(ob) & eb != ob
  if (!is.null(expected_build)) {
    expected <- offline_build(expected_build)
    declared <- offline_build(declared_build)
    if (length(declared) == 1L && !is.na(declared)) {
      eb[is.na(eb)] <- declared
      ob[is.na(ob)] <- declared
    }
    if (length(expected) != 1L || is.na(expected) || !expected %in% c("GRCh37", "GRCh38"))
      build_mismatch <- rep(TRUE, nrow(e))
    else build_mismatch <- build_mismatch | is.na(eb) | is.na(ob) | eb != expected | ob != expected
  }
  counts$n_coordinate_mismatch <- sum(mismatch)
  counts$n_coordinate_missing <- sum(missing_coord)
  counts$n_build_mismatch <- sum(build_mismatch)
  proxy <- offline_is_proxy(e, "exposure")
  oi <- match(e$SNP, o$SNP)
  proxy <- proxy | offline_is_proxy(o, "outcome")[oi]
  counts$n_proxy <- sum(proxy)
  valid <- is.finite(e$beta.exposure) & is.finite(e$se.exposure) & e$se.exposure > 0 &
    e$beta.exposure != 0
  oi <- match(e$SNP, o$SNP)
  valid <- valid & is.finite(o$beta.outcome[oi]) & is.finite(o$se.outcome[oi]) & o$se.outcome[oi] > 0
  if ("mr_keep.exposure" %in% names(e)) valid <- valid & e$mr_keep.exposure %in% TRUE
  if ("mr_keep.outcome" %in% names(o)) valid <- valid & o$mr_keep.outcome[oi] %in% TRUE
  counts$n_invalid_association <- sum(!valid)
  e <- e[!mismatch & !missing_coord & !build_mismatch & !proxy & valid, , drop = FALSE]
  if (nrow(e)) {
    # Action 3 discards every palindromic A/T or C/G variant, even with EAF.
    result <- suppressMessages(TwoSampleMR::harmonise_data(e, o, action = 3))
    counts$n_palindromic <- sum(result$palindromic %in% TRUE)
    counts$n_incompatible <- sum(result$remove %in% TRUE & !result$palindromic %in% TRUE)
    counts$n_usable <- sum(result$mr_keep %in% TRUE)
  } else result <- data.frame()
  attr(result, "counts") <- counts
  attr(result, "reason") <- if (counts$n_usable) "" else paste(
    "No usable harmonised variants;", paste(names(counts)[counts[1, ] > 0], collapse = ", "))
  result
}

estimate_mr <- function(dat) {
  result <- data.frame(b = NA_real_, se = NA_real_, pval = NA_real_, lo = NA_real_,
                       hi = NA_real_, nsnp = 0L, method = NA_character_,
                       status = "unavailable", reason = "No usable instruments",
                       F_min = NA_real_, F_mean = NA_real_, n_weak_F = NA_integer_)
  if (!nrow(dat)) {
    if (!is.null(attr(dat, "reason"))) result$reason <- attr(dat, "reason")
    return(result)
  }
  d <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
  if (!nrow(d)) {
    result$reason <- attr(dat, "reason")
    return(result)
  }
  stopifnot(!anyDuplicated(d$SNP))
  strength <- (d$beta.exposure / d$se.exposure)^2
  result$F_min <- min(strength)
  result$F_mean <- mean(strength)
  result$n_weak_F <- sum(strength < 10)
  fit <- if (nrow(d) == 1L) {
    result$method <- "Wald ratio"
    TwoSampleMR::mr_wald_ratio(d$beta.exposure, d$beta.outcome,
                               d$se.exposure, d$se.outcome, TwoSampleMR::default_parameters())
  } else {
    result$method <- "IVW (multiplicative random effects)"
    TwoSampleMR::mr_ivw(d$beta.exposure, d$beta.outcome,
                        d$se.exposure, d$se.outcome)
  }
  result$nsnp <- nrow(d)
  if (!is.finite(fit$b) || !is.finite(fit$se) || fit$se <= 0) {
    result$reason <- "MR estimate or standard error is non-finite"
    return(result)
  }
  result$b <- fit$b
  result$se <- fit$se
  result$pval <- fit$pval
  result$lo <- fit$b - qnorm(0.975) * fit$se
  result$hi <- fit$b + qnorm(0.975) * fit$se
  result$status <- "available"
  result$reason <- ""
  result
}

compute_shared_comparators <- function(inputs, window_bp = 1e6, mhc = FALSE) {
  pairs <- canonical_pairs(inputs)
  mediator_ids <- unique(pairs$mediator_id)
  loci <- known_lead_loci(inputs)
  comparators <- setNames(lapply(mediator_ids, function(mediator_id) {
    targets <- unique(pairs$target_snp[pairs$mediator_id == mediator_id])
    selected <- tryCatch(select_instruments(inputs, mediator_id, targets,
                                            window_bp, mhc), error = identity)
    if (inherits(selected, "error")) {
      comparator <- estimate_mr(data.frame())
      comparator$reason <- conditionMessage(selected)
      counts <- data.frame(n_saved = sum(inputs$int_inst$id.exposure == mediator_id),
                           n_duplicates_removed = NA_integer_, n_duplicate_conflicts = NA_integer_,
                           n_excluded_target = NA_integer_, n_excluded_region = NA_integer_,
                           n_excluded_mhc = NA_integer_, n_excluded_coordinate_missing = NA_integer_,
                           n_excluded_coordinate_mismatch = NA_integer_, n_retained = NA_integer_,
                           retained_SNPs = "", retained_SNPs_before_outcome_matching = "",
                           n_retained_before_outcome_matching = NA_integer_, n_exclusion_loci = nrow(loci),
                           excluded_locus_snps = paste(loci$target_snp, collapse = ";"))
      hc <- attr(harmonise_local(inputs$int_inst[FALSE, ], inputs$outgwasf), "counts")
    } else {
      comp_h <- harmonise_local(selected$instruments, inputs$outgwasf)
      comparator <- estimate_mr(comp_h)
      counts <- selected$counts
      hc <- attr(comp_h, "counts")
    }
    if (inherits(selected, "error")) comp_h <- data.frame()
    usable <- if (nrow(comp_h)) comp_h$SNP[comp_h$mr_keep %in% TRUE] else character()
    counts$usable_SNPs <- paste(usable, collapse = ";")
    counts$retained_SNPs_after_harmonisation <- counts$usable_SNPs
    counts$n_retained_after_harmonisation <- length(usable)
    list(estimate = comparator, counts = counts, harmonisation_counts = hc,
         selection = if (inherits(selected, "error")) NULL else selected,
         window_bp = window_bp, mhc = mhc)
  }), mediator_ids)
  comparators
}

compute_comparisons <- function(inputs, window_bp = 1e6, mhc = FALSE, comparators = NULL) {
  pairs <- canonical_pairs(inputs)
  if (is.null(comparators)) comparators <- compute_shared_comparators(inputs, window_bp, mhc)
  for (id in unique(pairs$mediator_id)) {
    shared <- comparators[[id]]
    if (is.null(shared) || shared$window_bp != window_bp || shared$mhc != mhc)
      stop("Missing or incompatible shared comparator for ", id)
  }
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    exact <- inputs$int_chd[inputs$int_chd$id.outcome == pair$mediator_id &
                            inputs$int_chd$SNP == pair$target_snp, , drop = FALSE]
    # Original GBC PheWAS targets carried no coordinate columns. Retain that
    # compatibility using validated source lead coordinates; normalized inputs
    # must provide their own association coordinates and never use this fallback.
    if (is.null(inputs$lead_inventory) && !any(c("chr.outcome", "chr", "pos.outcome", "pos") %in% names(exact))) {
      exact$chr.outcome <- pair$target_chr
      exact$pos.outcome <- pair$target_pos
    }
    h <- harmonise_local(offline_as_exposure(exact), inputs$outgwasf,
      expected_build = if (is.null(inputs$lead_inventory)) NULL else offline_outcome_metadata(inputs)$build,
      declared_build = inputs$coordinate_build)
    target <- estimate_mr(h)
    # A total G->Y association remains interpretable even when G->M is ambiguous.
    raw_out <- inputs$outgwasf[inputs$outgwasf$SNP == pair$target_snp, , drop = FALSE]
    raw_dd <- offline_deduplicate(raw_out, "outcome")
    raw_out <- raw_dd$data
    raw_valid <- nrow(raw_out) == 1L && is.finite(raw_out$beta.outcome) &&
      is.finite(raw_out$se.outcome) && raw_out$se.outcome > 0 &&
      !offline_is_proxy(raw_out, "outcome") &&
      !is.na(raw_out$effect_allele.outcome) && nzchar(raw_out$effect_allele.outcome) &&
      !is.na(raw_out$other_allele.outcome) && nzchar(raw_out$other_allele.outcome) &&
      raw_out$effect_allele.outcome != raw_out$other_allele.outcome &&
      (!"mr_keep.outcome" %in% names(raw_out) || raw_out$mr_keep.outcome %in% TRUE)
    if (raw_valid) {
      rc <- offline_coord(raw_out, "outcome")
      raw_valid <- offline_valid_coords(rc$chr, rc$pos) &&
        !is.na(pair$target_chr) && is.finite(pair$target_pos) &&
        rc$chr == pair$target_chr && rc$pos == pair$target_pos
    }
    if (raw_valid && !is.null(inputs$lead_inventory)) {
      metadata <- offline_outcome_metadata(inputs)
      rb <- offline_data_build(raw_out, "outcome")
      if (is.na(rb)) rb <- offline_build(inputs$coordinate_build)
      expected <- offline_build(metadata$build)
      raw_valid <- length(rb) == 1L && length(expected) == 1L &&
        !is.na(rb) && !is.na(expected) && rb == expected
    }
    pair$total_association_status <- if (raw_valid) "available" else "unavailable"
    pair$total_association_reason <- if (raw_valid) "" else if (length(raw_dd$conflicts))
      "Conflicting duplicate outcome associations" else "Missing or invalid exact outcome association, alleles, coordinates or build"
    pair$effect_allele <- if (raw_valid) raw_out$effect_allele.outcome else NA_character_
    pair$a <- pair$a_se <- NA_real_
    pair$t <- if (raw_valid) raw_out$beta.outcome else NA_real_
    pair$t_se <- if (raw_valid) raw_out$se.outcome else NA_real_
    if (target$nsnp == 1L && target$status == "available") {
      stopifnot(nrow(h[h$mr_keep %in% TRUE, ]) == 1L,
                h$SNP[h$mr_keep %in% TRUE] == pair$target_snp)
      hh <- h[h$mr_keep %in% TRUE, ]
      pair$effect_allele <- hh$effect_allele.exposure
      pair$a <- hh$beta.exposure
      pair$a_se <- hh$se.exposure
      pair$t <- hh$beta.outcome
      pair$t_se <- hh$se.outcome
      pair$total_association_status <- "available"
      pair$total_association_reason <- ""
    }
    pair$target_status <- target$status
    pair$target_reason <- target$reason
    for (field in c("b", "se", "pval", "lo", "hi", "nsnp", "method", "F_min", "F_mean", "n_weak_F")) {
      pair[[paste0("target_", field)]] <- target[[field]]
    }
    shared <- comparators[[pair$mediator_id]]
    comparator <- shared$estimate
    counts <- shared$counts
    hc <- shared$harmonisation_counts
    for (field in c("b", "se", "pval", "lo", "hi", "nsnp", "method", "F_min", "F_mean", "n_weak_F", "status", "reason")) {
      pair[[paste0("comparator_", field)]] <- comparator[[field]]
    }
    pair$q <- pair$q_pval <- NA_real_
    pair$status <- if (target$status == "available" && comparator$status == "available") "available" else "unavailable"
    pair$reason <- paste(c(if (target$status != "available") paste("Target:", target$reason),
                          if (comparator$status != "available") paste("Comparator:", comparator$reason)), collapse = "; ")
    if (pair$status == "available") {
      pair$q <- (target$b - comparator$b)^2 / (target$se^2 + comparator$se^2)
      pair$q_pval <- pchisq(pair$q, df = 1, lower.tail = FALSE)
    }
    pair$window_bp <- window_bp
    pair$mhc_exclusion <- mhc
    cbind(pair, counts, hc)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  attr(result, "comparators") <- comparators
  result
}

comparison_sensitivity <- function(inputs) {
  wider <- compute_comparisons(inputs, window_bp = 2e6, mhc = FALSE)
  wider$sensitivity <- "+/-2 Mb"
  mhc <- compute_comparisons(inputs, window_bp = 2e6, mhc = TRUE)
  mhc$sensitivity <- "+/-2 Mb plus extended MHC (chr6:25-34 Mb)"
  rbind(wider, mhc)
}
