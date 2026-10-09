# Local, reproducible input reconstruction for the CHD offline analysis.
# No network access or package installation is performed here.

chd_mediator_files <- c(
  "ieu-b-5138" = "ieu-b-5138.vcf.gz",
  "ukb-d-30780_raw" = "ukb-d-30780_raw.vcf.gz",
  "ukb-b-10831" = "ukb-b-10831.vcf.gz",
  "ukb-b-5779" = "ukb-b-5779.vcf.gz",
  "ukb-b-2303" = "ukb-b-2303.vcf.gz"
)

chd_mediator_labels <- c("ieu-b-5138" = "SBP", "ukb-d-30780_raw" = "LDL",
                         "ukb-b-10831" = "Smoking", "ukb-b-5779" = "Alcohol",
                         "ukb-b-2303" = "BMI")

chd_lp_to_p <- function(lp) {
  lp <- as.numeric(lp)
  if (any(is.nan(lp))) stop("VCF LP values cannot be NaN", call. = FALSE)
  present <- !is.na(lp)
  if (any(!is.finite(lp[present])) || any(lp[present] < 0)) {
    stop("VCF LP values must be finite, non-negative, or missing", call. = FALSE)
  }
  p <- rep(NA_real_, length(lp))
  p[present] <- 10^(-lp[present])
  if (any(!is.finite(p[present])) || any(p[present] < 0 | p[present] > 1)) {
    stop("Converted VCF P values must be finite and in [0, 1]", call. = FALSE)
  }
  p
}

chd_require_packages <- function() {
  required <- c("gwasvcf", "VariantAnnotation", "Rsamtools", "GenomicRanges")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing local R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

chd_vcf_header <- function(path) {
  tabix <- Rsamtools::TabixFile(path)
  header <- VariantAnnotation::scanVcfHeader(tabix)
  seqinfo <- VariantAnnotation::seqinfo(header)
  text <- paste(capture.output(show(header)), capture.output(show(seqinfo)), collapse = "\n")
  genome <- paste(unique(stats::na.omit(methods::slot(seqinfo, "genome"))), collapse = ";")
  list(header = text, build = if (grepl("HG19/GRCh37|GRCh37|b37", genome, ignore.case = TRUE)) {
    "GRCh37"
  } else if (grepl("GRCh38|hg38", text, ignore.case = TRUE)) {
    "GRCh38"
  } else NA_character_)
}

chd_query_vcf <- function(path, loci, build, role, id) {
  if (all(c("target_chr", "target_pos") %in% names(loci))) {
    names(loci)[match(c("target_chr", "target_pos"), names(loci))] <- c("chr", "pos")
  }
  loci <- unique(loci[is.finite(loci$pos) & !is.na(loci$chr) & nzchar(loci$chr), c("chr", "pos")])
  if (!nrow(loci)) return(data.frame())
  ranges <- data.frame(chrom = as.character(loci$chr), start = as.integer(loci$pos),
                       end = as.integer(loci$pos), stringsAsFactors = FALSE)
  variants <- gwasvcf::query_gwas(path, chrompos = ranges, build = build)
  if (!length(variants)) return(data.frame())
  dat <- as.data.frame(gwasvcf::vcf_to_tibble(variants))
  dat$chr <- sub("^chr", "", as.character(dat$seqnames), ignore.case = TRUE)
  dat$pos <- as.numeric(dat$start)
  dat$SNP <- as.character(dat$rsid)
  dat$effect_allele <- as.character(dat$ALT)
  dat$other_allele <- as.character(dat$REF)
  dat$beta <- as.numeric(dat$ES)
  dat$se <- as.numeric(dat$SE)
  dat$pval <- chd_lp_to_p(dat$LP)
  dat$eaf <- as.numeric(dat$AF)
  dat$samplesize <- as.numeric(dat$SS)
  dat$study_variant_id <- as.character(dat$ID)
  dat$build <- build
  dat$source_file <- basename(path)
  dat$role <- role
  dat$gwas_id <- id
  dat
}

chd_unique_locus_table <- function(dat, chr, pos, snp) {
  x <- data.frame(target_snp = as.character(dat[[snp]]),
                  target_chr = sub("^chr", "", as.character(dat[[chr]]), ignore.case = TRUE),
                  target_pos = suppressWarnings(as.numeric(dat[[pos]])),
                  stringsAsFactors = FALSE)
  x <- x[!is.na(x$target_snp) & nzchar(x$target_snp), , drop = FALSE]
  groups <- split(seq_len(nrow(x)), x$target_snp)
  out <- lapply(groups, function(ix) {
    coords <- unique(x[ix, c("target_chr", "target_pos"), drop = FALSE])
    coords <- coords[!is.na(coords$target_chr) & nzchar(coords$target_chr) &
                       is.finite(coords$target_pos), , drop = FALSE]
    if (nrow(coords) > 1L) stop("Conflicting saved coordinates for ", x$target_snp[ix[1]], call. = FALSE)
    data.frame(target_snp = x$target_snp[ix[1]],
               target_chr = if (nrow(coords)) coords$target_chr[1] else NA_character_,
               target_pos = if (nrow(coords)) coords$target_pos[1] else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

chd_assoc_schema <- function(associations, id, label, role = c("exposure", "outcome"), units = NA_character_) {
  role <- match.arg(role)
  suffix <- paste0(".", role)
  n <- nrow(associations)
  data.frame(
    SNP = as.character(associations$SNP),
    setNames(list(rep(as.character(id), n)), paste0("id", suffix)),
    setNames(list(rep(as.character(label), n)), role),
    setNames(list(as.character(associations$effect_allele)), paste0("effect_allele", suffix)),
    setNames(list(as.character(associations$other_allele)), paste0("other_allele", suffix)),
    setNames(list(as.numeric(associations$beta)), paste0("beta", suffix)),
    setNames(list(as.numeric(associations$se)), paste0("se", suffix)),
    setNames(list(as.numeric(associations$pval)), paste0("pval", suffix)),
    setNames(list(as.numeric(associations$eaf)), paste0("eaf", suffix)),
    setNames(list(as.numeric(associations$samplesize)), paste0("samplesize", suffix)),
    setNames(list(as.character(associations$chr)), paste0("chr", suffix)),
    setNames(list(as.numeric(associations$pos)), paste0("pos", suffix)),
    setNames(list(rep("GRCh37", n)), paste0("build", suffix)),
    setNames(list(rep(units, n)), paste0("units", suffix)),
    setNames(list(if ("mr_keep" %in% names(associations)) as.logical(associations$mr_keep) else rep(TRUE, n)), paste0("mr_keep", suffix)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

reconstruct_chd_inputs <- function(input_dir = ".") {
  chd_require_packages()
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  rdata <- file.path(input_dir, "chd_phewas_result.Rdata")
  if (!file.exists(rdata)) stop("Missing CHD input Rdata: ", rdata, call. = FALSE)
  saved <- new.env(parent = emptyenv())
  load(rdata, envir = saved)
  required <- c("chd", "int_chd", "int_inst", "phewas_result")
  missing <- setdiff(required, ls(saved))
  if (length(missing)) stop("Rdata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  leads <- as.data.frame(saved$chd)
  pairs <- as.data.frame(saved$int_chd)
  instruments <- as.data.frame(saved$int_inst)
  historical_phewas <- as.data.frame(saved$phewas_result)
  if (!all(c("rsid", "chr", "position") %in% names(leads))) stop("CHD lead table lacks rsid/chr/position", call. = FALSE)
  if (!all(c("SNP", "id.outcome") %in% names(pairs))) stop("Saved pair identities lack SNP/id.outcome", call. = FALSE)
  if (!all(c("SNP", "id.exposure", "chr.exposure", "pos.exposure") %in% names(instruments))) stop("Saved instrument identities lack SNP/id/coordinates", call. = FALSE)

  lead_inventory <- chd_unique_locus_table(leads, "chr", "position", "rsid")
  lead_inventory$build <- "GRCh37"
  if (anyNA(lead_inventory$target_chr) || any(!is.finite(lead_inventory$target_pos))) {
    stop("The 41-lead coordinate inventory contains missing coordinates", call. = FALSE)
  }
  instrument_loci <- chd_unique_locus_table(instruments, "chr.exposure", "pos.exposure", "SNP")
  pair_coords <- lead_inventory[match(pairs$SNP, lead_inventory$target_snp), , drop = FALSE]
  if (anyNA(pair_coords$target_chr) || any(!is.finite(pair_coords$target_pos))) {
    stop("A retained mediator target lacks coordinates in the saved CHD lead table", call. = FALSE)
  }
  outcome_loci <- rbind(
    data.frame(target_snp = instrument_loci$target_snp, target_chr = instrument_loci$target_chr,
               target_pos = instrument_loci$target_pos),
    data.frame(target_snp = lead_inventory$target_snp, target_chr = lead_inventory$target_chr,
               target_pos = lead_inventory$target_pos),
    data.frame(target_snp = pairs$SNP, target_chr = pair_coords$target_chr,
               target_pos = pair_coords$target_pos)
  )
  outcome_loci <- chd_unique_locus_table(outcome_loci, "target_chr", "target_pos", "target_snp")
  if (anyNA(outcome_loci$target_chr) || any(!is.finite(outcome_loci$target_pos))) stop("Outcome query inventory has missing coordinates", call. = FALSE)

  mediator_tables <- list()
  headers <- list()
  for (id in names(chd_mediator_files)) {
    path <- file.path(input_dir, chd_mediator_files[[id]])
    if (!file.exists(path) || !file.exists(paste0(path, ".tbi"))) stop("Missing indexed mediator VCF: ", path, call. = FALSE)
    headers[[id]] <- chd_vcf_header(path)
    if (!identical(headers[[id]]$build, "GRCh37")) stop("Cannot establish GRCh37 build from ", basename(path), call. = FALSE)
    mediator_tables[[id]] <- chd_query_vcf(path, pair_coords[, c("target_chr", "target_pos")],
                                           "GRCh37", "mediator", id)
  }
  outcome_path <- file.path(input_dir, "ukb-d-I9_IHD.vcf.gz")
  if (!file.exists(outcome_path) || !file.exists(paste0(outcome_path, ".tbi"))) stop("Missing indexed CHD outcome VCF", call. = FALSE)
  headers[["ukb-d-I9_IHD"]] <- chd_vcf_header(outcome_path)
  if (!identical(headers[["ukb-d-I9_IHD"]]$build, "GRCh37")) stop("Cannot establish GRCh37 build from CHD outcome VCF", call. = FALSE)

  # Match by saved ID and position. Keep every matching VCF row so the shared
  # analysis layer can flag conflicting or repeated association records.
  int_chd_rows <- list()
  for (i in seq_len(nrow(pairs))) {
    id <- as.character(pairs$id.outcome[i])
    vcf <- mediator_tables[[id]]
    if (is.null(vcf)) stop("No local mediator VCF mapped for ", id, call. = FALSE)
    hit <- vcf[vcf$SNP == pairs$SNP[i] & vcf$chr == pair_coords$target_chr[i] &
                 vcf$pos == pair_coords$target_pos[i], , drop = FALSE]
    if (!nrow(hit)) hit <- data.frame(SNP = pairs$SNP[i], chr = pair_coords$target_chr[i],
                                       pos = pair_coords$target_pos[i], effect_allele = NA_character_,
                                       other_allele = NA_character_, beta = NA_real_, se = NA_real_,
                                       pval = NA_real_, eaf = NA_real_, samplesize = NA_real_,
                                       build = "GRCh37")
    hit$id.outcome <- id
    hit$outcome <- as.character(instruments$exposure[match(id, instruments$id.exposure)])
    hit$source_trait.outcome <- hit$outcome
    hit$originalname.outcome <- unname(chd_mediator_labels[id])
    hit$mr_keep.outcome <- is.finite(hit$beta) & is.finite(hit$se) & hit$se > 0
    hit$mr_keep <- hit$mr_keep.outcome
    int_chd_rows[[i]] <- chd_assoc_schema(hit, id, hit$outcome[1], "outcome", "Saved mediator GWAS units; specific transformation/unit unavailable")
    int_chd_rows[[i]]$originalname.outcome <- hit$originalname.outcome[1]
    int_chd_rows[[i]]$source_trait.outcome <- hit$source_trait.outcome[1]
    int_chd_rows[[i]]$vcf_source <- chd_mediator_files[[id]]
  }
  int_chd <- do.call(rbind, int_chd_rows)
  rownames(int_chd) <- NULL

  # Reconstruct each saved instrument identity from its local mediator VCF.
  int_inst_rows <- vector("list", nrow(instruments))
  for (id in names(chd_mediator_files)) {
    idx <- which(instruments$id.exposure == id)
    if (!length(idx)) next
    inst_vcf <- chd_query_vcf(file.path(input_dir, chd_mediator_files[[id]]),
                              instrument_loci[instrument_loci$target_snp %in% instruments$SNP[idx],
                                              c("target_chr", "target_pos"), drop = FALSE],
                              "GRCh37", "exposure", id)
    mapped <- lapply(idx, function(i) {
      hit <- inst_vcf[inst_vcf$SNP == instruments$SNP[i] & inst_vcf$chr == instruments$chr.exposure[i] &
                        inst_vcf$pos == instruments$pos.exposure[i], , drop = FALSE]
      if (!nrow(hit)) hit <- data.frame(SNP = instruments$SNP[i], chr = instruments$chr.exposure[i],
                                         pos = instruments$pos.exposure[i], effect_allele = NA_character_,
                                         other_allele = NA_character_, beta = NA_real_, se = NA_real_,
                                         pval = NA_real_, eaf = NA_real_, samplesize = NA_real_, build = "GRCh37")
      hit$id.exposure <- id
      hit$exposure <- instruments$exposure[i]
      hit$mr_keep.exposure <- is.finite(hit$beta) & is.finite(hit$se) & hit$se > 0
      hit$mr_keep <- hit$mr_keep.exposure
      row <- chd_assoc_schema(hit, id, hit$exposure[1], "exposure", "Saved mediator GWAS units; specific transformation/unit unavailable")
      row$vcf_source <- chd_mediator_files[[id]]
      row
    })
    int_inst_rows[idx] <- mapped
  }
  int_inst <- do.call(rbind, int_inst_rows)
  rownames(int_inst) <- NULL

  out_vcf <- chd_query_vcf(outcome_path, outcome_loci[, c("target_chr", "target_pos")],
                           "GRCh37", "outcome", "ukb-d-I9_IHD")
  out_vcf <- out_vcf[out_vcf$SNP %in% outcome_loci$target_snp, , drop = FALSE]
  outcome_expected <- outcome_loci[match(out_vcf$SNP, outcome_loci$target_snp), , drop = FALSE]
  mismatch <- out_vcf$chr != outcome_expected$target_chr | out_vcf$pos != outcome_expected$target_pos
  if (any(mismatch)) stop("Outcome VCF coordinate conflicts with saved coordinates for ",
                          paste(unique(out_vcf$SNP[mismatch]), collapse = ", "), call. = FALSE)
  outgwasf <- chd_assoc_schema(out_vcf, "ukb-d-I9_IHD", "CHD", "outcome",
                               "Saved GWAS units; scale unavailable")
  outgwasf$chr.outcome <- out_vcf$chr
  outgwasf$pos.outcome <- out_vcf$pos
  outgwasf$build.outcome <- "GRCh37"
  outgwasf$vcf_source <- basename(outcome_path)
  outcome_metadata <- list(name = "CHD", id = "ukb-d-I9_IHD",
                           units = "Saved GWAS units; scale unavailable", build = "GRCh37",
                           study_type = "Continuous", vcf_header = headers[["ukb-d-I9_IHD"]]$header)

  # Validate coordinates/build against every available association row.
  if (any(int_inst$build.exposure != "GRCh37", na.rm = TRUE) ||
      any(int_chd$build.outcome != "GRCh37", na.rm = TRUE) ||
      any(outgwasf$build.outcome != "GRCh37", na.rm = TRUE)) stop("Association build mismatch", call. = FALSE)
  saved_inst_keys <- unique(paste(instruments$id.exposure, instruments$SNP,
                                  instruments$chr.exposure, instruments$pos.exposure, sep = "\r"))
  rebuilt_inst_keys <- unique(paste(int_inst$id.exposure, int_inst$SNP,
                                    int_inst$chr.exposure, int_inst$pos.exposure, sep = "\r"))
  if (!setequal(saved_inst_keys, rebuilt_inst_keys)) stop("Reconstructed instrument identities differ from the saved inventory", call. = FALSE)

  vcf_paths <- c(file.path(input_dir, unname(chd_mediator_files)), outcome_path)
  files <- c(rdata, vcf_paths, paste0(vcf_paths, ".tbi"))
  provenance <- data.frame(path = basename(files), md5 = unname(tools::md5sum(files)),
                           size_bytes = as.numeric(file.info(files)$size), stringsAsFactors = FALSE)
  provenance$build <- c(NA_character_, rep("GRCh37", length(vcf_paths)), rep("GRCh37", length(vcf_paths)))
  provenance$header <- c(NA_character_, vapply(headers[names(chd_mediator_files)], `[[`, character(1), "header"),
                         headers[["ukb-d-I9_IHD"]]$header, rep(NA_character_, length(vcf_paths)))
  provenance$index_type <- c(NA_character_, rep("Tabix .tbi; indexed coordinate query via gwasvcf/Rsamtools", length(vcf_paths)),
                             rep("Tabix index paired with source VCF", length(vcf_paths)))
  extraction_counts <- data.frame(
    source = c("saved_pair_identities", "exact_mediator_associations", "saved_instrument_identities",
               "instrument_vcf_associations", "lead_inventory", "outcome_query_loci", "outcome_vcf_associations",
               "outcome_leads_found", "outcome_instruments_found", "outcome_target_snps_found"),
    queried = c(nrow(pairs), nrow(pairs), nrow(instruments), nrow(instruments), nrow(lead_inventory),
                nrow(outcome_loci), nrow(outcome_loci), nrow(lead_inventory),
                length(unique(instruments$SNP)), length(unique(pairs$SNP))),
    found = c(nrow(pairs), sum(is.finite(int_chd$beta.outcome)), nrow(instruments),
              sum(is.finite(int_inst$beta.exposure)), nrow(lead_inventory), nrow(outcome_loci), nrow(outgwasf),
              sum(lead_inventory$target_snp %in% outgwasf$SNP), sum(unique(instruments$SNP) %in% outgwasf$SNP),
              sum(unique(pairs$SNP) %in% outgwasf$SNP)),
    stringsAsFactors = FALSE
  )
  list(int_chd = int_chd, int_inst = int_inst, outgwasf = outgwasf,
       phewas_result = historical_phewas, lead_inventory = lead_inventory,
       outcome_metadata = outcome_metadata, coordinate_build = "GRCh37",
       provenance = provenance, extraction_counts = extraction_counts,
       inventory_counts = data.frame(n_saved_pair_rows = nrow(pairs),
                                     n_saved_instrument_rows = nrow(instruments),
                                     n_full_outcome_rows = NA_integer_,
                                     n_retained_outcome_rows = nrow(outgwasf),
                                     n_original_phewas_rows = nrow(historical_phewas)))
}

# `reconstruct_chd_inputs()` is the fresh-input path used by reports. The
# snapshot loader is for an explicit cached handoff and supports `refresh=TRUE`.
load_chd_inputs <- function(input_dir = ".", snapshot_path = "offline/data/chd-inputs.rds",
                            refresh = FALSE) {
  path <- if (grepl("^(/|~)", snapshot_path)) path.expand(snapshot_path) else file.path(input_dir, snapshot_path)
  if (!refresh && file.exists(path)) return(readRDS(path))
  inputs <- reconstruct_chd_inputs(input_dir)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(inputs, path, compress = "xz")
  inputs
}
