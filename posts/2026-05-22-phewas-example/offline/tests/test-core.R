# Run from the repository root with testthat::test_file().
core_dir <- if (dir.exists("posts/2026-05-22-phewas-example/offline/R")) {
  "posts/2026-05-22-phewas-example/offline/R"
} else if (dir.exists("../R")) "../R" else "R"
source(file.path(core_dir, "offline-core.R"))
source(file.path(core_dir, "offline-pooled.R"))
source(file.path(core_dir, "offline-inference.R"))

core_fixture <- function() {
  ex <- data.frame(id.exposure = "mediator-1", exposure = "Trait M",
    SNP = c("lead", "boundary", "outside", "background", "unpaired-lead",
            "near-unpaired", "unknown", "mismatch", "no-outcome"),
    chr.exposure = c(1, 1, 1, 3, 2, 2, NA, 3, 4),
    pos.exposure = c(10e6, 11e6, 12e6, 30e6, 20e6, 19e6, NA, 31e6, 40e6),
    beta.exposure = c(.2, .1, .3, .4, .5, .6, .7, .8, .9),
    se.exposure = .02, effect_allele.exposure = "A", other_allele.exposure = "G",
    eaf.exposure = .2, pval.exposure = 1e-9, mr_keep.exposure = TRUE,
    build.exposure = "GRCh37", stringsAsFactors = FALSE)
  ou <- ex[!ex$SNP %in% c("unpaired-lead", "no-outcome"), ]
  names(ou) <- sub("\\.exposure$", ".outcome", names(ou))
  names(ou)[names(ou) == "exposure"] <- "outcome"
  ou$id.outcome <- "outcome-z"
  ou$outcome <- "Example continuous outcome"
  ou$beta.outcome <- c(-.1, .05, .15, .19, .3, .35, .4)
  ou$se.outcome <- .04
  ou$chr.outcome[ou$SNP == "unknown"] <- 3
  ou$pos.outcome[ou$SNP == "unknown"] <- 35e6
  ou$pos.outcome[ou$SNP == "mismatch"] <- 31e6 + 1
  ou$effect_allele.outcome[1] <- "G"
  ou$other_allele.outcome[1] <- "A"
  ou$eaf.outcome[1] <- .8
  target <- ex[1, ]
  names(target) <- sub("\\.exposure$", ".outcome", names(target))
  names(target)[names(target) == "exposure"] <- "outcome"
  target$originalname.outcome <- "Trait M"
  target$units.outcome <- "mmol/L"
  list(int_chd = rbind(target, target), int_inst = rbind(ex, ex[3, ]),
       outgwasf = ou, phewas_result = data.frame(rsid = c("lead", "unpaired-lead")),
       lead_inventory = data.frame(target_snp = c("lead", "unpaired-lead"),
         target_chr = c(1, 2), target_pos = c(10e6, 20e6), build = "GRCh37"),
       coordinate_build = "GRCh37",
       outcome_metadata = list(name = "Example continuous outcome", id = "outcome-z",
                               units = "standard deviations", build = "GRCh37"))
}

testthat::test_that("full independent lead inventory excludes inclusive boundaries before matching", {
  x <- core_fixture()
  selected <- select_instruments(x, "mediator-1", "lead")
  testthat::expect_setequal(selected$instruments$SNP, c("outside", "background", "no-outcome"))
  testthat::expect_equal(selected$counts$n_excluded_target, 2L)
  testthat::expect_equal(selected$counts$n_excluded_region, 2L)
  testthat::expect_equal(selected$counts$n_excluded_coordinate_missing, 1L)
  testthat::expect_equal(selected$counts$n_excluded_coordinate_mismatch, 1L)
  testthat::expect_equal(selected$counts$n_duplicates_removed, 1L)
  testthat::expect_equal(selected$counts$n_exclusion_loci, 2L)
  testthat::expect_false("unpaired-lead" %in% x$outgwasf$SNP)
  testthat::expect_match(selected$counts$excluded_locus_snps, "unpaired-lead")
  # The inclusive region exclusion is established even if the outcome coordinate conflicts.
  x$outgwasf$pos.outcome[x$outgwasf$SNP == "boundary"] <- 100
  selected <- select_instruments(x, "mediator-1", "lead")
  testthat::expect_equal(selected$excluded_instruments$exclusion_reason[
    selected$excluded_instruments$SNP == "boundary"], "region")
})

testthat::test_that("missing/conflicting lead coordinates and builds fail explicitly", {
  x <- core_fixture()
  x$lead_inventory$target_pos[2] <- NA_real_
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "missing lead coordinates")
  bad <- compute_comparisons(x)
  testthat::expect_equal(bad$comparator_status, "unavailable")
  testthat::expect_match(bad$comparator_reason, "missing lead coordinates")
  x <- core_fixture()
  x$lead_inventory <- rbind(x$lead_inventory, transform(x$lead_inventory[1, ], target_pos = 10e6 + 1))
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "conflicting")
  x <- core_fixture()
  x$outcome_metadata$build <- "GRCh38"
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "Conflicting coordinate builds")
  x <- core_fixture()
  x$outcome_metadata$build <- NA_character_
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "build declaration")
})

testthat::test_that("exact identity, allele reversal and TwoSampleMR reference formulas are retained", {
  x <- core_fixture()
  r <- compute_comparisons(x)
  testthat::expect_equal(nrow(r), 1L)
  testthat::expect_equal(r$pair_id, "mediator-1::lead")
  testthat::expect_equal(r$target_nsnp, 1L)
  testthat::expect_equal(r$a, .2)
  testthat::expect_equal(r$t, .1)
  testthat::expect_equal(r$target_b, .5)
  testthat::expect_equal(r$effect_allele, "A")
  # These are the unchanged pinned-reference Wald and multiplicative-IVW fits.
  wald <- TwoSampleMR::mr_wald_ratio(.2, .1, .02, .04, TwoSampleMR::default_parameters())
  ivw <- TwoSampleMR::mr_ivw(c(.3, .4), c(.15, .19), c(.02, .02), c(.04, .04))
  testthat::expect_equal(r$target_se, wald$se)
  testthat::expect_equal(r$comparator_b, ivw$b)
  testthat::expect_equal(r$comparator_se, ivw$se)
  testthat::expect_equal(r$q, (.5 - ivw$b)^2 / (wald$se^2 + ivw$se^2))
  testthat::expect_equal(r$mr_effect_scale, "standard deviations per mmol/L")
  reversed <- x
  reversed$int_chd$beta.outcome <- -x$int_chd$beta.outcome
  reversed$int_chd$effect_allele.outcome <- "G"
  reversed$int_chd$other_allele.outcome <- "A"
  rev <- compute_comparisons(reversed)
  testthat::expect_equal(rev$a, -r$a)
  testthat::expect_equal(rev$t, -r$t)
  testthat::expect_equal(rev$target_b, r$target_b)
  testthat::expect_equal(rev$q, r$q)
  testthat::expect_equal(rev$target_se, r$target_se)
})

testthat::test_that("one shared fit is reused for every pair and pooled output", {
  x <- core_fixture()
  second <- x$int_chd[1, ]
  second$SNP <- "second-target"
  second$chr.outcome <- 5
  second$pos.outcome <- 50e6
  x$int_chd <- rbind(x$int_chd, second)
  x$lead_inventory <- rbind(x$lead_inventory,
    data.frame(target_snp = "second-target", target_chr = 5, target_pos = 50e6, build = "GRCh37"))
  y <- x$outgwasf[1, ]
  y$SNP <- second$SNP
  y$chr.outcome <- 5
  y$pos.outcome <- 50e6
  x$outgwasf <- rbind(x$outgwasf, y)
  comparisons <- annotate_comparisons(compute_comparisons(x))
  pooled <- compute_pooled(x, comparisons = comparisons)
  testthat::expect_equal(nrow(comparisons), 2L)
  testthat::expect_equal(length(unique(comparisons$comparator_b)), 1L)
  testthat::expect_equal(comparisons$comparator_b, rep(pooled$b, 2L))
  testthat::expect_equal(comparisons$comparator_se, rep(pooled$se, 2L))
  testthat::expect_equal(comparisons$comparator_nsnp, rep(pooled$nsnp, 2L))
  testthat::expect_equal(comparisons$usable_SNPs, rep(pooled$usable_SNPs, 2L))
  testthat::expect_equal(pooled$n_retained_before_outcome_matching, 3L)
  testthat::expect_equal(pooled$n_retained_after_harmonisation, 2L)
  testthat::expect_setequal(strsplit(pooled$usable_SNPs, ";", fixed = TRUE)[[1]], c("outside", "background"))
  testthat::expect_match(pooled$retained_SNPs_before_outcome_matching, "no-outcome")
  testthat::expect_error(compute_pooled(x, window_bp = 2e6, comparisons = comparisons), "settings")
  # A marker proves pooled output uses the supplied fit, rather than recomputing it.
  shared <- attr(comparisons, "comparators")
  shared[[1]]$estimate$b <- 123
  attr(comparisons, "comparators") <- shared
  testthat::expect_equal(compute_pooled(x, comparisons = comparisons)$b, 123)
})

testthat::test_that("duplicates, palindrome ambiguity and absent associations preserve inventory", {
  x <- core_fixture()
  x$int_chd$beta.outcome[2] <- .4
  r <- compute_comparisons(x)
  testthat::expect_equal(nrow(r), 1L)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_match(r$target_reason, "Conflicting duplicate")
  testthat::expect_true(is.na(r$a))
  testthat::expect_true(is.finite(r$t))
  x <- core_fixture()
  x$int_chd$other_allele.outcome <- "T"
  x$outgwasf$effect_allele.outcome[1] <- "A"
  x$outgwasf$other_allele.outcome[1] <- "T"
  r <- compute_comparisons(x)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_match(r$target_reason, "n_palindromic")
  testthat::expect_true(is.na(r$a))
  testthat::expect_true(is.finite(r$t))
  x <- core_fixture()
  x$outgwasf <- x$outgwasf[x$outgwasf$SNP != "lead", ]
  r <- compute_comparisons(x)
  testthat::expect_equal(nrow(r), 1L)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_true(is.na(r$t))
  testthat::expect_equal(r$comparator_status, "available")
  x <- core_fixture()
  x$int_inst$beta.exposure[nrow(x$int_inst)] <- .8
  selected <- select_instruments(x, "mediator-1", "lead")
  testthat::expect_false("outside" %in% selected$instruments$SNP)
  testthat::expect_equal(selected$duplicate_conflict_snps, "outside")
})

testthat::test_that("outcome proxies cannot supply exact or independent total effects", {
  x <- core_fixture()
  x$outgwasf$proxy.outcome <- x$outgwasf$SNP == "lead"
  r <- compute_comparisons(x)
  testthat::expect_match(r$target_reason, "n_proxy")
  testthat::expect_true(is.na(r$a) && is.na(r$t))
  x <- core_fixture()
  x$outgwasf$target_snp.outcome <- x$outgwasf$SNP
  x$outgwasf$target_snp.outcome[1] <- "different-exact-snp"
  r <- compute_comparisons(x)
  testthat::expect_true(is.na(r$t))
  testthat::expect_match(r$target_reason, "n_proxy")
  x <- core_fixture()
  x$int_chd$proxy.outcome <- TRUE
  r <- compute_comparisons(x)
  testthat::expect_true(is.na(r$a) && is.finite(r$t))
  testthat::expect_match(r$target_reason, "n_proxy")
})

testthat::test_that("BH uses available unique pairs and explicitly reports denominators", {
  x <- compute_comparisons(core_fixture())
  y <- x[rep(1, 3), ]
  y$target_snp <- c("one", "two", "three")
  y$pair_id <- paste(y$mediator_id, y$target_snp, sep = "::")
  y$q_pval <- c(.01, .04, NA_real_)
  result <- annotate_comparisons(y)
  testthat::expect_equal(result$q_padj, c(.02, .04, NA_real_))
  testthat::expect_equal(result$q_family_tested, rep(2L, 3))
  testthat::expect_equal(result$q_family_unavailable, rep(1L, 3))
  testthat::expect_error(annotate_comparisons(rbind(y, y[1, ])), "unique mediator/target")
  testthat::expect_equal(bh_available(c(.01, .04, -1, 2, Inf)), c(.02, .04, NA_real_, NA_real_, NA_real_))
})

testthat::test_that("legacy GBC fixture preserves pinned reference behavior", {
  x <- core_fixture()
  # The legacy inventory gets unpaired-lead coordinates from its outcome row.
  extra <- x$outgwasf[1, ]
  extra$SNP <- "unpaired-lead"
  extra$chr.outcome <- 2
  extra$pos.outcome <- 20e6
  x$outgwasf <- rbind(x$outgwasf, extra)
  x$outgwasf$id.outcome <- "gbc"
  x$outgwasf$outcome <- "GBC"
  x$lead_inventory <- x$outcome_metadata <- x$coordinate_build <- NULL
  r <- compute_comparisons(x)
  # Numerical golden values from pinned GBC scientific functions on this fixture.
  testthat::expect_equal(r$target_b, .5)
  testthat::expect_equal(r$target_se, .2)
  testthat::expect_equal(r$comparator_b, .484)
  testthat::expect_equal(r$comparator_se, .08)
  testthat::expect_equal(r$comparator_nsnp, 2L)
  testthat::expect_setequal(strsplit(r$usable_SNPs, ";", fixed = TRUE)[[1]], c("outside", "background"))
  testthat::expect_equal(known_gbc_loci(x), known_lead_loci(x))
})

# Compact real GBC associations from the pinned source bundle.
# Golden estimates were computed with the exact source commit.
gbc_reference_fixture <- function() {
list(int_chd = structure(list(SNP = "rs35604230", id.outcome = "finn-b-DM_RETINOPATHY",
    beta.outcome = -0.1745, se.outcome = 0.0271, effect_allele.outcome = "A",
    other_allele.outcome = "G", eaf.outcome = 0.06064, pval.outcome = 1.15601e-10,
    mr_keep.outcome = TRUE, originalname.outcome = "Diabetic retinopathy",
    outcome = "Diabetic retinopathy"), row.names = c(NA, -1L), class = "data.frame"),
    int_inst = structure(list(SNP = c("rs9275207", "rs73410774",
    "rs2854051", "rs10490924", "rs9277194", "rs6679677", "rs10922098",
    "rs3130151", "rs115380430", "rs34872471", "rs7260507", "rs689"
    ), id.exposure = c("finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY",
    "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY",
    "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY",
    "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY", "finn-b-DM_RETINOPATHY",
    "finn-b-DM_RETINOPATHY"), beta.exposure = c(0.5505, 0.6878,
    0.4692, 0.1802, 0.1144, 0.1357, -0.0849, -0.0816, 0.2537,
    0.094, -0.0763, 0.0888), se.exposure = c(0.0204, 0.0294,
    0.0371, 0.0151, 0.0137, 0.0181, 0.0129, 0.0134, 0.0431, 0.0161,
    0.0133, 0.016), effect_allele.exposure = c("G", "T", "T",
    "T", "T", "A", "T", "C", "C", "C", "C", "T"), other_allele.exposure = c("A",
    "G", "G", "G", "C", "C", "C", "T", "A", "T", "A", "A"), chr.exposure = c("6",
    "6", "6", "10", "6", "1", "1", "6", "6", "10", "19", "11"
    ), pos.exposure = c(32657710L, 32789739L, 32184050L, 124214448L,
    33023894L, 114303808L, 196664651L, 33115024L, 29849257L,
    114754071L, 41947625L, 2182224L), eaf.exposure = c(0.1206,
    0.05438, 0.03244, 0.2435, 0.4209, 0.1472, 0.563, 0.4346,
    0.02356, 0.2013, 0.3785, 0.7908), pval.exposure = c(5.5847e-160,
    1.22462e-120, 1.14499e-36, 6.85804e-33, 5.63508e-17, 7.45075e-14,
    5.19159e-11, 1.085e-09, 4.04203e-09, 4.81305e-09, 1.04501e-08,
    2.62398e-08), mr_keep.exposure = c(TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), exposure = c("Diabetic retinopathy",
    "Diabetic retinopathy", "Diabetic retinopathy", "Diabetic retinopathy",
    "Diabetic retinopathy", "Diabetic retinopathy", "Diabetic retinopathy",
    "Diabetic retinopathy", "Diabetic retinopathy", "Diabetic retinopathy",
    "Diabetic retinopathy", "Diabetic retinopathy")), row.names = c(NA,
    -12L), class = "data.frame"), outgwasf = structure(list(SNP = c("rs60541816",
    "rs67274221", "rs6733452", "rs15734", "rs2312959", "rs35604230",
    "rs77155094", "rs2131242", "rs2854051", "rs7260507", "rs9275207",
    "rs3130151", "rs6679677", "rs689", "rs9277194", "rs10922098",
    "rs34872471", "rs10490924"), id.outcome = c("f7Cgzv", "f7Cgzv",
    "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv",
    "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv",
    "f7Cgzv", "f7Cgzv", "f7Cgzv", "f7Cgzv"), beta.outcome = c(-0.543693227391608,
    0.269721747074788, 0.84070474491756, 0.266432899412607, 0.394538967974611,
    0.459195430529132, 0.4218632539178, -0.347998408459836, 0.116003675756306,
    0.0467881614987589, -0.0988263543068949, 0.0448777587780855,
    -0.142024500029418, 0.0398935636616766, 0.0256677467485778,
    -0.0259333820265044, 0.0073727543294131, -0.001300845733048
    ), se.outcome = c(0.0655456853925391, 0.0407018115653192,
    0.13069834073862, 0.0439050492284711, 0.0686001952409779,
    0.0808152860449297, 0.077044596492064, 0.0691326132463761,
    0.0749533979909975, 0.0396847695482637, 0.0851502588684111,
    0.0402448637931589, 0.19090312266697, 0.0543571327825631,
    0.0404384010314384, 0.0430277857559159, 0.0436926052225217,
    0.0406960000830267), effect_allele.outcome = c("A", "A",
    "A", "A", "A", "A", "T", "A", "T", "C", "G", "T", "A", "A",
    "T", "C", "C", "T"), other_allele.outcome = c("G", "G", "G",
    "G", "C", "G", "C", "G", "G", "A", "A", "C", "C", "T", "C",
    "T", "T", "G"), chr.outcome = c(7L, 3L, 2L, 3L, 16L, 6L,
    7L, 11L, 6L, 19L, 6L, 6L, 1L, 11L, 6L, 1L, 10L, 10L), pos.outcome = c(87127861L,
    56888845L, 44094845L, 8615489L, 8763255L, 31768305L, 107474366L,
    69839305L, 32184050L, 41947625L, 32657710L, 33115024L, 114303808L,
    2182224L, 33023894L, 196664651L, 114754071L, 124214448L),
        eaf.outcome = c(0.1291, 0.3733, 0.01892, 0.2626, 0.07619,
        0.0831, 0.05823, 0.1045, 0.0817, 0.4352, 0.06816, 0.4582,
        0.0114, 0.1495, 0.4208, 0.3004, 0.2924, 0.3532), pval.outcome = c(1e-200,
        4.685e-11, 9.069e-11, 1.233e-09, 9.324e-11, 1.152e-08,
        2.14e-08, 9.862e-07, 0.08987, 0.2077, 0.4169, 0.2481,
        0.3094, 0.3787, 0.5034, 0.6144, 0.907, 0.9777), mr_keep.outcome = c(TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), outcome = c("GBC",
        "GBC", "GBC", "GBC", "GBC", "GBC", "GBC", "GBC", "GBC",
        "GBC", "GBC", "GBC", "GBC", "GBC", "GBC", "GBC", "GBC",
        "GBC")), row.names = c(NA, -18L), class = "data.frame"),
    phewas_result = structure(list(rsid = c("rs15734", "rs2131242",
    "rs2312959", "rs35604230", "rs60541816", "rs67274221", "rs6733452",
    "rs77155094")), class = "data.frame", row.names = c(NA, -8L
    )))
}


testthat::test_that("real GBC saved associations match the exact pinned implementation", {
  x <- gbc_reference_fixture()
  r <- compute_comparisons(x)
  testthat::expect_equal(r$target_snp, "rs35604230")
  testthat::expect_equal(r$target_b, -2.63149243856236, tolerance = 1e-12)
  testthat::expect_equal(r$target_se, 0.463124848395013, tolerance = 1e-12)
  testthat::expect_equal(r$comparator_b, 0.00659895577204082, tolerance = 1e-12)
  testthat::expect_equal(r$comparator_se, 0.149999472465065, tolerance = 1e-12)
  testthat::expect_equal(r$comparator_nsnp, 7L)
  testthat::expect_equal(r$status, "available")
  testthat::expect_equal(r$retained_SNPs,
    "rs73410774;rs10490924;rs9277194;rs6679677;rs10922098;rs3130151;rs115380430;rs34872471;rs7260507;rs689")
  p <- compute_pooled(x, comparisons = r)
  testthat::expect_equal(p$b, r$comparator_b)
  testthat::expect_equal(p$se, r$comparator_se)
  testthat::expect_equal(p$nsnp, r$comparator_nsnp)
})

testthat::test_that("total availability survives missing or mismatched mediator coordinates", {
  x <- core_fixture()
  x$int_chd$pos.outcome <- 10e6 + 1
  r <- compute_comparisons(x)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_equal(r$comparator_status, "unavailable")
  testthat::expect_match(r$target_reason, "coordinate_mismatch")
  testthat::expect_equal(r$total_association_status, "available")
  testthat::expect_true(is.finite(r$t))
  x <- core_fixture()
  x$int_chd$pos.outcome <- NA_real_
  r <- compute_comparisons(x)
  testthat::expect_match(r$target_reason, "coordinate_missing")
  testthat::expect_true(is.na(r$a) && is.finite(r$t))
  x <- core_fixture()
  x$outgwasf$build.outcome[1] <- "GRCh38"
  r <- compute_comparisons(x)
  testthat::expect_true(is.na(r$t))
  testthat::expect_match(r$target_reason, "n_build_mismatch")
})

testthat::test_that("unsupported coordinates are rejected conservatively", {
  x <- core_fixture()
  x$lead_inventory$target_chr[2] <- "unknown-contig"
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "missing lead coordinates")
  x <- core_fixture()
  x$int_inst$pos.exposure[x$int_inst$SNP == "background"] <- 30e6 + .5
  selected <- select_instruments(x, "mediator-1", "lead")
  testthat::expect_false("background" %in% selected$instruments$SNP)
  x <- core_fixture()
  x$outgwasf$pos.outcome[x$outgwasf$SNP == "outside"] <- NA_real_
  p <- compute_pooled(x)
  testthat::expect_equal(p$n_retained_after_harmonisation, 1L)
  testthat::expect_match(p$usable_SNPs, "^background$")
})

testthat::test_that("invalid association fields stay unavailable without dropping pairs", {
  x <- core_fixture()
  x$int_chd$se.outcome <- 0
  r <- compute_comparisons(x)
  testthat::expect_equal(nrow(r), 1L)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_match(r$target_reason, "invalid_association")
  testthat::expect_true(is.finite(r$t))
  x <- core_fixture()
  x$outgwasf$se.outcome[1] <- NA_real_
  r <- compute_comparisons(x)
  testthat::expect_true(is.na(r$t))
  testthat::expect_equal(r$total_association_status, "unavailable")
  x <- core_fixture()
  x$coordinate_build <- NULL
  x$int_inst$build.exposure <- NULL
  testthat::expect_error(select_instruments(x, "mediator-1", "lead"), "Missing association coordinate build")
})

testthat::test_that("agreement between associations cannot override the declared build", {
  x <- core_fixture()
  x$int_chd$build.outcome <- "GRCh38"
  x$outgwasf$build.outcome[1] <- "GRCh38"
  r <- compute_comparisons(x)
  testthat::expect_equal(r$target_status, "unavailable")
  testthat::expect_match(r$target_reason, "n_build_mismatch")
  testthat::expect_true(is.na(r$a) && is.na(r$t))
})
