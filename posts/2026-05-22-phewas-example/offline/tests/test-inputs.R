source(testthat::test_path("../R/chd-inputs.R"))

testthat::test_that("VCF LP values convert to valid P values and preserve missingness", {
  testthat::expect_equal(chd_lp_to_p(c(0, 1, 2)), c(1, 0.1, 0.01))
  testthat::expect_true(is.na(chd_lp_to_p(NA_real_)))
  testthat::expect_error(chd_lp_to_p(-0.1), "non-negative")
  testthat::expect_error(chd_lp_to_p(Inf), "finite")
  testthat::expect_error(chd_lp_to_p(NaN), "NaN")
  testthat::expect_equal(chd_lp_to_p(1000), 0)
  testthat::expect_true(is.na(chd_lp_to_p(NA_real_)))
})

testthat::test_that("lead coordinate inventories reject conflicting or incomplete sources", {
  clean <- data.frame(id = c("a", "b"), chr = c("1", "1"), pos = c(10, 10),
                      snp = c("rs1", "rs1"))
  got <- chd_unique_locus_table(clean, "chr", "pos", "snp")
  testthat::expect_equal(nrow(got), 1L)
  testthat::expect_equal(got$target_pos, 10)

  conflict <- rbind(clean[1, ], transform(clean[2, ], pos = 11))
  testthat::expect_error(chd_unique_locus_table(conflict, "chr", "pos", "snp"),
                         "Conflicting saved coordinates")
})

testthat::test_that("association schema preserves duplicate records and carries build", {
  rows <- data.frame(SNP = c("rs1", "rs1"), effect_allele = c("A", "A"),
                     other_allele = c("G", "G"), beta = c(0.2, 0.3), se = c(0.1, 0.1),
                     pval = c(0.04, 0.003), eaf = c(0.2, 0.2), samplesize = c(100, 100),
                     chr = c("1", "1"), pos = c(10, 10), mr_keep = c(TRUE, TRUE))
  got <- chd_assoc_schema(rows, "id1", "trait", "exposure", "units")
  testthat::expect_equal(nrow(got), 2L)
  testthat::expect_equal(got$build.exposure, c("GRCh37", "GRCh37"))
  testthat::expect_equal(got$beta.exposure, c(0.2, 0.3))
})

testthat::test_that("the saved CHD input file is the intended result artifact", {
  saved <- new.env(parent = emptyenv())
  load("../../chd_phewas_result.Rdata", envir = saved)
  testthat::expect_equal(nrow(saved$chd), 41L)
  testthat::expect_equal(nrow(saved$int_chd), 23L)
  testthat::expect_equal(nrow(saved$int_inst), 1002L)
  testthat::expect_false(identical(tools::md5sum("../../chd_phewas_result.Rdata"),
                                   tools::md5sum("../../chd_phewas_results.Rdata")))
})
