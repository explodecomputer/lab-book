testthat::test_that("source-derived signed components, availability and covariance scenarios agree", {
fixture <- data.frame(mediator_id='m1', target_snp='rs1', mediator='Mediator',
                      mediator_scale='SD', effect_allele='A', a=0.2, a_se=0.02,
                      t=0.1, t_se=0.01, comparator_b=0.5, comparator_se=0.05,
                      status='complete', reason=NA_character_)
result <- compute_decomposition(fixture)
stopifnot(nrow(result) == 1L, result$decomposition_status == 'complete',
          isTRUE(all.equal(result$i, 0.1)), isTRUE(all.equal(result$d, 0)),
          isTRUE(all.equal(result$t, result$i + result$d)),
          isTRUE(all.equal(result$i_se, sqrt(0.5^2 * 0.02^2 + 0.2^2 * 0.05^2))))
cat('Complete mediation and uncertainty: OK\n')
residual <- fixture
residual$t <- 0.16
opposing <- fixture
opposing$comparator_b <- -0.5
reversed <- fixture
reversed$a <- -fixture$a
reversed$t <- -fixture$t
signed <- compute_decomposition(rbind(residual, opposing, reversed))
stopifnot(isTRUE(all.equal(signed$d, c(0.06, 0.2, 0))),
          signed$i[2] < 0, signed$t[3] == -result$t,
          signed$i[3] == -result$i, signed$i_se[3] == result$i_se)
uncertain <- fixture
uncertain$a_se <- uncertain$a_se * 2
uncertain$comparator_se <- uncertain$comparator_se * 2
uncertain$t_se <- uncertain$t_se * 2
wide <- compute_decomposition(uncertain)
stopifnot(wide$i_se > result$i_se, wide$d_se > result$d_se,
          wide$i_hi - wide$i_lo > result$i_hi - result$i_lo)
cat('Residual, opposing signs, allele reversal and increasing uncertainty: OK\n')
unavailable <- fixture
unavailable$comparator_b <- NA_real_
unavailable$status <- 'unavailable'
unavailable$reason <- 'No other usable instruments'
missing <- compute_decomposition(rbind(fixture, unavailable))
stopifnot(nrow(missing) == 2L, missing$total_status[2] == 'complete',
          missing$indirect_status[2] == 'unavailable',
          is.na(missing$i[2]), is.na(missing$d[2]),
          missing$total_reason[2] %in% NA_character_,
          grepl('No other usable instruments', missing$indirect_reason[2]))
cat('Unavailable components preserve independently usable total: OK\n')
sensitivity <- decomposition_covariance_sensitivity(rbind(fixture, unavailable))
stopifnot(nrow(sensitivity) == 6L, identical(sort(unique(sensitivity$rho)), c(-0.3, 0, 0.3)),
          all(sensitivity$covariance_psd), all(sensitivity$t_se[1:3] == fixture$t_se))
zero <- subset(sensitivity, target_snp == 'rs1' & rho == 0 & decomposition_status == 'complete')
stopifnot(isTRUE(all.equal(zero$i_se, result$i_se)),
          isTRUE(all.equal(zero$d_se, result$d_se)),
          all(is.na(sensitivity$i_se[sensitivity$decomposition_status == 'unavailable'])))
# At rho=0.3, positive a,b increase indirect variance; shared t covariance
# reduces residual variance relative to simply adding Var(t)+Var(i).
positive <- subset(sensitivity, rho == 0.3 & decomposition_status == 'complete')
stopifnot(positive$i_se > result$i_se,
          positive$d_se^2 < fixture$t_se^2 + positive$i_se^2)
cat('Labelled PSD covariance sensitivity and zero-covariance primary consistency: OK\n')
missing_total <- fixture
missing_total$t <- NA_real_
partial <- compute_decomposition(missing_total)
stopifnot(partial$indirect_status == 'complete', is.finite(partial$i),
          partial$total_status == 'unavailable', is.na(partial$t), is.na(partial$d))
invalid_se <- fixture
invalid_se$a_se <- -0.02
invalid <- compute_decomposition(invalid_se)
stopifnot(invalid$indirect_status == 'unavailable', is.na(invalid$i),
          invalid$total_status == 'complete')
empty <- compute_decomposition(fixture[FALSE, ])
stopifnot(nrow(empty) == 0L, nrow(decomposition_covariance_sensitivity(fixture[FALSE, ])) == 0L)
# The one shared comparator is reused in each retained pair, without summing.
pairs <- rbind(fixture, fixture)
pairs$target_snp[2] <- 'rs2'
pairs$a[2] <- 0.3
pair_results <- compute_decomposition(pairs)
stopifnot(identical(pair_results$target_snp, pairs$target_snp),
          isTRUE(all.equal(pair_results$i, c(0.1, 0.15))))
cat('Partial availability, invalid SEs, empty inventory and shared mediator comparator: OK\n')
for (rho in unique(sensitivity$rho)) {
  correlation <- matrix(rho, 3, 3)
  diag(correlation) <- 1
  covariance <- diag(c(fixture$a_se, fixture$comparator_se, fixture$t_se)) %*%
    correlation %*% diag(c(fixture$a_se, fixture$comparator_se, fixture$t_se))
  stopifnot(min(eigen(covariance, symmetric=TRUE, only.values=TRUE)$values) >= -1e-12)
  # Independent scalar calculation includes all three covariance terms.
  gi <- c(fixture$comparator_b, fixture$a, 0)
  gd <- c(-fixture$comparator_b, -fixture$a, 1)
  scenario <- sensitivity[sensitivity$rho == rho & sensitivity$decomposition_status == 'complete', ]
  stopifnot(isTRUE(all.equal(scenario$i_se^2, as.numeric(t(gi) %*% covariance %*% gi))),
            isTRUE(all.equal(scenario$d_se^2, as.numeric(t(gd) %*% covariance %*% gd))))
}
cat('Independent covariance eigenvalue and gradient validation: OK\n')

available_input <- fixture
available_input$comparator_status <- 'available'
stopifnot(compute_decomposition(available_input)$decomposition_status == 'complete')
failed_comparator <- available_input
failed_comparator$comparator_status <- 'unavailable'
failed_comparator$comparator_reason <- 'Failed instrument check'
stopifnot(compute_decomposition(failed_comparator)$indirect_status == 'unavailable')
cat('Exact-hit comparator availability contract: OK\n')

# Outcome metadata must control units without assumptions about case/control coding.
units <- fixture
units$outcome_name <- 'CHD'
units$outcome_units <- 'Saved GWAS units; scale unavailable'
stopifnot(grepl('scale unavailable', compute_decomposition(units)$effect_scale),
          !grepl('GBC|log odds', compute_decomposition(units)$effect_scale),
          !grepl('GBC|log odds', compute_decomposition(fixture)$effect_scale))

  testthat::expect_true(TRUE)
})
