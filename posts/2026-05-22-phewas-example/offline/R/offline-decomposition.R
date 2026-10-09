# Adapted from MRCIEU/gbc-mr 3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
# Per-target decomposition in the recorded outcome GWAS units. a and b describe the
# same mediator scale; a and t must already share the recorded effect allele.
# First-order delta method, zero covariance working assumption in the primary
# analysis. Summary associations do not identify Cov(a,b), Cov(a,t), Cov(b,t).
# With parameter order (a,b,t), gradients are (b,a,0) and (-b,-a,1).
# Thus Var(i)=b^2 Var(a)+a^2 Var(b)+2ab Cov(a,b); Var(d)=Var(t)+
# Var(i)-2b Cov(a,t)-2a Cov(b,t). The residual is model-based, not an
# identified natural direct effect. Binary traits, non-collapsibility,
# nonlinearity, instrument validity and population differences limit causal
# interpretation. Separate mediators' indirect effects must not be summed.

.decomposition_interval <- function(estimate, se) {
  z <- qnorm(0.975)
  p <- ifelse(se > 0, 2 * pnorm(-abs(estimate / se)),
              ifelse(estimate == 0, 1, 0))
  list(se = se, lo = estimate - z * se, hi = estimate + z * se, pval = p)
}

compute_decomposition <- function(comparisons) {
  required <- c('mediator_id', 'target_snp', 'mediator', 'mediator_scale',
                'effect_allele', 'a', 'a_se', 't', 't_se',
                'comparator_b', 'comparator_se')
  missing <- setdiff(required, names(comparisons))
  if (length(missing)) stop('Missing decomposition columns: ', paste(missing, collapse = ', '))
  out <- comparisons
  usable_number <- function(x) is.numeric(x) & is.finite(x)
  usable_text <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
  total_ok <- usable_number(out$t) & usable_number(out$t_se) & out$t_se >= 0 &
    usable_text(out$effect_allele)
  indirect_ok <- usable_number(out$a) & usable_number(out$a_se) & out$a_se >= 0 &
    usable_number(out$comparator_b) & usable_number(out$comparator_se) & out$comparator_se >= 0 &
    usable_text(out$effect_allele) & usable_text(out$mediator_scale)
  if ('comparator_status' %in% names(out))
    indirect_ok <- indirect_ok & !is.na(out$comparator_status) & out$comparator_status %in% c('available', 'complete')
  out$i <- out$a * out$comparator_b
  out$d <- out$t - out$i
  out$i_se <- sqrt(out$comparator_b^2 * out$a_se^2 + out$a^2 * out$comparator_se^2)
  out$d_se <- sqrt(out$t_se^2 + out$i_se^2)
  for (component in c('t', 'i', 'd')) {
    se <- if (component == 't') out$t_se else out[[paste0(component, '_se')]]
    interval <- .decomposition_interval(out[[component]], se)
    for (field in names(interval)) out[[paste0(component, '_', field)]] <- interval[[field]]
  }
  out$total_status <- out$indirect_status <- out$residual_status <- out$decomposition_status <- rep('complete', nrow(out))
  out$total_reason <- out$indirect_reason <- out$residual_reason <- out$decomposition_reason <- rep(NA_character_, nrow(out))
  extra_reason <- rep('', nrow(out))
  for (column in intersect(c('reason', 'exact_reason', 'target_reason', 'comparator_reason'), names(out))) {
    use <- !is.na(out[[column]]) & nzchar(as.character(out[[column]]))
    extra_reason[use] <- paste(extra_reason[use], as.character(out[[column]][use]), sep = '; ')
  }
  component_ok <- list(total = total_ok, indirect = indirect_ok,
                       residual = total_ok & indirect_ok)
  base_reason <- c(total = 'Missing/unusable total association, standard error or effect allele',
                   indirect = 'Missing/unusable mediator association, shared mediator comparator, standard error, mediator scale or effect allele',
                   residual = 'Total and indirect components are both required')
  for (component in names(component_ok)) {
    ok <- component_ok[[component]]
    out[[paste0(component, '_status')]][!ok] <- 'unavailable'
    out[[paste0(component, '_reason')]][!ok] <- paste0(base_reason[[component]], extra_reason[!ok])
    prefix <- c(total = 't', indirect = 'i', residual = 'd')[[component]]
    for (column in c(prefix, paste0(prefix, c('_se', '_lo', '_hi', '_pval'))))
      out[[column]][!ok] <- NA_real_
  }
  incomplete <- !(total_ok & indirect_ok)
  out$decomposition_status[incomplete] <- 'unavailable'
  out$decomposition_reason[incomplete] <- vapply(which(incomplete), function(j) {
    paste(na.omit(c(out$total_reason[j], out$indirect_reason[j])), collapse = '; ')
  }, character(1))
  out$covariance_assumption <- rep('Zero Cov(a,b), Cov(a,t), Cov(b,t); unidentifiable from saved summary data', nrow(out))
  out$uncertainty_method <- rep('First-order delta method; approximate normal 95% intervals', nrow(out))
  units <- if ('outcome_units' %in% names(out)) as.character(out$outcome_units) else
    rep('Saved outcome GWAS units; scale unavailable', nrow(out))
  units[is.na(units) | !nzchar(trimws(units))] <- 'Saved outcome GWAS units; scale unavailable'
  out$effect_scale <- if (nrow(out)) paste0(units, ' per target effect allele') else character(0)
  out$residual_label <- rep('Model-based residual component relative to this mediator', nrow(out))
  out
}

# Prespecified hypothetical equicorrelation sensitivity, parameter order a,b,t.
# The 3x3 matrix has eigenvalues 1-rho (twice), 1+2*rho: the prespecified
# values -0.3, 0, +0.3 are positive definite. This is not an empirical estimate
# of sample overlap/covariance and does not identify any causal natural effect.
# Reversal of allele coding in an empirical covariance would also transform
# its entries; these fixed hypothetical correlations are separate scenarios.
decomposition_covariance_sensitivity <- function(comparisons) {
  primary <- compute_decomposition(comparisons)
  if (!nrow(primary)) {
    primary$rho <- numeric(0)
    primary$sensitivity <- character(0)
    primary$covariance_psd <- logical(0)
    return(primary)
  }
  rows <- rep(seq_len(nrow(primary)), each = 3L)
  out <- primary[rows, , drop = FALSE]
  rownames(out) <- NULL
  out$rho <- rep(c(-0.3, 0, 0.3), times = nrow(primary))
  out$sensitivity <- paste0('Hypothetical equicorrelation rho=', out$rho)
  out$covariance_psd <- 1 - out$rho >= 0 & 1 + 2 * out$rho >= 0
  out$covariance_assumption <- paste0('Hypothetical Corr(a,b)=Corr(a,t)=Corr(b,t)=',
                                     out$rho, '; empirical covariances unidentifiable')
  indirect_ok <- out$indirect_status == 'complete'
  residual_ok <- out$residual_status == 'complete'
  vi <- out$comparator_b^2 * out$a_se^2 + out$a^2 * out$comparator_se^2 +
    2 * out$a * out$comparator_b * out$rho * out$a_se * out$comparator_se
  vd <- out$t_se^2 + vi - 2 * out$comparator_b * out$rho * out$a_se * out$t_se -
    2 * out$a * out$rho * out$comparator_se * out$t_se
  # Valid covariance matrices imply nonnegative variances; guard tiny rounding.
  out$i_se[indirect_ok] <- sqrt(pmax(0, vi[indirect_ok]))
  out$d_se[residual_ok] <- sqrt(pmax(0, vd[residual_ok]))
  for (component in c('i', 'd')) {
    interval <- .decomposition_interval(out[[component]], out[[paste0(component, '_se')]])
    for (field in names(interval)) out[[paste0(component, '_', field)]] <- interval[[field]]
  }
  out
}
