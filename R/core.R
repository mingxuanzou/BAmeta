# Bias-aware meta-analysis for independent study-level estimates.
# The implementation follows Zou and Ding's treatment-control BAMA model.

.bama_numeric <- function(x, name, length_required = NULL) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      !length(x) || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a nonempty vector of finite numbers.", call. = FALSE)
  }
  if (!is.null(length_required) && length(x) != length_required) {
    stop(name, " must have length ", length_required, ".", call. = FALSE)
  }
  unname(as.double(x))
}

.bama_check_fit <- function(object) {
  if (!inherits(object, "bama_fit")) {
    stop("object must be a fit returned by ba_meta().", call. = FALSE)
  }
  invisible(object)
}

.bama_bias <- function(b) {
  b <- .bama_numeric(b, "b")
  if (any(b < 0)) stop("b must be nonnegative.", call. = FALSE)
  b
}

# Subset sums for one half of the meet-in-the-middle partition calculation.
.bama_subset_sums <- function(x) {
  sums <- 0
  for (value in x) sums <- c(sums, sums + value)
  sums
}

.bama_partition <- function(weights) {
  # If one weight is at least one half, placing it against all other weights
  # is optimal; this also avoids cancellation for strongly unequal weights.
  largest <- which.max(weights)
  if (weights[largest] >= 0.5) {
    small <- sum(weights[-largest])
    return(list(kappa = 4 * weights[largest] * small,
                partition_residual = weights[largest] - small))
  }
  split <- floor(length(weights) / 2)
  left <- .bama_subset_sums(weights[seq_len(split)])
  right <- sort(.bama_subset_sums(weights[-seq_len(split)]))
  target <- 0.5 - left
  index <- findInterval(target, right)
  below <- pmax(1L, index)
  above <- pmin(length(right), index + 1L)
  candidates <- c(left + right[below], left + right[above])
  best <- candidates[which.min(abs(candidates - 0.5))]
  # 4 p (1-p) is algebraically 1 - (1-2p)^2 but more stable near p=0.
  list(kappa = min(1, max(0, 4 * best * (1 - best))),
       partition_residual = abs(1 - 2 * best))
}

.bama_wald <- function(estimate, se, null, level) {
  critical <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  statistic <- (estimate - null) / se
  list(estimate = estimate, se = se, t = statistic,
       p = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
       lower = estimate - critical * se,
       upper = estimate + critical * se)
}

# Solve for L = mu + offset rather than L itself. The offset always lies
# between z_(1-alpha) and z_(1-alpha/2), including when mu is large.
.bama_folded_quantile <- function(mu, alpha) {
  upper <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  if (mu == 0) return(upper)
  lower <- stats::qnorm(alpha, lower.tail = FALSE)
  tail_minus_alpha <- function(offset) {
    stats::pnorm(offset, lower.tail = FALSE) +
      stats::pnorm(-2 * mu - offset) - alpha
  }
  if (tail_minus_alpha(lower) <= 0) return(mu + lower)
  if (tail_minus_alpha(upper) >= 0) return(mu + upper)
  mu + stats::uniroot(tail_minus_alpha, c(lower, upper),
                     tol = 1e-12)$root
}

.bama_significance_threshold <- function(t_abs, J, alpha) {
  critical <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  if (t_abs <= critical) return(0)
  lower <- max(-t_abs, stats::qnorm(alpha / 2))
  upper <- stats::qnorm(alpha)
  tail_minus_alpha <- function(offset) {
    stats::pnorm(offset) + stats::pnorm(-2 * t_abs - offset) - alpha
  }
  if (tail_minus_alpha(upper) <= 0) offset <- upper else {
    offset <- stats::uniroot(tail_minus_alpha, c(lower, upper),
                            tol = 1e-12)$root
  }
  max(0, (t_abs + offset) / sqrt(J))
}

.bama_q_tail <- function(Q, df, ncp) {
  if (ncp == 0) return(stats::pchisq(Q, df = df, lower.tail = FALSE))
  if (is.infinite(ncp)) return(1)
  stats::pchisq(Q, df = df, ncp = ncp, lower.tail = FALSE)
}

.bama_compatibility_threshold <- function(Q, J, kappa, alpha) {
  if (.bama_q_tail(Q, J - 1L, 0) >= alpha) return(0)
  # Invert in sqrt(ncp); the final scaling is applied only after inversion.
  f <- function(x) .bama_q_tail(Q, J - 1L, x * x) - alpha
  upper <- sqrt(max(1, Q))
  while (f(upper) < 0) {
    upper <- 2 * upper
    if (!is.finite(upper)) {
      stop("Could not bracket the Q-test threshold at this numeric scale.",
           call. = FALSE)
    }
  }
  stats::uniroot(f, c(0, upper), tol = 1e-11)$root / sqrt(J * kappa)
}

ba_meta <- function(estimate, se = NULL, variance = NULL, study = NULL,
                    level = 0.95, null = 0,
                    q_method = c("auto", "exact", "conservative"),
                    exact_max = 36L) {
  estimate <- .bama_numeric(estimate, "estimate")
  J <- length(estimate)
  if (J < 2L) stop("At least two independent studies are required.", call. = FALSE)
  if (is.null(se) == is.null(variance)) {
    stop("Supply exactly one of se and variance.", call. = FALSE)
  }
  if (!is.null(se)) {
    se <- .bama_numeric(se, "se", J)
    if (any(se <= 0)) stop("All standard errors must be positive.", call. = FALSE)
    variance <- se^2
  } else variance <- .bama_numeric(variance, "variance", J)
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("All variances must be positive and finite at double precision.",
         call. = FALSE)
  }
  level <- .bama_numeric(level, "level", 1L)
  if (level <= 0.5 || level >= 1) {
    stop("level must be strictly between 0.5 and 1.", call. = FALSE)
  }
  null <- .bama_numeric(null, "null", 1L)
  exact_max <- .bama_numeric(exact_max, "exact_max", 1L)
  if (exact_max < 2 || exact_max > 36 || exact_max != floor(exact_max)) {
    stop("exact_max must be an integer between 2 and 36.", call. = FALSE)
  }
  q_method <- match.arg(q_method)
  requested_method <- q_method
  if (q_method == "auto") {
    q_method <- if (J <= exact_max) "exact" else "conservative"
    if (q_method == "conservative") {
      message("Using the conservative Q-test bound J*b^2 because J exceeds exact_max.")
    }
  }
  if (q_method == "exact" && J > exact_max) {
    stop("Exact partition calculation requires J <= exact_max (at most 36). ",
         "Use q_method = 'conservative' for larger analyses.", call. = FALSE)
  }
  if (is.null(study)) study <- paste("Study", seq_len(J))
  if (length(study) != J || anyNA(study)) {
    stop("study must contain one nonmissing identifier per estimate.", call. = FALSE)
  }
  study <- as.character(study)
  if (any(!nzchar(trimws(study))) || anyDuplicated(study)) {
    stop("Study identifiers must be nonempty and unique.", call. = FALSE)
  }

  # Relative precisions avoid overflow from taking reciprocals of small
  # variances. Reject values beyond floating-point resolution explicitly.
  minimum_variance <- min(variance)
  relative <- minimum_variance / variance
  if (any(relative == 0)) {
    stop("The range of study variances exceeds double-precision resolution.",
         call. = FALSE)
  }
  total <- sum(relative)
  weights <- relative / total
  fe_se <- sqrt(minimum_variance / total)
  fe_estimate <- sum(weights * estimate)
  residual <- estimate - fe_estimate
  Q <- sum((residual / sqrt(variance))^2)
  # sum_{j<k} relative[j] relative[k] is stable when one study dominates.
  cross <- sum(relative[-1L] * cumsum(relative)[-J])
  denominator <- 2 * cross / total
  tau2 <- max(0, Q - J + 1) * minimum_variance / denominator
  if (!is.finite(fe_estimate) || !is.finite(fe_se) || fe_se <= 0 ||
      !is.finite(Q) || !is.finite(tau2)) {
    stop("Input magnitudes exceed the supported numerical range; rescale the effects and standard errors.",
         call. = FALSE)
  }
  re_variance <- variance + tau2
  if (any(!is.finite(re_variance))) {
    stop("Random-effects variances overflow; rescale the effects and standard errors.",
         call. = FALSE)
  }
  re_relative <- min(re_variance) / re_variance
  re_weights <- re_relative / sum(re_relative)
  re_estimate <- sum(re_weights * estimate)
  re_se <- sqrt(min(re_variance) / sum(re_relative))
  dagger_se <- sqrt(fe_se^2 + tau2 * sum(weights^2))

  partition <- if (q_method == "exact") .bama_partition(weights) else {
    list(kappa = 1, partition_residual = NA_real_)
  }
  if (partition$kappa <= 0 || !is.finite(partition$kappa)) {
    stop("The precision weights are too unequal for a stable Q-test calculation.",
         call. = FALSE)
  }
  fe <- .bama_wald(fe_estimate, fe_se, null, level)
  re <- .bama_wald(re_estimate, re_se, null, level)
  re$tau2 <- tau2
  if (!is.finite(fe$t) || !is.finite(re$t) || !is.finite(re_se) ||
      !is.finite(dagger_se) || re_se <= 0 || dagger_se <= 0) {
    stop("The standardized estimates exceed the supported numerical range.",
         call. = FALSE)
  }
  result <- list(
    call = match.call(), J = J, estimate = estimate, variance = variance,
    study = study, weights = weights, level = level, null = null,
    scale = sqrt(J) * fe_se, fe = fe, re = re,
    fe_dagger = .bama_wald(fe_estimate, dagger_se, null, level),
    Q = Q, Q_df = J - 1L, Q_p = .bama_q_tail(Q, J - 1L, 0),
    q = c(list(method = q_method, requested_method = requested_method,
               exact_max = exact_max), partition)
  )
  class(result) <- "bama_fit"
  alpha <- 1 - level
  b_significance <- .bama_significance_threshold(abs(fe$t), J, alpha)
  b_compatibility <- .bama_compatibility_threshold(Q, J, partition$kappa, alpha)
  optimistic <- bama_ci(result, b_compatibility)
  result$thresholds <- list(
    b_significance = b_significance, b_compatibility = b_compatibility,
    approx_significance = abs(fe$t) / sqrt(J),
    approx_compatibility = sqrt(tau2) / result$scale,
    optimistic_lower = optimistic$lower, optimistic_upper = optimistic$upper,
    optimistic_valid = FALSE,
    compatible_significant = b_compatibility < b_significance,
    q_method = q_method
  )
  result
}

bama_ci <- function(object, b) {
  .bama_check_fit(object)
  b <- .bama_bias(b)
  mu <- sqrt(object$J) * b
  if (any(!is.finite(mu))) {
    stop("b is too large to represent sqrt(J)*b at double precision.", call. = FALSE)
  }
  critical <- vapply(mu, .bama_folded_quantile, numeric(1),
                     alpha = 1 - object$level)
  radius <- object$fe$se * critical
  data.frame(b = b, estimate = object$fe$estimate, se = object$fe$se,
             critical = critical, lower = object$fe$estimate - radius,
             upper = object$fe$estimate + radius, level = object$level)
}

bama_pvalue <- function(object, b) {
  .bama_check_fit(object)
  b <- .bama_bias(b)
  mu <- sqrt(object$J) * b
  t_abs <- abs(object$fe$t)
  pmin(1, stats::pnorm(mu - t_abs) + stats::pnorm(-mu - t_abs))
}

bama_q_pvalue <- function(object, b) {
  .bama_check_fit(object)
  b <- .bama_bias(b)
  ncp <- (sqrt(object$J * object$q$kappa) * b)^2
  vapply(ncp, function(x) .bama_q_tail(object$Q, object$Q_df, x), numeric(1))
}

bama_thresholds <- function(object) {
  .bama_check_fit(object)
  object$thresholds
}

bama_diagnostics <- function(object) {
  .bama_check_fit(object)
  # Recompute each deletion using relative precisions. This avoids subtracting
  # nearly equal numerator/denominator values for a dominant study.
  deleted <- vapply(seq_len(object$J), function(j) {
    keep <- seq_len(object$J) != j
    relative <- min(object$variance[keep]) / object$variance[keep]
    sum(relative / sum(relative) * object$estimate[keep])
  }, numeric(1))
  shift <- deleted - object$fe$estimate
  data.frame(
    study = object$study, estimate = object$estimate,
    se = sqrt(object$variance), weight = object$weights,
    b_j = abs(object$estimate - object$fe$estimate) / object$scale,
    leave_one_out_estimate = deleted, leave_one_out_shift = shift,
    standardized_shift = abs(shift) / object$scale,
    stringsAsFactors = FALSE
  )
}
