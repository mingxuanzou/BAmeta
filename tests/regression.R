# Independent numerical regression and mathematical invariants.
# Uses only base/recommended R packages, so R CMD check needs no test runner.
library(BaMeta)

close_to <- function(actual, expected, tolerance = 1e-8, label = "numeric result") {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(!is.finite(expected)) ||
      any(abs(actual - expected) > tolerance * pmax(1, abs(expected)))) {
    stop(label, " differs from its independently derived reference.",
         call. = FALSE)
  }
  invisible(TRUE)
}

must_error <- function(expr) {
  if (!inherits(tryCatch(force(expr), error = identity), "error")) {
    stop("Invalid input did not produce an error.", call. = FALSE)
  }
  invisible(TRUE)
}

# These constants were calculated independently with SciPy, using the saved
# reconstructed data and exhaustive enumeration of all subsets for Q's bound.
# They also reproduce every rounded entry in the manuscript's Panels A and B.
reference <- list(
  TAD = c(J = 15, fe = -13.078240344768409, fe_se = 0.16157232874894792,
          Q = 228.40884666143654, Q_p = 8.196601378459031e-41,
          tau2 = 6.302037053934121, re = -12.602806188759894,
          re_se = 0.679838689768677, dagger_se = 0.8549078017605278,
          residue = 1.7051536205570095e-05,
          b_min = 3.3578418436061153, b_max = 20.474839315565237,
          min_approx = 4.011694583720947, max_approx = 20.899538695840707,
          optimistic_lower = -15.445229586007635,
          optimistic_upper = -10.711251103529182,
          max_bj = 7.476010727470836),
  PAT = c(J = 8, fe = -11.266538291013026, fe_se = 0.6370243454615554,
          Q = 15.215489053406849, Q_p = 0.03333424725340622,
          tau2 = 4.00806697418175, re = -11.529917321433087,
          re_se = 0.9883617047972263, dagger_se = 1.039240778517284,
          residue = 0.0015152572586160096,
          b_min = 0.27026936016518155, b_max = 5.6714714631920256,
          min_approx = 1.1111340620055707, max_approx = 6.253015040030363,
          optimistic_lower = -12.805876341568432,
          optimistic_upper = -9.72720024045762,
          max_bj = 2.3124606623268553)
)

dat <- williams_data(format = "data.frame")
stopifnot(nrow(dat) == 23L, !any(grepl("Arns", dat$study)))
for (group in names(reference)) {
  d <- dat[dat$class == group, , drop = FALSE]
  fit <- bama(d$estimate, variance = d$variance, study = d$study,
                 q_method = "exact")
  th <- bama_thresholds(fit)
  dg <- bama_diagnostics(fit)
  observed <- c(J = fit$J, fe = fit$fe$estimate, fe_se = fit$fe$se,
                Q = fit$Q, Q_p = fit$Q_p, tau2 = fit$re$tau2,
                re = fit$re$estimate, re_se = fit$re$se,
                dagger_se = fit$fe_dagger$se,
                residue = fit$q$partition_residual,
                b_min = th$b_compatibility, b_max = th$b_significance,
                min_approx = th$approx_compatibility,
                max_approx = th$approx_significance,
                optimistic_lower = th$optimistic_lower,
                optimistic_upper = th$optimistic_upper,
                max_bj = max(dg$b_j))
  stopifnot(identical(names(observed), names(reference[[group]])))
  close_to(observed, reference[[group]], label = paste(group, "results"))
  # Relative accuracy is meaningful for the very small TAD central Q tail.
  close_to(fit$Q_p / reference[[group]]["Q_p"], 1, 1e-10, "Q tail")

  ci <- bama_ci(fit, c(0, th$b_compatibility, th$b_significance))
  close_to(c(ci$lower[1], ci$upper[1]), c(fit$fe$lower, fit$fe$upper))
  close_to(ci$upper[3], fit$null, label = "CI touching the null")
  close_to(bama_pvalue(fit, th$b_significance), 1 - fit$level)
  close_to(bama_q_pvalue(fit, th$b_compatibility), 1 - fit$level)
  close_to(bama_q_pvalue(fit, 0), fit$Q_p)
  stopifnot(identical(th$optimistic_valid, FALSE),
            identical(th$compatible_significant,
                      th$b_compatibility < th$b_significance))

  grid <- seq(0, th$b_significance * 1.1, length.out = 50)
  curve <- bama_ci(fit, grid)
  stopifnot(all(diff(curve$upper) >= 0), all(diff(curve$lower) <= 0),
            all(diff(bama_q_pvalue(fit, grid)) >= -1e-12))

  # Check deleted estimates by explicit refitting, separately from the
  # algebraic leverage identity used to interpret standardized shifts.
  deleted <- vapply(seq_len(nrow(d)), function(j) {
    bama(d$estimate[-j], variance = d$variance[-j])$fe$estimate
  }, numeric(1))
  close_to(dg$leave_one_out_estimate, deleted)
  close_to(dg$standardized_shift,
           fit$weights / (1 - fit$weights) * dg$b_j)
}

# Check the exact weighted partition against an independent, intentionally
# small brute-force sign enumeration. These cases include equal odd/even J,
# dominant weights, an exact balance, and unequal weights without a balance.
weight_cases <- list(rep(1 / 4, 4), rep(1 / 5, 5), c(0.8, 0.2),
                     c(0.6, 0.1, 0.1, 0.2), c(0.1, 0.2, 0.3, 0.4),
                     c(1, 2, 4, 8, 16) / 31,
                     c(0.13, 0.17, 0.21, 0.23, 0.26))
for (weights in weight_cases) {
  J <- length(weights)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), J)))
  residue <- min(abs(as.vector(signs %*% weights)))
  fit <- bama(seq_len(J), variance = 1 / weights, q_method = "exact")
  close_to(fit$q$partition_residual, residue, 1e-12, "partition residue")
  close_to(fit$q$kappa, 1 - residue^2, 1e-12, "partition coefficient")
}

# Direct upper-tail calculations agree with the Q function. A conservative
# noncentrality bound must give larger p-values and a smaller lower threshold.
y <- c(-3, -1, 0.4, 2.1, 4.8)
se <- c(0.3, 0.5, 0.6, 0.8, 1.1)
fit <- bama(y, se = se)
conservative <- bama(y, se = se, q_method = "conservative")
grid <- c(0, 0.1, 0.5, 1, 2, 4)
close_to(bama_q_pvalue(fit, grid),
         stats::pchisq(fit$Q, df = length(y) - 1,
                       ncp = length(y) * fit$q$kappa * grid^2,
                       lower.tail = FALSE))
stopifnot(all(bama_q_pvalue(conservative, grid) >=
              bama_q_pvalue(fit, grid) - 1e-12),
          conservative$thresholds$b_compatibility <=
            fit$thresholds$b_compatibility + 1e-12,
          identical(conservative$q$method, "conservative"))

invariant_values <- function(fit) {
  c(fit$Q, fit$Q_p, fit$q$kappa,
    unlist(fit$thresholds[c("b_significance", "b_compatibility",
                           "approx_significance", "approx_compatibility")]))
}

# Change the unit, move the origin together with the tested null, and reorder
# studies. These transformations should preserve all dimensionless results.
scaled <- bama(7.3 * y, se = 7.3 * se)
shifted <- bama(y + 12.5, se = se, null = 12.5)
order <- c(4, 1, 5, 2, 3)
permuted <- bama(y[order], se = se[order])
for (other in list(scaled, shifted, permuted)) {
  close_to(invariant_values(other), invariant_values(fit))
}
close_to(bama_diagnostics(scaled)$b_j, bama_diagnostics(fit)$b_j)
close_to(bama_diagnostics(shifted)$b_j, bama_diagnostics(fit)$b_j)
close_to(bama_diagnostics(permuted)$b_j,
         bama_diagnostics(fit)$b_j[order])
close_to(scaled$fe$estimate, 7.3 * fit$fe$estimate)
close_to(scaled$re$tau2, 7.3^2 * fit$re$tau2)

# A sign change reflects the intervals and preserves absolute thresholds.
reversed <- bama(-y, se = se)
close_to(invariant_values(reversed), invariant_values(fit))
original_ci <- bama_ci(fit, grid)
reversed_ci <- bama_ci(reversed, grid)
close_to(reversed_ci$lower, -original_ci$upper)
close_to(reversed_ci$upper, -original_ci$lower)

# The b=0 nested model and non-significant limit must not need root-finding.
zero <- bama(rep(0, 4), se = rep(1, 4))
stopifnot(zero$Q == 0, zero$Q_p == 1, zero$re$tau2 == 0,
          zero$thresholds$b_significance == 0,
          zero$thresholds$b_compatibility == 0,
          !zero$thresholds$compatible_significant,
          all(bama_diagnostics(zero)$b_j == 0),
          all(bama_q_pvalue(zero, c(0, 1, 10)) == 1))
close_to(bama_ci(zero, 0)$critical, stats::qnorm(0.975))
nonsignificant <- bama(c(-0.2, 0, 0.2), se = rep(1, 3))
stopifnot(nonsignificant$thresholds$b_significance == 0)

# Strong signals and confidence levels close to either allowed endpoint
# exercise the folded-normal inversion independently through its tail equation.
for (level in c(0.51, 0.95, 0.999999)) {
  strong <- bama(c(80, 81, 79), se = c(1, 1.5, 2), level = level)
  boundary <- strong$thresholds$b_significance
  close_to(bama_pvalue(strong, boundary), 1 - level, 1e-10)
  close_to(bama_ci(strong, boundary)$lower, 0, 1e-8)
  mu <- sqrt(strong$J) * c(0, 0.01, 1, 100)
  critical <- bama_ci(strong, mu / sqrt(strong$J))$critical
  tail <- stats::pnorm(mu - critical) + stats::pnorm(-mu - critical)
  close_to(tail / (1 - level), rep(1, length(mu)), 1e-7)
}

# Validation must report unsupported/ambiguous inputs rather than silently
# coercing values, dropping studies, or selecting an uncertainty measure.
must_error(bama(1, se = 1))
must_error(bama(numeric(), se = numeric()))
must_error(bama(c(1, 2), se = 1))
must_error(bama(c(1, NA_real_), se = c(1, 1)))
must_error(bama(c(1, Inf), se = c(1, 1)))
must_error(bama(c("1", "2"), se = c(1, 1)))
must_error(bama(c(1, 2), se = c(0, 1)))
must_error(bama(c(1, 2), variance = c(-1, 1)))
must_error(bama(c(1, 2), se = c(1, 1), variance = c(1, 1)))
must_error(bama(c(1, 2)))
must_error(bama(c(1, 2), se = c(1, 1), level = 0.5))
must_error(bama(c(1, 2), se = c(1, 1), level = 1))
must_error(bama(c(1, 2), se = c(1, 1), null = NA_real_))
must_error(bama(c(1, 2), se = c(1, 1), study = c("x", "x")))
must_error(bama(c(1, 2), se = c(1, 1), study = c("x", "")))
must_error(bama(c(1, 2, 3), se = c(1, 1, 1), q_method = "exact",
                   exact_max = 2))
must_error(bama_ci(fit, -1))
must_error(bama_ci(fit, NA_real_))
must_error(bama_q_pvalue(fit, Inf))
must_error(bama_pvalue(fit, numeric()))

# Data preparation and fitting are separate. The documented vector, list,
# data-frame and grouped-input forms must produce the same analyses.
inputs <- williams_data()
fits <- bama(inputs)
stopifnot(identical(names(inputs), c("Open-label TADs", "PAT")),
          identical(names(fits), names(inputs)),
          all(vapply(inputs, function(x) identical(names(x),
                     c("estimate", "se", "study")), logical(1))))
legacy <- williams_example()
for (i in seq_along(inputs)) {
  d <- inputs[[i]]
  direct <- bama(estimate = d$estimate, se = d$se, study = d$study)
  single <- bama(williams_data(names(inputs)[i]))
  framed <- bama(as.data.frame(d))
  variance_input <- bama(list(estimate = d$estimate, variance = d$se^2,
                             study = d$study))
  for (other in list(direct, single, framed, variance_input, legacy[[i]])) {
    close_to(invariant_values(other), invariant_values(fits[[i]]))
    close_to(other$fe$estimate, fits[[i]]$fe$estimate)
    stopifnot(identical(other$study, fits[[i]]$study))
  }
  stopifnot(identical(fits[[i]]$label, names(inputs)[i]),
            identical(single$label, names(inputs)[i]))
  close_to(c(fits[[i]]$thresholds$b_significance,
             fits[[i]]$thresholds$b_compatibility),
           reference[[i]][c("b_max", "b_min")])
}
stopifnot(identical(williams_data("TAD"),
                    williams_data("Open-label TADs")),
          nrow(williams_data("PAT", format = "data.frame")) == 8L)
aliased <- ba_meta(inputs[["PAT"]])
close_to(invariant_values(aliased), invariant_values(fits[["PAT"]]))
configured <- bama(inputs, level = 0.90, null = -1,
                   q_method = "conservative", exact_max = 2L)
for (i in seq_along(inputs)) {
  expected <- bama(inputs[[i]], level = 0.90, null = -1,
                   q_method = "conservative", exact_max = 2L)
  close_to(invariant_values(configured[[i]]), invariant_values(expected))
  stopifnot(configured[[i]]$level == 0.90, configured[[i]]$null == -1,
            configured[[i]]$q$exact_max == 2,
            configured[[i]]$q$method == "conservative")
}
must_error(bama(inputs, q_method = "exact", exact_max = 2L))
must_error(bama(inputs[[1]], se = inputs[[1]]$se))
must_error(bama(inputs, study = "override"))
must_error(bama(list()))
must_error(bama(unname(inputs)))
must_error(bama(list(estimate = c(1, 2), se = c(1, 1), se = c(2, 2))))
must_error(bama(list(estimate = c(1, 2), se = c(1, 1), extra = 1)))
must_error(bama(list(estimate = c(1, 2), se = c(1, 1), variance = c(1, 1))))
must_error(bama(list(estimates = c(1, 2), se = c(1, 1))))
must_error(bama(list(estimate = list(1, 2), se = c(1, 1))))
must_error(bama(list(PAT = inputs[[2]], PAT = inputs[[2]])))
must_error(bama(setNames(inputs, c("PAT", ""))))
must_error(bama(setNames(inputs, c("PAT", NA_character_))))
must_error(bama(list(PAT = inputs[[2]], other = list(se = c(1, 1)))))
must_error(williams_data("unknown"))
must_error(williams_data(c("TAD", "PAT")))
must_error(williams_data(NA_character_))

cat("All BaMeta numerical regression and invariant checks passed.\n")
