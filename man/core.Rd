\name{ba_meta}
\alias{ba_meta}
\alias{bama_ci}
\alias{bama_pvalue}
\alias{bama_q_pvalue}
\alias{bama_thresholds}
\alias{bama_diagnostics}
\title{Bias-Aware Meta-Analysis of Independent Study Estimates}
\description{
Fits the treatment-control bias-aware meta-analysis (BAMA) model of Zou and
Ding, with fixed, unknown study biases bounded on a common scale. Computes
FEMA, DerSimonian--Laird REMA, FEMA with REMA-adjusted variance,
sensitivity intervals, significance and compatibility thresholds, and
study-level dispersion and leave-one-study-out diagnostics.
}
\usage{
ba_meta(estimate, se = NULL, variance = NULL, study = NULL,
  level = 0.95, null = 0,
  q_method = c("auto", "exact", "conservative"), exact_max = 36L)

bama_ci(object, b)

bama_pvalue(object, b)

bama_q_pvalue(object, b)

bama_thresholds(object)

bama_diagnostics(object)
}
\arguments{
  \item{estimate}{Finite numeric vector of study-level effect estimates,
  with at least two entries. One entry must represent one independent study.}
  \item{se}{Positive finite standard errors. Supply exactly one of
  \code{se} and \code{variance}.}
  \item{variance}{Positive finite sampling variances.}
  \item{study}{Optional unique, nonmissing study identifiers.}
  \item{level}{Confidence level, strictly between 0.5 and 1. The same
  significance level \eqn{1-\mathrm{level}} is used for model checking.}
  \item{null}{Finite scalar null effect. The paper uses zero; all
  significance calculations here use the supplied null value.}
  \item{q_method}{Computation of the maximum noncentrality for the Q-test.
  \code{"exact"} solves the weighted number-partitioning problem;
  \code{"conservative"} uses the upper bound \eqn{Jb^2};
  \code{"auto"} uses the exact method when \eqn{J} is at most
  \code{exact_max}, otherwise announces the conservative bound.}
  \item{exact_max}{Largest number of studies for exact partitioning;
  an integer from 2 to 36. The upper limit prevents exponential memory growth.}
  \item{object}{A \code{bama_fit} object returned by \code{ba_meta}.}
  \item{b}{Nonempty vector of finite, nonnegative standardized bias bounds.}
}
\details{
Let \eqn{J} be the number of studies, \eqn{\sigma_j^2} their sampling
variances, \eqn{\lambda_j=\sigma_j^{-2}/\sum_k\sigma_k^{-2}},
\eqn{\widehat\tau_{FE}=\sum_j\lambda_j\widehat\tau_j}, and
\eqn{\sigma_{FE}^2=1/\sum_j\sigma_j^{-2}}.
BAMA bounds each fixed study bias by
\eqn{|\Delta_j|\leq b\sqrt J\sigma_{FE}}.
For a prespecified \eqn{b}, the confidence interval is centered at
\eqn{\widehat\tau_{FE}} with radius
\eqn{\sigma_{FE}L_\alpha(\sqrt Jb)}, where \eqn{L_\alpha(x)} is the
\eqn{1-\alpha} quantile of \eqn{|N(x,1)|} and
\eqn{\alpha=1-\mathrm{level}}. At \eqn{b=0} this is the FEMA interval.

The maximum bias for significance is the smallest \eqn{b} at which this
interval includes \code{null}; it is zero when FEMA is already nonsignificant.
\code{bama_pvalue} returns the corresponding worst-case two-sided target
test p-value. \code{bama_q_pvalue} is a different test: it checks compatibility
of the specified BAMA model using Cochran's Q, \eqn{J-1} degrees of freedom,
and maximum noncentrality
\eqn{Jb^2\{1-\min_{s_j\in\{-1,1\}}(\sum_j\lambda_js_j)^2\}}.
The minimum bias for compatibility is the smallest \eqn{b} for which this
p-value reaches \eqn{\alpha}. It is zero when the ordinary Q-test does
not reject FEMA.

Exact partitioning uses meet-in-the-middle subset sums, requiring exponential
memory and time in half the study count. The conservative method replaces
the bracketed factor by 1. Its Q p-values are at least the exact values, and
its compatibility threshold can be smaller. The selected method is recorded
in the fit and threshold results; conservative results must not be labeled
as exact. The significance threshold and BAMA confidence intervals are
unaffected by this computational choice.

The REMA variance component is the nonnegative DerSimonian--Laird estimate.
The FEMA-dagger variance is
\eqn{\sigma_{FE}^2+\widehat\sigma_B^2\sum_j\lambda_j^2}.
Normal critical values are used for all three conventional intervals;
Hartung--Knapp adjustments are not applied.

The diagnostic \eqn{b_j=|\widehat\tau_j-\widehat\tau_{FE}|/
(\sqrt J\sigma_{FE})} measures observed dispersion, not an identified
absolute bias. It differs from the standardized absolute leave-one-out shift
by the multiplier \eqn{\lambda_j/(1-\lambda_j)}. Both use the full-data
reference scale. Neither can detect a bias common to every study or establish
an upper confidence bound on the BAMA bias level.

The optimistic interval evaluates the sensitivity interval at the
data-dependent compatibility threshold. It is a descriptive interval and
does not have a nominal coverage guarantee. The compatibility threshold
does not provide the upper bias bound needed for the BAMA guarantee.

The finite-sample BAMA guarantees assume independent Gaussian estimates
with known variances and the specified fixed-bias bound. Using estimated
standard errors relies on the corresponding large-sample approximation.
Missing, nonfinite, or duplicate study inputs are rejected rather than
silently omitted. Multiple dependent estimates from the same trial must
not be entered as independent studies.
}
\value{
\code{ba_meta} returns an object of class \code{bama_fit} containing the
input vectors, normalized precision \code{weights}, \code{J},
\code{level}, \code{null}, and reference \code{scale}. Components
\code{fe}, \code{re}, and \code{fe_dagger} are lists containing
\code{estimate}, \code{se}, \code{t}, \code{p}, \code{lower}, and
\code{upper}; \code{re} also contains \code{tau2}. The fit includes
\code{Q}, \code{Q_df}, \code{Q_p}, and a \code{q} list with the
selected \code{method}, \code{kappa} (so noncentrality is
\eqn{J\,\mathrm{kappa}\,b^2}), and exact \code{partition_residual}
(\code{NA} for the conservative method).

\code{bama_ci} returns a data frame with \code{b}, \code{estimate},
\code{se}, \code{critical}, \code{lower}, \code{upper}, and \code{level}.
\code{bama_pvalue} and \code{bama_q_pvalue} return numeric vectors.

\code{bama_thresholds} returns a list with \code{b_significance},
\code{b_compatibility}, their first-order \code{approx_significance}
and \code{approx_compatibility}, \code{optimistic_lower},
\code{optimistic_upper}, \code{optimistic_valid} (always \code{FALSE}),
\code{compatible_significant}, and \code{q_method}. The same list is
stored in the fit's \code{thresholds} component.

\code{bama_diagnostics} returns a data frame containing \code{study},
\code{estimate}, \code{se}, \code{weight}, \code{b_j},
\code{leave_one_out_estimate}, the signed \code{leave_one_out_shift}
(deleted estimate minus full estimate), and \code{standardized_shift}
(absolute shift divided by the full-data reference scale).
}
\references{
Zou, M. and Ding, P. (2026). Bias-aware meta-analysis: beyond fixed-effect
and random-effects meta-analyses. Manuscript.

DerSimonian, R. and Laird, N. (1986). Meta-analysis in clinical trials.
\emph{Controlled Clinical Trials}, \bold{7}, 177--188.
}
\examples{
fit <- ba_meta(c(-1.2, -0.8, -1.5, -0.9),
               se = c(0.2, 0.3, 0.25, 0.2))
fit$fe
bama_ci(fit, b = c(0, 0.5, 1))
bama_thresholds(fit)
bama_diagnostics(fit)
bama_q_pvalue(fit, b = c(0, 0.5, 1))
}
