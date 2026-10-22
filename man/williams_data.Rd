\name{williams_data}
\alias{williams_data}
\alias{williams_example}
\title{Reconstructed Williams Example}
\description{Load the saved study summaries used in the BAMA manuscript or fit the two groups separately.}
\usage{
williams_data()
williams_example(level = 0.95, null = 0, q_method = "exact")
}
\arguments{
\item{level}{Confidence level, strictly between 0.5 and 1.}
\item{null}{Null effect, zero by default.}
\item{q_method}{Method for model-checking noncentrality; see \code{ba_meta}.}
}
\value{\code{williams_data()} returns a data frame with 23 rows and columns \code{class}, \code{study}, \code{estimate}, \code{variance}, \code{se}, \code{n_published_arms}, and \code{aggregation}. \code{williams_example()} returns a named list of two \code{bama_fit} objects. Each fit retains its treatment label so that individual summaries display Open-label TADs or PAT.}
\details{The data contain 15 open-label traditional-antidepressant (TAD) trials and 8 psychedelic-assisted-therapy (PAT) trials. The three-arm Arns et al. (2016) trial is excluded. Negative changes indicate improvement. Published posterior medians and 95\% credible intervals from Figure 2A of Williams et al. are treated as approximately Gaussian study summaries; standard errors were reconstructed as interval width divided by \code{2 * qnorm(0.975)}. These are within-arm changes, not treatment-versus-control causal contrasts. This example reproduces the BAMA illustrative reanalysis, not the original Bayesian analysis.}
\references{Williams, Z. J., Barnett, H., and Szigeti, B. (2026). Psychedelic therapy vs antidepressants for the treatment of depression under equal unblinding conditions: a systematic review and meta-analysis. JAMA Psychiatry.}
\examples{
dat <- williams_data()
fits <- williams_example()
bama_tables(fits)$panel_b
}
