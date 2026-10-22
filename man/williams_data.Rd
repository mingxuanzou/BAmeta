\name{williams_data}
\alias{williams_data}
\alias{williams_example}
\title{Reconstructed Williams Study Inputs}
\description{Load the reconstructed study summaries used in the BAMA manuscript as inputs to \code{bama}.}
\usage{
williams_data(group = NULL, format = c("list", "data.frame"))
williams_example(level = 0.95, null = 0, q_method = "exact")
}
\arguments{
\item{group}{Treatment group: \code{NULL} returns both groups; \code{"TAD"}
(or \code{"Open-label TADs"}) and \code{"PAT"} select one group.}
\item{format}{\code{"list"} returns study inputs ready for \code{bama};
\code{"data.frame"} returns the original data frame with provenance columns.}
\item{level}{Confidence level for the legacy \code{williams_example} wrapper,
strictly between 0.5 and 1.}
\item{null}{Null effect for the legacy wrapper, zero by default.}
\item{q_method}{Model-checking computation for the legacy wrapper; see \code{bama}.}
}
\value{
By default, \code{williams_data()} returns a named list with entries
\code{"Open-label TADs"} and \code{"PAT"}. Each entry is a list containing
numeric vectors \code{estimate} and \code{se}, a character vector \code{study},
and a treatment-label attribute. Selecting one group returns that group's
input list directly. These are data inputs, not fitted objects.

With \code{format = "data.frame"}, the function returns a data frame with
columns \code{class}, \code{study}, \code{estimate}, \code{variance},
\code{se}, \code{n_published_arms}, and \code{aggregation}. There are 23 rows
when both groups are selected. To fit this raw form, select one treatment
and the \code{estimate}, \code{se}, and \code{study} columns first.

\code{williams_example()} is retained for compatibility and returns a named
list of two \code{bama_fit} objects. New code should use
\code{bama(williams_data())}, setting fitting options in \code{bama}.
Each fitted group retains its treatment label for summaries and plots.
}
\details{The data contain 15 open-label traditional-antidepressant (TAD) trials and 8 psychedelic-assisted-therapy (PAT) trials. The three-arm Arns et al. (2016) trial is excluded. Negative changes indicate improvement. Published posterior medians and 95\% credible intervals from Figure 2A of Williams et al. are treated as approximately Gaussian study summaries; standard errors were reconstructed as interval width divided by \code{2 * qnorm(0.975)}. These are within-arm changes, not treatment-versus-control causal contrasts. This example reproduces the BAMA illustrative reanalysis, not the original Bayesian analysis.}
\references{Williams, Z. J., Barnett, H., and Szigeti, B. (2026). Psychedelic therapy vs antidepressants for the treatment of depression under equal unblinding conditions: a systematic review and meta-analysis. JAMA Psychiatry.}
\examples{
fit <- bama(williams_data())
summary(fit[["Open-label TADs"]])
summary(fit[["PAT"]])
bama_tables(fit)$panel_b

# Fit one treatment group.
pat <- bama(williams_data("PAT"))

# Inspect the input vectors or the original data frame.
williams_data("PAT")
raw <- williams_data(format = "data.frame")
}
\seealso{\code{\link{bama}}}
