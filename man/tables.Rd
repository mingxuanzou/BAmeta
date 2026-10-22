\name{bama_tables}
\alias{bama_tables}
\alias{bama_report}
\alias{summary.bama_fit}
\alias{print.bama_fit}
\alias{print.summary.bama_fit}
\title{Paper-Style BAMA Tables and Reports}
\description{
Construct the numerical and formatted counterparts of Panels A and B of the
BAMA paper. Panel A reports FEMA, REMA, and FEMA with REMA-adjusted
variance. Panel B reports the two bias thresholds, their first-order
approximations, and the descriptive optimistic interval and p-value.
}
\usage{
bama_tables(fits, labels = NULL)

bama_report(fits, file = NULL,
            format = c("html", "latex", "markdown"),
            labels = NULL, digits = 3)

\method{summary}{bama_fit}(object, ...)
\method{print}{bama_fit}(x, digits = 3, ...)
\method{print}{summary.bama_fit}(x, digits = 3, ...)
}
\arguments{
  \item{fits}{A \code{bama_fit}, or a nonempty named list of such fits.}
  \item{labels}{Optional nonempty character labels, one per fit. By default,
  list names are used, followed by stored treatment labels if available;
  remaining unnamed entries receive labels such as \code{"Analysis 1"}.}
  \item{file}{Optional output path. The containing directory must exist.}
  \item{format}{Output format. HTML is a complete self-contained document.
  LaTeX is a table fragment requiring \code{booktabs} to be loaded in
  the containing document. Markdown is plain Markdown text.}
  \item{digits}{Number of digits after the decimal point, from 1 to 10.
  Smaller p-values are displayed as less than \eqn{10^{-\mathrm{digits}}}.}
  \item{object, x}{A fitted object or its summary, as appropriate.}
  \item{...}{Additional arguments. Currently unused by \code{summary};
  forwarded by the fitted-object print method to the summary print method.}
}
\value{
\code{bama_tables} returns a list with three data frames:
\describe{
  \item{\code{panel_a}}{One row per procedure and analysis, with numeric
  estimates, standard errors, test statistics, interval endpoints, p-values,
  confidence levels and null values.}
  \item{\code{panel_b}}{One row per analysis, with optimistic interval endpoints,
  \code{optimistic_p_value},
  significance and compatibility thresholds and approximations, the Q
  calculation method, confidence level, null value, and
  \code{optimistic_valid = FALSE}.}
  \item{\code{heterogeneity}}{The number of studies, Q statistic, degrees of
  freedom, central Q-test p-value, DerSimonian--Laird between-study variance
  and standard deviation, common reference scale, between-study to
  harmonic-mean within-study variance ratio, and Q calculation method.}
}
\code{bama_report} returns the rendered text invisibly and, if \code{file} is
supplied, writes that same text as UTF-8. \code{summary} returns the same
numerical tables with class \code{summary.bama_fit}. Print methods return
their input invisibly.
}
\details{
P-values display significance codes, with a legend below Panel A,
using the conventions of \code{stats::printCoefmat}. Codes are based
on the unrounded p-values; numerical tables retain numeric p-values.

FEMA with REMA-adjusted variance, denoted FEMA-dagger in the paper, retains the
FEMA estimate and adds estimated between-study variation to its variance.
REMA uses DerSimonian--Laird estimation and a normal critical value.
Compatibility results are labelled exact only when the fit used the exact
partition calculation. Conservative calculations use noncentrality
\eqn{Jb^2}; their thresholds can be smaller than the exact thresholds.

The optimistic interval uses the estimated minimum bias for compatibility.
The optimistic p-value evaluates \code{bama_pvalue} at the same threshold,
using the fitted null value. Because that p-value is nondecreasing in the
bias level, it is the smallest target-test p-value among bias levels not
rejected by the Q-test. It depends on the model-checking level
\eqn{1-\mathrm{level}}. A value below that level is equivalent to exclusion
of the fitted null by the optimistic interval.
Both quantities are descriptive: the interval need not attain nominal
coverage, and the optimistic p-value does not generally control type I
error. Significance codes do not provide such a guarantee. The two bias
thresholds are not confidence limits for the true bias. When comparing fits
with different levels or null values, consult the \code{level} and \code{null}
columns returned by \code{bama_tables}.
Labels are escaped for the selected output format.
}
\examples{
x <- ba_meta(c(-0.4, -0.7, -0.2, -0.5), se = c(0.1, 0.2, 0.15, 0.12))
tab <- bama_tables(x)
tab$panel_a
tab$panel_b
cat(bama_report(x, format = "markdown"))
html <- tempfile(fileext = ".html")
bama_report(x, file = html)
unlink(html)
}
\seealso{\code{\link{ba_meta}}, \code{\link{bama_diagnostics}}}
