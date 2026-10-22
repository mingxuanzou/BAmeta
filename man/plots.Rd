\name{bama_plots}
\alias{bama_plot_sensitivity}
\alias{bama_plot_q}
\alias{bama_plot_calibration}
\title{Paper-Style BAMA Sensitivity and Calibration Plots}
\description{
Draw sensitivity confidence limits, model-checking p-values, or descriptive
histograms of standardized study differences using base R graphics.
The defaults follow the displays accompanying the BAMA paper.
}
\usage{
bama_plot_sensitivity(fits, b = NULL, labels = NULL,
  panel_labels = NULL, ncol = 2L, col = "#1F77B4", lwd = 2)

bama_plot_q(fits, b = NULL, labels = NULL, panel_labels = NULL,
  ncol = 2L, col = "#1F77B4", lwd = 2)

bama_plot_calibration(fits, breaks = NULL, labels = NULL,
  panel_labels = NULL, binwidth = NULL, xlim = NULL,
  ncol = 2L, shared_y = TRUE, col = "#557FA3",
  max_col = "#23384A", threshold_col = "#AA4B35")
}
\arguments{
\item{fits}{A \code{bama_fit} returned by \code{bama}, or a nonempty
list of fits. Each fit receives its own panel.}
\item{b}{A finite nonnegative numeric grid with at least two distinct values,
or a list containing one such grid per fit. The default grids include the
relevant threshold. Supplied grids are sorted and deduplicated.}
\item{labels}{One character title per fit. Defaults to the list names, or
generic analysis labels when names are absent.}
\item{panel_labels}{One character panel label per fit. The default is
\code{(a)}, \code{(b)}, and so on for multiple fits; the single-fit default
is an empty label. Supply \code{c("A", "B")} for upper-case labels.}
\item{ncol}{Maximum number of panels per row. Use a sufficiently wide
graphics device, for example width 10.2 and height 4.65 inches for two panels.}
\item{col}{Line color for sensitivity and model-checking plots, or bar
color for histograms.}
\item{lwd}{Line width for sensitivity and model-checking curves.}
\item{breaks}{A \code{hist} break specification common to every histogram,
or a list with one specification per fit. Defaults to \code{"FD"}.
Explicit break vectors must increase and include all observations.}
\item{binwidth}{An alternative to \code{breaks}: a positive width common
to every histogram or one width per fit. Bins start at zero. The Williams
example uses \code{c(1, 0.5)} for TAD and PAT.}
\item{xlim}{Histogram x limits, either an increasing nonnegative numeric
pair or a list with one pair per fit. The default includes every histogram
bin and the maximum bias for significance. Explicit limits may hide a line.}
\item{shared_y}{Use a common count-axis range across histograms?}
\item{max_col}{Color for the observed maximum and individual-study rug marks.}
\item{threshold_col}{Color for the significance threshold in histograms.}
}
\details{
Sensitivity panels show both confidence limits as solid curves, the null
value as a dotted horizontal line, and the maximum bias for significance
as a dotted vertical line. A maximum bias of zero means that the FEMA
interval already contains the null value.

Model-checking panels show \eqn{p_Q(b)}, the significance level as a dotted
horizontal line, and the minimum bias for compatibility as a dotted vertical
line. A conservative noncentrality bound is explicitly labeled when used.

Calibration histograms count independent studies without precision weighting.
They plot \eqn{b_j = |\hat\tau_j - \hat\tau_{FE}| /
(\sqrt{J}\sigma_{FE})}, with rug marks for individual studies, a dotted
observed maximum, and a dashed significance threshold. The horizontal scales
may differ across panels. These diagnostics describe dispersion around the
pooled estimate; they do not estimate an upper bound on bias and cannot
detect common bias. The actual leave-one-study-out pooled shift is separately
returned by \code{bama_diagnostics}.

The functions draw on the current device and restore its graphical parameters
on exit. Use \code{pdf()} or \code{png()} before calling a plotting function
and \code{dev.off()} afterward to save a figure.
}
\value{
Invisibly, a named list of plotted data frames for sensitivity and p-value
plots. Calibration returns a list with named \code{diagnostics} data frames
and \code{histograms} containing the break points and unweighted counts.
}
\examples{
fit <- bama(c(-1.2, -0.7, -1.6, -0.9), se = c(0.3, 0.4, 0.5, 0.3))
bama_plot_sensitivity(fit)
bama_plot_q(fit)
bama_plot_calibration(fit, binwidth = 0.5)
}
\seealso{\code{\link{bama}}, \code{\link{bama_diagnostics}}}
