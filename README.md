# BaMeta

**Bias-aware meta-analysis for independent study estimates.** This repository
implements the treatment-control methods in Mingxuan Zou and Peng Ding,
*Bias-aware meta-analysis: beyond fixed-effect and random-effects meta-analyses*
(2026 manuscript).

The R package and GitHub repository are named `BaMeta`. The main fitting
function is `bama()`.

Maintainer: Mingxuan Zou <mingxuanzou@uchicago.edu>.
See [release notes](docs/release.md) for publication and CRAN checks.

## Installation

Requires R >= 4.1.0. There are no third-party package dependencies and no
compiled code. Install from GitHub:

```r
install.packages("remotes")  # Needed once to install from GitHub
remotes::install_github("mingxuanzou/BaMeta")
library(BaMeta)
```

Alternatively, from the root of a downloaded or cloned repository:

```sh
R CMD INSTALL .
```

The publication preparation notes are in [docs/release.md](docs/release.md).

## Your own study summaries

Supply one estimate and its standard error (or variance) per independent study:

```r
library(BaMeta)

fit <- bama(
  estimate = c(-0.30, -0.50, -0.20, -0.40),
  se = c(0.10, 0.20, 0.15, 0.12),
  study = c("Study 1", "Study 2", "Study 3", "Study 4")
)
summary(fit)
bama_ci(fit, b = c(0, 0.5, 1, 2))
bama_thresholds(fit)
bama_diagnostics(fit)
```

`bama()` also accepts a named list or data frame containing `estimate`,
`se` (or `variance`), and optionally `study`. Use `variance = ...` instead
of `se = ...` when your input contains variances.
The default confidence level is 0.95 and the null value is zero. Invalid or
missing values are rejected; studies are never silently dropped. At least two
independent studies are required. Correlated multi-arm contrasts must not be
entered as if they were independent studies.

## Paper-style output

`williams_data()` returns a named list of the two treatment groups. Each group
contains `estimate`, `se`, and `study` vectors, in the same form as the inputs
above. Pass these data to `bama()` to fit both groups separately:

```r
fit <- bama(williams_data())
summary(fit[["Open-label TADs"]])
summary(fit[["PAT"]])

# Raw numeric tables, suitable for further analysis or CSV export
tables <- bama_tables(fit)
tables$panel_a
tables$panel_b

# Manuscript-style Panels A and B
bama_report(fit, file = "results.html")
bama_report(fit, file = "results.tex", format = "latex")

# Figures
bama_plot_sensitivity(fit)
bama_plot_calibration(fit, binwidth = c(1, 0.5))
bama_plot_q(fit)
```

To fit one treatment group, use:

```r
fit_pat <- bama(williams_data("PAT"))
summary(fit_pat)
```

You can inspect or edit the input vectors before fitting:

```r
dat <- williams_data("PAT")
dat$estimate
dat$se
dat$study
fit_pat <- bama(dat)
```

For the complete source data and provenance columns, use
`williams_data(format = "data.frame")`.

**Panel A** reports FEMA, REMA (DerSimonian-Laird), and the FEMA estimator with
REMA-adjusted variance (FEMA-dagger): estimate, standard error, normal-reference
test statistic, confidence interval, and p-value. Its notes report study counts,
Cochran's Q, degrees of freedom, Q p-values, and between-study variance.

**Panel B** reports the optimistic interval and optimistic p-value, the maximum bias for significance,
and the minimum bias for compatibility, alongside their first-order
approximations. Exact and conservative Q calculations are explicitly identified.

The optimistic p-value evaluates `bama_pvalue()` at the minimum bias for
compatibility, using the fitted null value and model-checking level. It is
the smallest target-test p-value among the compatible bias levels. Like the
optimistic interval, it is descriptive and has no general type I error
guarantee, because compatibility does not establish an upper bias bound.
The numeric value is returned as `bama_tables(fit)$panel_b$optimistic_p_value`.

![Panels A and B](docs/figures/panels.png)

Sensitivity curves have solid blue confidence bounds and dotted null and
significance-threshold lines. Calibration histograms use study counts, individual
rug marks, a dotted observed maximum, and a dashed significance threshold.

![Sensitivity curves](docs/figures/sensitivity.png)

![Calibration histograms](docs/figures/calibration.png)

PDF output uses standard R graphics devices:

```r
pdf("calibration.pdf", width = 10.2, height = 4.65)
bama_plot_calibration(fit, binwidth = c(1, 0.5))
dev.off()
```

The HTML report is self-contained. The LaTeX output is a fragment requiring
`\usepackage{booktabs}`. Every plotting function returns its underlying data
invisibly and restores the previous graphical settings.

## Reproduce the paper example

After installation, from the repository root:

```sh
Rscript scripts/reproduce_williams.R outputs
```

This writes Panels A/B in HTML, LaTeX, Markdown and CSV; sensitivity, Q and
calibration figures in PDF and PNG; study-level calibration/leave-one-out
diagnostics; curve values; and R session information. The same script is
available after installation under
`system.file("examples", "reproduce_williams.R", package = "BaMeta")`.

The included data reproduce these reference results:

| Group | Studies | FEMA estimate | Maximum bias for significance | Minimum bias for compatibility | Maximum observed b_j |
|---|---:|---:|---:|---:|---:|
| Open-label TADs | 15 | -13.078 | 20.475 | 3.358 | 7.476 |
| PAT | 8 | -11.267 | 5.671 | 0.270 | 2.312 |

The Williams data are saved reconstructions of posterior medians and 95%
credible intervals from Figure 2A, with standard errors approximated from
interval widths. They represent within-arm HAM-D changes, not direct
treatment-control causal effects. Arns et al. (2016) is excluded. See
[data provenance](inst/extdata/README.md).

## Citation and license

Please cite Zou and Ding's BAMA manuscript and the package version used.
`citation("BaMeta")` prints the software citation. Code is distributed under
the [MIT license](LICENSE.md). The data documentation retains attribution to
the source publication.
