# Base-R graphics for the main-paper sensitivity and calibration displays.

.bama_plot_fits <- function(fits, labels, panel_labels) {
  if (inherits(fits, "bama_fit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), what = "bama_fit"))) {
    stop("fits must be a bama_fit or a nonempty list of bama_fit objects.",
         call. = FALSE)
  }
  n <- length(fits)
  if (is.null(labels)) {
    labels <- names(fits)
    if (is.null(labels)) labels <- rep("", n)
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- if (n == 1L) "Meta-analysis" else
      paste("Analysis", seq_len(n))[missing]
  }
  if (!is.character(labels) || length(labels) != n || anyNA(labels))
    stop("labels must contain one character label per fit.", call. = FALSE)
  if (is.null(panel_labels)) {
    panel_labels <- if (n == 1L) "" else if (n <= 26L)
      paste0("(", letters[seq_len(n)], ")") else paste0("(", seq_len(n), ")")
  }
  if (!is.character(panel_labels) || length(panel_labels) != n ||
      anyNA(panel_labels))
    stop("panel_labels must contain one character label per fit.", call. = FALSE)
  names(fits) <- make.unique(labels)
  list(fits = fits, labels = labels, panels = panel_labels)
}

.bama_plot_setup <- function(n, ncol, calibration = FALSE) {
  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol > .Machine$integer.max || ncol != as.integer(ncol))
    stop("ncol must be a positive integer.", call. = FALSE)
  nc <- min(ncol, n)
  graphics::par(mfrow = c(ceiling(n / nc), nc),
                mar = c(4.6, 4.5, 3.6, 1.1),
                oma = if (calibration) c(2.9, 0, 2.5, 0) else c(0, 0, 0, 0),
                mgp = c(2.8, 0.7, 0), tcl = -0.25, las = 1,
                family = "sans", cex = 0.94, cex.axis = 0.92,
                cex.lab = 1, bty = "l", xaxs = "i", yaxs = "i")
}

.bama_plot_grid <- function(b, n, maxima) {
  if (is.null(b)) b <- lapply(maxima, function(hi) seq(0, hi, length.out = 301L))
  else if (!is.list(b)) b <- rep(list(b), n)
  if (length(b) != n)
    stop("b must be a numeric grid or a list with one grid per fit.", call. = FALSE)
  lapply(b, function(x) {
    if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x)) || any(x < 0))
      stop("Each b grid needs at least two finite, nonnegative values.", call. = FALSE)
    x <- sort(unique(x))
    if (length(x) < 2L) stop("Each b grid needs two distinct values.", call. = FALSE)
    x
  })
}

.bama_plot_ranges <- function(xlim, n, defaults) {
  if (is.null(xlim)) return(defaults)
  if (!is.list(xlim)) xlim <- rep(list(xlim), n)
  if (length(xlim) != n)
    stop("xlim must be a numeric pair or a list with one pair per fit.", call. = FALSE)
  lapply(xlim, function(x) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) ||
        x[1L] < 0 || x[2L] <= x[1L])
      stop("Each xlim must be an increasing, finite, nonnegative pair.", call. = FALSE)
    x
  })
}

.bama_plot_title <- function(label, panel, subtitle) {
  title <- paste(c(panel[nzchar(panel)], label), collapse = " ")
  graphics::mtext(title, side = 3, line = 1.9, adj = 0, cex = 1.08)
  graphics::mtext(subtitle, side = 3, line = 0.45, adj = 0, cex = 0.92)
}

#' Plot BAMA sensitivity intervals
#' @export
bama_plot_sensitivity <- function(fits, b = NULL, labels = NULL,
                                  panel_labels = NULL, ncol = 2L,
                                  col = "#1F77B4", lwd = 2) {
  spec <- .bama_plot_fits(fits, labels, panel_labels)
  fits <- spec$fits
  n <- length(fits)
  maxima <- vapply(fits, function(x)
    max(1, 1.1 * x$thresholds$b_significance), numeric(1))
  grids <- .bama_plot_grid(b, n, maxima)
  curves <- Map(bama_ci, fits, grids)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  .bama_plot_setup(n, ncol)
  for (i in seq_len(n)) {
    fit <- fits[[i]]
    d <- curves[[i]]
    yr <- range(d$lower, d$upper, fit$null)
    pad <- max(diff(yr) * 0.04, .Machine$double.eps)
    graphics::plot(d$b, d$lower, type = "n", xlab = expression("Bias level " * b),
                   ylab = paste0(format(100 * fit$level, trim = TRUE),
                                 "% confidence interval"),
                   ylim = yr + c(-pad, pad), xlim = range(d$b), bty = "o")
    graphics::abline(h = fit$null, col = col, lty = 3, lwd = 1.2)
    graphics::abline(v = fit$thresholds$b_significance, col = col, lty = 3,
                     lwd = 1.2)
    graphics::lines(d$b, d$lower, col = col, lwd = lwd)
    graphics::lines(d$b, d$upper, col = col, lwd = lwd)
    subtitle <- bquote(bar(b)(.(formatC(1 - fit$level, digits = 3, format = "fg"))) ==
                        .(formatC(fit$thresholds$b_significance, digits = 3, format = "f")))
    .bama_plot_title(spec$labels[i], spec$panels[i], subtitle)
  }
  invisible(curves)
}

#' Plot BAMA model-checking p-values
#' @export
bama_plot_q <- function(fits, b = NULL, labels = NULL, panel_labels = NULL,
                        ncol = 2L, col = "#1F77B4", lwd = 2) {
  spec <- .bama_plot_fits(fits, labels, panel_labels)
  fits <- spec$fits
  n <- length(fits)
  maxima <- vapply(fits, function(x)
    max(1, 1.5 * x$thresholds$b_compatibility,
        1.25 * x$thresholds$approx_compatibility), numeric(1))
  grids <- .bama_plot_grid(b, n, maxima)
  curves <- Map(function(x, g) data.frame(b = g, p = bama_q_pvalue(x, g)),
                fits, grids)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  .bama_plot_setup(n, ncol)
  for (i in seq_len(n)) {
    fit <- fits[[i]]
    d <- curves[[i]]
    graphics::plot(d$b, d$p, type = "l", col = col, lwd = lwd,
                   xlab = expression("Bias level " * b), ylab = expression(p[Q](b)),
                   ylim = c(0, 1), xlim = range(d$b))
    graphics::abline(h = 1 - fit$level, v = fit$thresholds$b_compatibility,
                     lty = 3, col = "#AA4B35", lwd = 1.2)
    method <- if (identical(fit$q$method, "conservative")) " (conservative)" else ""
    subtitle <- bquote(underline(b)(.(formatC(1 - fit$level, digits = 3, format = "fg"))) ==
                        .(formatC(fit$thresholds$b_compatibility, digits = 3, format = "f")) *
                        .(method))
    .bama_plot_title(spec$labels[i], spec$panels[i], subtitle)
  }
  invisible(curves)
}

.bama_hist_breaks <- function(values, breaks, binwidth) {
  if (!is.null(binwidth)) {
    if (!is.numeric(binwidth) || length(binwidth) != 1L ||
        !is.finite(binwidth) || binwidth <= 0)
      stop("binwidth must be finite and strictly positive.", call. = FALSE)
    bins <- max(1, ceiling(max(values) / binwidth))
    if (!is.finite(bins) || bins > 1e5)
      stop("binwidth would require more than 100,000 bins.", call. = FALSE)
    return(seq(0, bins * binwidth, length.out = bins + 1L))
  }
  if (is.null(breaks)) breaks <- if (max(values) == 0) c(0, 1) else "FD"
  if (is.numeric(breaks) && length(breaks) > 1L &&
      (any(!is.finite(breaks)) || any(diff(breaks) <= 0) ||
       min(breaks) > min(values) || max(breaks) < max(values)))
    stop("Explicit histogram breaks must increase and include every b_j.", call. = FALSE)
  breaks
}

#' Plot descriptive standardized study differences
#' @export
bama_plot_calibration <- function(fits, breaks = NULL, labels = NULL,
                                  panel_labels = NULL, binwidth = NULL,
                                  xlim = NULL, ncol = 2L, shared_y = TRUE,
                                  col = "#557FA3", max_col = "#23384A",
                                  threshold_col = "#AA4B35") {
  spec <- .bama_plot_fits(fits, labels, panel_labels)
  fits <- spec$fits
  n <- length(fits)
  if (!is.null(breaks) && !is.null(binwidth))
    stop("Supply breaks or binwidth, not both.", call. = FALSE)
  if (!is.logical(shared_y) || length(shared_y) != 1L || is.na(shared_y))
    stop("shared_y must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(binwidth) && !length(binwidth) %in% c(1L, n))
    stop("binwidth must have length one or the number of fits.", call. = FALSE)
  if (!is.list(breaks)) breaks <- rep(list(breaks), n)
  if (length(breaks) != n)
    stop("breaks must be a common specification or a list with one per fit.", call. = FALSE)
  widths <- if (is.null(binwidth)) rep(list(NULL), n) else
    as.list(rep(binwidth, length.out = n))
  diagnostics <- lapply(fits, bama_diagnostics)
  histograms <- lapply(seq_len(n), function(i) {
    values <- diagnostics[[i]]$b_j
    graphics::hist(values,
                   breaks = .bama_hist_breaks(values, breaks[[i]], widths[[i]]),
                   plot = FALSE, include.lowest = TRUE, right = FALSE)
  })
  defaults <- lapply(seq_len(n), function(i) {
    hi <- max(histograms[[i]]$breaks, fits[[i]]$thresholds$b_significance)
    c(0, max(0.1, hi * 1.085))
  })
  ranges <- .bama_plot_ranges(xlim, n, defaults)
  ymax <- vapply(histograms, function(h) max(1, max(h$counts)), numeric(1))
  if (shared_y) ymax[] <- max(ymax)
  ymax <- ymax * 1.25
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  .bama_plot_setup(n, ncol, calibration = TRUE)
  for (i in seq_len(n)) {
    fit <- fits[[i]]
    values <- diagnostics[[i]]$b_j
    h <- histograms[[i]]
    maximum <- max(values)
    threshold <- fit$thresholds$b_significance
    graphics::plot(NA_real_, NA_real_, type = "n", xlim = ranges[[i]],
                   ylim = c(0, ymax[i]), axes = FALSE,
                   xlab = expression("Standardized absolute difference, " * b[j]),
                   ylab = "Number of studies")
    yticks <- pretty(c(0, ymax[i]), n = 5)
    yticks <- yticks[yticks >= 0 & yticks <= ymax[i] & yticks == floor(yticks)]
    graphics::abline(h = yticks, col = "#E2E6E9", lwd = 0.8)
    graphics::rect(utils::head(h$breaks, -1L), 0, utils::tail(h$breaks, -1L), h$counts,
                   col = col, border = "white", lwd = 0.8)
    graphics::axis(1)
    graphics::axis(2, at = yticks)
    graphics::box(bty = "l")
    graphics::segments(values, -0.025 * ymax[i], values, -0.055 * ymax[i],
                       col = max_col, xpd = NA)
    graphics::abline(v = maximum, col = max_col, lty = 3, lwd = 1.4)
    graphics::abline(v = threshold, col = threshold_col, lty = 2, lwd = 1.4)
    span <- diff(ranges[[i]])
    if (maximum >= ranges[[i]][1L] && maximum <= ranges[[i]][2L]) {
      right <- maximum > ranges[[i]][1L] + 0.64 * span
      graphics::text(maximum + if (right) -0.02 * span else 0.02 * span,
                     0.91 * ymax[i],
                     bquote(max[j]~b[j] == .(formatC(maximum, digits = 3, format = "f"))),
                     adj = c(as.integer(right), 0.5), col = max_col, cex = 0.89)
    }
    if (threshold >= ranges[[i]][1L] && threshold <= ranges[[i]][2L]) {
      right <- threshold > ranges[[i]][1L] + 0.5 * span
      graphics::text(threshold + if (right) -0.02 * span else 0.02 * span,
                     0.83 * ymax[i],
                     bquote(bar(b)(.(formatC(1 - fit$level, digits = 3, format = "fg"))) ==
                              .(formatC(threshold, digits = 3, format = "f"))),
                     adj = c(as.integer(right), 0.5), col = threshold_col, cex = 0.89)
    }
    .bama_plot_title(spec$labels[i], spec$panels[i], bquote(J == .(fit$J)))
  }
  graphics::mtext("Study differences are descriptive; they cannot detect a bias shared by all studies.",
                  side = 1, outer = TRUE, line = 1.1, cex = 0.78, col = "#303D49")
  # Draw a single legend in the outer top margin without opening another page.
  graphics::par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4), oma = rep(0, 4))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::legend(0.5, 0.99, xjust = 0.5, yjust = 1, horiz = TRUE, bty = "n",
                    legend = c("Study counts", "Observed maximum", "Significance threshold"),
                    col = c(col, max_col, threshold_col),
                    pch = c(15, NA_integer_, NA_integer_),
                    lty = c(0, 3, 2), lwd = c(0, 1.4, 1.4),
                    pt.cex = 1.3, cex = 0.84, x.intersp = 0.65)
  invisible(list(diagnostics = diagnostics, histograms = histograms))
}
