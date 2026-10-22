.bama_table_input <- function(fits, labels = NULL) {
  if (inherits(fits, "bama_fit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), "bama_fit"))) {
    stop("fits must be a bama_fit or a nonempty list of bama_fit objects.", call. = FALSE)
  }
  if (is.null(labels)) {
    labels <- names(fits)
    if (is.null(labels)) labels <- rep("", length(fits))
    blank <- is.na(labels) | !nzchar(labels)
    labels[blank] <- vapply(fits[blank], function(fit) {
      if (is.null(fit$label)) "" else fit$label
    }, character(1))
    blank <- is.na(labels) | !nzchar(labels)
    labels[blank] <- paste("Analysis", which(blank))
  }
  if (!is.character(labels) || length(labels) != length(fits) ||
      anyNA(labels) || any(!nzchar(labels))) {
    stop("labels must contain one nonempty character label per fit.", call. = FALSE)
  }
  list(fits = fits, labels = labels)
}

#' Tables corresponding to Panels A and B of the BAMA paper.
bama_tables <- function(fits, labels = NULL) {
  input <- .bama_table_input(fits, labels)
  panel_a <- panel_b <- heterogeneity <- vector("list", length(input$fits))
  for (i in seq_along(input$fits)) {
    x <- input$fits[[i]]
    label <- input$labels[[i]]
    methods <- list(FEMA = x$fe, REMA = x$re, FEMA_dagger = x$fe_dagger)
    panel_a[[i]] <- do.call(rbind, lapply(names(methods), function(method) {
      r <- methods[[method]]
      data.frame(treatment = label, procedure = method,
                 estimate = r$estimate, se = r$se, statistic = r$t,
                 lower = r$lower, upper = r$upper, p_value = r$p,
                 level = x$level, null = x$null, stringsAsFactors = FALSE)
    }))
    h <- x$thresholds
    panel_b[[i]] <- data.frame(
      treatment = label, optimistic_lower = h$optimistic_lower,
      optimistic_upper = h$optimistic_upper,
      optimistic_p_value = bama_pvalue(x, h$b_compatibility),
      b_significance = h$b_significance,
      approx_significance = h$approx_significance,
      b_compatibility = h$b_compatibility,
      approx_compatibility = h$approx_compatibility,
      q_method = x$q$method, level = x$level, null = x$null,
      optimistic_valid = FALSE, stringsAsFactors = FALSE)
    heterogeneity[[i]] <- data.frame(
      treatment = label, J = x$J, Q = x$Q, df = x$J - 1L,
      p_Q_0 = stats::pchisq(x$Q, df = x$J - 1L, lower.tail = FALSE),
      tau2 = x$re$tau2, sigma_B = sqrt(x$re$tau2),
      reference_scale = x$scale,
      variance_ratio = x$re$tau2 / x$scale^2,
      q_method = x$q$method, level = x$level, null = x$null,
      stringsAsFactors = FALSE)
  }
  out <- lapply(list(panel_a = panel_a, panel_b = panel_b,
                    heterogeneity = heterogeneity), function(rows) {
    d <- do.call(rbind, rows)
    rownames(d) <- NULL
    d
  })
  out
}

.bama_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 1L || digits > 10L || digits != floor(digits)) {
    stop("digits must be an integer between 1 and 10.", call. = FALSE)
  }
  as.integer(digits)
}

.bama_number <- function(x, digits) {
  # Avoid reporting negative zero after rounding.
  x[is.finite(x) & abs(x) < 0.5 * 10^(-digits)] <- 0
  formatC(x, format = "f", digits = digits)
}

.bama_p_text <- function(x, digits) {
  ifelse(x < 10^(-digits), paste0("<", .bama_number(10^(-digits), digits)),
         .bama_number(x, digits))
}

.bama_p_stars <- function(x, digits) {
  codes <- trimws(as.character(stats::symnum(
    x, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " "))))
  values <- .bama_p_text(x, digits)
  ifelse(nzchar(codes), paste(values, codes), values)
}

.bama_signif_legend <- function() {
  "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"
}

.bama_interval_text <- function(lower, upper, digits) {
  paste0("[", .bama_number(lower, digits), ", ", .bama_number(upper, digits), "]")
}

.bama_display_tables <- function(tables, digits) {
  a <- tables$panel_a
  b <- tables$panel_b
  levels <- unique(c(a$level, b$level))
  ci_label <- if (length(levels) == 1L) {
    paste0(format(100 * levels, trim = TRUE), "% CI")
  } else "CI"
  exact <- all(b$q_method == "exact")
  aout <- data.frame(
    Treatment = a$treatment, Procedure = a$procedure,
    Estimate = .bama_number(a$estimate, digits),
    SE = .bama_number(a$se, digits),
    `t-statistic` = .bama_number(a$statistic, digits),
    CI = .bama_interval_text(a$lower, a$upper, digits),
    `p-value` = .bama_p_stars(a$p_value, digits),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(aout)[6L] <- ci_label
  aout$Treatment[seq_len(nrow(aout)) %% 3L != 1L] <- ""
  bout <- data.frame(
    Treatment = b$treatment,
    Optimistic = .bama_interval_text(b$optimistic_lower, b$optimistic_upper, digits),
    `Optimistic p-value` = .bama_p_stars(b$optimistic_p_value, digits),
    `Significance: exact` = .bama_number(b$b_significance, digits),
    `Significance: approximation` = .bama_number(b$approx_significance, digits),
    Compatibility = .bama_number(b$b_compatibility, digits),
    `Compatibility: approximation` = .bama_number(b$approx_compatibility, digits),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(bout)[2L] <- "Optimistic interval"
  names(bout)[6L] <- if (exact) "Compatibility: exact" else "Compatibility: computed"
  list(panel_a = aout, panel_b = bout, ci_label = ci_label,
       compatibility_label = if (exact) "Exact" else "Computed")
}

.bama_escape <- function(x, format) {
  x <- as.character(x)
  if (format == "html") {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    return(gsub("'", "&#39;", x, fixed = TRUE))
  }
  if (format == "latex") {
    replacements <- c("\\" = "\\textbackslash{}", "{" = "\\{", "}" = "\\}",
                      "$" = "\\$", "&" = "\\&", "#" = "\\#", "%" = "\\%",
                      "_" = "\\_", "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}",
                      "<" = "\\textless{}", ">" = "\\textgreater{}")
    return(vapply(strsplit(x, "", fixed = TRUE), function(chars) {
      hit <- chars %in% names(replacements)
      chars[hit] <- unname(replacements[chars[hit]])
      paste(chars, collapse = "")
    }, character(1)))
  }
  for (mark in c("\\", "|", "*", "_", intToUtf8(96), "[", "]")) {
    x <- gsub(mark, paste0("\\", mark), x, fixed = TRUE)
  }
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\r\n|\r|\n", " ", x)
  x
}

.bama_notes <- function(tables, digits) {
  h <- tables$heterogeneity
  group_notes <- vapply(seq_len(nrow(h)), function(i) {
    paste0(h$treatment[i], ": J = ", h$J[i], "; Q = ", .bama_number(h$Q[i], digits),
           " on ", h$df[i], " degrees of freedom; p_Q(0) ",
           if (h$p_Q_0[i] < 10^(-digits)) "" else "= ", .bama_p_text(h$p_Q_0[i], digits),
           "; between-study variance = ", .bama_number(h$tau2[i], digits), ".")
  }, character(1))
  c(group_notes,
    "FEMA_dagger uses the FEMA point estimate with the REMA-adjusted variance.",
    if (any(h$q_method != "exact")) {
      "For analyses using the conservative Q calculation, the noncentrality is J b^2. Its compatibility threshold is conservative and can be smaller than the exact threshold."
    } else character(0))
}

.bama_html_table <- function(d, caption, top = NULL) {
  header <- paste0("<tr>", paste0("<th scope=\"col\">", .bama_escape(names(d), "html"),
                                   "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(d)), function(i) {
    cells <- .bama_escape(unlist(d[i, ], use.names = FALSE), "html")
    cells <- sub("FEMA_dagger", "FEMA<sub>&dagger;</sub>", cells, fixed = TRUE)
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<div class=\"table-scroll\"><table><caption>", caption,
         "</caption><thead>", paste(top, collapse = ""), header,
         "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table></div>")
}

.bama_tex_break <- function() paste(rep(intToUtf8(92), 2L), collapse = "")

.bama_latex_rows <- function(d, group_size = NULL) {
  vapply(seq_len(nrow(d)), function(i) {
    cells <- .bama_escape(unlist(d[i, ], use.names = FALSE), "latex")
    cells[cells == "FEMA\\_dagger"] <- "FEMA$_\\dagger$"
    space <- if (!is.null(group_size) && i > 1L &&
                 (i - 1L) %% group_size == 0L) "\\addlinespace\n" else ""
    paste0(space, paste(cells, collapse = " & "), " ", .bama_tex_break())
  }, character(1))
}

#' Render paper-style tables as HTML, LaTeX or Markdown.
bama_report <- function(fits, file = NULL, format = c("html", "latex", "markdown"),
                        labels = NULL, digits = 3) {
  format <- match.arg(format)
  digits <- .bama_digits(digits)
  if (!is.null(file) && (!is.character(file) || length(file) != 1L ||
                        is.na(file) || !nzchar(file))) {
    stop("file must be NULL or one nonempty path.", call. = FALSE)
  }
  tab <- bama_tables(fits, labels)
  display <- .bama_display_tables(tab, digits)
  a <- display$panel_a
  b <- display$panel_b
  notes <- .bama_notes(tab, digits)
  if (format == "html") {
    names(b)[4:7] <- c("Exact", "Approximation", display$compatibility_label, "Approximation")
    top <- paste0("<tr><th colspan=\"3\"></th><th colspan=\"2\">Maximum bias for significance</th>",
                  "<th colspan=\"2\">Minimum bias for compatibility</th></tr>")
    out <- paste(c(
      "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>BAMA analysis</title>",
      "<style>body{max-width:1160px;margin:36px auto;padding:0 24px;color:#17212b;background:#fff;font:16px/1.5 Georgia,serif}h1{font-size:24px;font-weight:normal}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:14px;margin:20px 0 12px}caption{text-align:left;font-weight:bold;font-size:17px;padding:10px 0}thead{border-top:2px solid #17212b;border-bottom:1px solid #17212b}tbody{border-bottom:2px solid #17212b}th,td{padding:7px 10px;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}td:nth-child(2),th:nth-child(2){text-align:left}thead tr:first-child:not(:last-child) th{white-space:normal;text-align:center}p.note{font-size:13px;margin:6px 0;color:#34414a}.table-scroll{overflow-x:auto}footer{margin-top:24px;border-top:1px solid #ccd2d6;padding-top:12px}</style></head><body>",
      "<h1>Bias-aware meta-analysis</h1>",
      .bama_html_table(a, "Panel A. FEMA and REMA"),
      paste0("<p class=\"note\">", .bama_escape(.bama_signif_legend(), "html"), "</p>"),
      paste0("<p class=\"note\">", .bama_escape(notes[seq_len(nrow(tab$heterogeneity))], "html"), "</p>"),
      .bama_html_table(b, "Panel B. BAMA", top),
      "<footer>",
      paste0("<p class=\"note\">", .bama_escape(notes[-seq_len(nrow(tab$heterogeneity))], "html"), "</p>"),
      "</footer></body></html>"), collapse = "\n")
  } else if (format == "latex") {
    names(b)[4:7] <- c("Exact", "Approximation", display$compatibility_label, "Approximation")
    out <- paste(c(
      "% Requires \\usepackage{booktabs}. This is a table fragment.",
      "\\begin{table}[htbp]", "\\centering", "\\small",
      "\\begin{tabular}{llrrrrr}",
      "\\multicolumn{7}{l}{\\textbf{Panel A. FEMA and REMA}}\\\\",
      "\\toprule",
      paste0(paste(.bama_escape(names(a), "latex"), collapse = " & "), " ", .bama_tex_break()),
      "\\midrule", .bama_latex_rows(a, group_size = 3L), "\\bottomrule", "\\end{tabular}",
      "\\par\\smallskip",
      paste0("{\\footnotesize\\raggedright ", .bama_escape(.bama_signif_legend(), "latex"), "\\par}"),
      paste0("{\\footnotesize\\raggedright ", .bama_escape(notes[seq_len(nrow(tab$heterogeneity))], "latex"), "\\par}"),
      "\\vspace{0.75em}", "\\setlength{\\tabcolsep}{5pt}",
      "\\begin{tabular}{lrrrrrr}",
      "\\multicolumn{7}{l}{\\textbf{Panel B. BAMA}}\\\\", "\\toprule",
      " & & & \\multicolumn{2}{c}{Maximum bias for significance} & \\multicolumn{2}{c}{Minimum bias for compatibility}\\\\",
      "\\cmidrule(lr){4-5}\\cmidrule(l){6-7}",
      paste0(paste(.bama_escape(names(b), "latex"), collapse = " & "), " ", .bama_tex_break()),
      "\\midrule", .bama_latex_rows(b), "\\bottomrule", "\\end{tabular}",
      "\\par\\smallskip",
      paste0("{\\footnotesize\\raggedright ", .bama_escape(notes[-seq_len(nrow(tab$heterogeneity))], "latex"), "\\par}"),
      "\\end{table}"), collapse = "\n")
  } else {
    md_table <- function(d) {
      row <- function(x) paste0("| ", paste(.bama_escape(x, "markdown"), collapse = " | "), " |")
      paste(c(row(names(d)), row(c(":---", rep("---:", ncol(d) - 1L))),
              apply(d, 1L, row)), collapse = "\n")
    }
    out <- paste("**Panel A. FEMA and REMA**", md_table(a),
                 .bama_escape(.bama_signif_legend(), "markdown"),
                 "**Panel B. BAMA**", md_table(b),
                 paste(.bama_escape(notes, "markdown"), collapse = "\n\n"), sep = "\n\n")
  }
  if (!is.null(file)) writeLines(enc2utf8(out), con = file, useBytes = TRUE)
  invisible(out)
}

summary.bama_fit <- function(object, ...) {
  out <- bama_tables(object)
  class(out) <- "summary.bama_fit"
  out
}

print.summary.bama_fit <- function(x, digits = 3, ...) {
  digits <- .bama_digits(digits)
  display <- .bama_display_tables(x, digits)
  cat("Panel A. FEMA and REMA\n")
  print(display$panel_a, row.names = FALSE, right = TRUE)
  cat("---\n", .bama_signif_legend(), "\n", sep = "")
  cat("\nPanel B. BAMA\n")
  print(display$panel_b, row.names = FALSE, right = TRUE)
  cat("\n", paste(.bama_notes(x, digits), collapse = "\n"), "\n", sep = "")
  invisible(x)
}

print.bama_fit <- function(x, digits = 3, ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
