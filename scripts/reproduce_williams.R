#!/usr/bin/env Rscript
library(BaMeta)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "outputs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fits <- bama(williams_data())
tables <- bama_tables(fits)
for (nm in names(tables))
  write.csv(tables[[nm]], file.path(out, paste0(nm, ".csv")), row.names = FALSE)
for (fmt in c("html", "latex", "markdown")) {
  extension <- switch(fmt, html = "html", latex = "tex", markdown = "md")
  bama_report(fits, file.path(out, paste0("panels.", extension)), format = fmt)
}

render <- function(stem, draw) {
  pdf(file.path(out, paste0(stem, ".pdf")), width = 10.2, height = 4.65,
      useDingbats = FALSE)
  tryCatch(draw(), finally = dev.off())
  png(file.path(out, paste0(stem, ".png")), width = 2448, height = 1116,
      res = 240)
  tryCatch(draw(), finally = dev.off())
}
grids <- list(seq(0, 22.2, length.out = 301), seq(0, 6.2, length.out = 301))
render("sensitivity", function() bama_plot_sensitivity(fits, b = grids))
render("calibration", function() bama_plot_calibration(
  fits, binwidth = c(1, 0.5), xlim = list(c(0, 22.2), c(0, 6.2))))
render("q_model_check", function() bama_plot_q(fits))

diagnostics <- do.call(rbind, lapply(seq_along(fits), function(i) {
  data.frame(group = names(fits)[i], bama_diagnostics(fits[[i]]))
}))
write.csv(diagnostics, file.path(out, "study_diagnostics.csv"), row.names = FALSE)
curves <- do.call(rbind, lapply(seq_along(fits), function(i) {
  data.frame(group = names(fits)[i], bama_ci(fits[[i]], grids[[i]]),
             target_p = bama_pvalue(fits[[i]], grids[[i]]),
             q_p = bama_q_pvalue(fits[[i]], grids[[i]]))
}))
write.csv(curves, file.path(out, "sensitivity_values.csv"), row.names = FALSE)
write.csv(williams_data(format = "data.frame"), file.path(out, "input_data.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Reproduced Williams tables and figures in", normalizePath(out), "\n")
