#' Load the reconstructed Williams example
#' @export
williams_data <- function() {
  utils::read.csv(system.file("extdata", "williams2026.csv", package = "BaMeta"),
                  stringsAsFactors = FALSE, check.names = FALSE)
}

#' Fit the two treatment groups in the paper's example
#' @export
williams_example <- function(level = 0.95, null = 0, q_method = "exact") {
  dat <- williams_data()
  fits <- lapply(c("TAD", "PAT"), function(group) {
    d <- dat[dat$class == group, , drop = FALSE]
    ba_meta(d$estimate, variance = d$variance, study = d$study,
            level = level, null = null, q_method = q_method)
  })
  names(fits) <- c("Open-label TADs", "PAT")
  for (i in seq_along(fits)) fits[[i]]$label <- names(fits)[i]
  fits
}
