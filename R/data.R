#' Load the reconstructed Williams example
#' @export
williams_data <- function(group = NULL, format = c("list", "data.frame")) {
  format <- match.arg(format)
  groups <- c("TAD", "PAT")
  labels <- c("Open-label TADs", "PAT")
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || is.na(group) ||
        !group %in% c("TAD", "Open-label TADs", "PAT")) {
      stop("group must be NULL, 'TAD', 'Open-label TADs', or 'PAT'.",
           call. = FALSE)
    }
    index <- if (group == "PAT") 2L else 1L
    groups <- groups[index]
    labels <- labels[index]
  }
  dat <- utils::read.csv(
    system.file("extdata", "williams2026.csv", package = "BAmeta"),
    stringsAsFactors = FALSE, check.names = FALSE)
  dat <- dat[dat$class %in% groups, , drop = FALSE]
  if (format == "data.frame") return(dat)
  inputs <- lapply(seq_along(groups), function(i) {
    d <- dat[dat$class == groups[i], , drop = FALSE]
    input <- list(estimate = d$estimate, se = d$se, study = d$study)
    attr(input, "label") <- labels[i]
    input
  })
  names(inputs) <- labels
  if (is.null(group)) inputs else inputs[[1L]]
}

#' Fit the two treatment groups in the paper's example
#' @export
williams_example <- function(level = 0.95, null = 0, q_method = "exact") {
  bama(williams_data(), level = level, null = null, q_method = q_method)
}
