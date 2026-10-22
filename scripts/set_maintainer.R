#!/usr/bin/env Rscript
# Run from the repository root with the intended PUBLIC maintainer email.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !grepl("^[^[:space:]<>@]+@[^[:space:]<>@]+\\.[^[:space:]<>@]+$", args[1L]))
  stop("Usage: Rscript scripts/set_maintainer.R your-public-email")
email <- args[1L]
d <- read.dcf("DESCRIPTION")
d[1L, "Authors@R"] <- paste0("person(\"Mingxuan\", \"Zou\", email = ",
                            encodeString(email, quote = '"'),
                            ", role = c(\"aut\", \"cre\"))")
d[1L, "Maintainer"] <- paste0("Mingxuan Zou <", email, ">")
write.dcf(d, file = "DESCRIPTION", width = 78)
cat("Updated public maintainer metadata.\n")
