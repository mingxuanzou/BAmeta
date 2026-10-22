# BAmeta 0.1.2

- Uses `BAmeta` consistently as the R package and GitHub repository name.
- Updates installation commands, package data lookup, documentation, citation,
  examples, tests and automated checks for the renamed package.
- Retains the `bama()` interface and the existing numerical methods.

# BAmeta 0.1.1

- Makes `bama()` the main fitting function, with direct vectors, named study
  inputs, or a named list of inputs for separate analyses.
- Changes `williams_data()` to return study inputs containing `estimate`, `se`,
  and `study`, ready for `bama(williams_data())`. Select one treatment with
  `williams_data("PAT")` or `williams_data("TAD")`; use
  `format = "data.frame"` for the original data frame with provenance columns.
- Retains `ba_meta()` and `williams_example()` as compatibility functions.
- Renames the GitHub repository and simplifies the README.

# BAmeta 0.1.0

- Implements the treatment-control BAMA methods in the October 2026 manuscript.
- Adds FEMA, REMA and FEMA-dagger summaries and manuscript-style Panels A and B.
- Adds exact weighted-partition Q calibration for moderate numbers of studies,
  with an explicitly labelled conservative method for larger analyses.
- Adds sensitivity curves, Q curves, calibration histograms and distinct
  leave-one-study-out diagnostics.
- Includes saved reconstructed Williams data and a reproducible example.

This package has a separate treatment-control interface from the earlier
network-focused BiasRobustMeta package.
