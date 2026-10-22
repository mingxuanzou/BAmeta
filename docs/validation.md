# Validation record

Validated on 21 October 2026 using R 4.3.3 on Linux, after renaming the package
to `BaMeta` and adding Mingxuan Zou <mingxuanzou@uchicago.edu> as maintainer.

- Independent Williams calculations reproduce the latest manuscript's rounded
  Panel A/B entries, both bias thresholds and the histogram maxima.
- Numerical regression tests check exact weighted partitions against exhaustive
  enumeration, confidence/p-value inversion, conservative Q ordering, unit and
  location invariance, explicit leave-one-out refits and invalid-input errors.
- All plot functions were rendered and visually inspected, including the
  all-zero standardized-residual case. Histogram counts are 15 and 8.
- HTML, Markdown and LaTeX table output was exercised; the LaTeX output was
  compiled and visually inspected.
- The installed example script generated CSV, HTML, LaTeX, Markdown, PNG and
  PDF outputs successfully.
- The final source tarball builds, installs and loads as `BaMeta`. Its installed
  metadata contains the real maintainer email supplied by the author.
- `R CMD check --as-cran` completed on the deliverable tarball with no errors
  or warnings and three notes: new submission, inability to verify current
  time, and skipped optional HTML validation/math rendering because `tidy`
  and the R package `V8` are unavailable in the local runtime.
- PDF manual compilation passed with `R_RD4PDF=times,hyper`, the documented
  R setting for using installed standard fonts when Inconsolata is unavailable.
- Both treatment summaries were rerun from the installed `BaMeta` package,
  including the optimistic p-values and significance codes. Their numerical
  results match the prior output.

Current R and R-devel checks on multiple operating systems remain pending.
The GitHub Actions workflow is configured for those checks. Their status is
reported at <https://github.com/mingxuanzou/ba_meta/actions>.

No annotated manuscript, private comments or unrelated source archive is
included. No CRAN submission has been made.
