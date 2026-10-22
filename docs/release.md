# Publication preparation

The installable R package is `BaMeta`, version 0.1.0. The maintainer is
Mingxuan Zou <mingxuanzou@uchicago.edu>. Build, check and install the source
package from the repository root:

```sh
R CMD build .
R CMD check --as-cran BaMeta_0.1.0.tar.gz
R CMD INSTALL BaMeta_0.1.0.tar.gz
Rscript scripts/reproduce_williams.R outputs
```

The public repository is <https://github.com/mingxuanzou/ba_meta>, retaining
the installable name `BaMeta`. The DESCRIPTION and README contain the repository
links and GitHub installation command. Private annotated manuscripts and
unrelated source archives are excluded.

The prepared GitHub Actions workflow runs CRAN-style checks on current R
for Linux, Windows and macOS, and on R-devel for Linux. It checks installation,
documentation, examples and numerical regression tests, excluding PDF manual
compilation. Check the actual run results at
<https://github.com/mingxuanzou/ba_meta/actions> before CRAN submission. A public
v0.1.0 tag/release should identify the exact checked commit used for the paper.

Before CRAN submission, run checks with current R and R-devel and on at least
two major operating systems. Upload the built source tarball through
<https://cran.r-project.org/submit.html>. The maintainer must accept CRAN's
confirmation email. No CRAN submission has been made.
