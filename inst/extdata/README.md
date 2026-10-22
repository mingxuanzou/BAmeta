# Williams example provenance

`williams2026.csv` preserves the saved study-level reconstruction used in the
BAMA manuscript. There are 15 open-label TAD studies and 8 PAT studies;
the three-arm Arns et al. (2016) trial is excluded. The groups are analysed
separately. Negative HAM-D changes indicate improvement.

Source: Williams, Z. J., Barnett, H., and Szigeti, B. (2026),
*Psychedelic therapy vs antidepressants for the treatment of depression under
equal unblinding conditions: a systematic review and meta-analysis*,
JAMA Psychiatry, Figure 2A.

The displayed posterior medians are used as estimates. Standard errors were
approximated by the displayed 95% credible-interval width divided by
`2 * qnorm(0.975)`. These are model-based within-arm changes, not raw independent
treatment-versus-control effect estimates. Treating them as independent
Gaussian study summaries is an illustrative approximation. The package does
not claim to reproduce the original Bayesian model or infer a direct PAT vs
TAD treatment contrast.

The numeric file is preserved from the existing BAMA analysis; no study values
were simulated or newly digitised for this package.
