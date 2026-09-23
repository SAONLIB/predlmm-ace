# Inference and interpretation

`ExactProfileML` optimizes variance proportions in unconstrained log-ratio
coordinates and profiles the total variance. For an interior ACE solution, the
observed profile-likelihood Hessian is inverted in log-ratio coordinates and
mapped to `(h2, c2, h2+c2)` with the full delta-method Jacobian.

The resulting standard errors condition on the selected rank and knots. They
do not account for rank selection or knot-selection uncertainty. Rank and knot
sensitivity should therefore be reported separately.

When `h2`, `c2`, or the residual proportion lies on or extremely near a
boundary, the regular quadratic approximation is unreliable. A missing or
invalid Hessian SE is an inferential diagnostic, not a failed point estimate.
Use a boundary-aware likelihood-ratio procedure, profile likelihood, or a
validated parametric bootstrap when boundary inference is scientifically
required.

The exact backend uses profile maximum likelihood, whereas GCTA-GREML uses
restricted maximum likelihood. Comparisons must identify this difference and
should use exactly matched participants and covariates whenever possible.
