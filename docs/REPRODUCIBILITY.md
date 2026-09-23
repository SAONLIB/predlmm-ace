# Reproducibility checklist

Record the following for every reported analysis:

1. genotype QC and GRM construction;
2. authoritative sample ID order;
3. phenotype and covariate definitions;
4. family/household label definition;
5. Nyström rank, knot-selection method, and random seed;
6. diagonal-correction diagnostics;
7. likelihood convention and optimizer convergence;
8. Hessian step and whether the solution is interior;
9. Slurm CPUs, memory request, elapsed time, and peak RSS;
10. software versions and exact repository commit.

For simulation coverage, use the primary 100-replicate run. Use the separate
20-replicate experiment only for paired rank and knot-strategy sensitivity.
