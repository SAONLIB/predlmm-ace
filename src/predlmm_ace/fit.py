#!/usr/bin/env python3
"""Fit PredLMM-ACE/AE to an aligned ABCD phenotype and metadata table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from predlmm_ace.model import ExactProfileML, LowRankACE


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--factor-dir", required=True)
    p.add_argument("--metadata", required=True, help="TSV with FID, IID and analysis variables")
    p.add_argument("--phenotype", required=True)
    p.add_argument("--family-column", required=True)
    p.add_argument("--continuous", default="", help="Comma-separated numeric covariates")
    p.add_argument("--categorical", default="", help="Comma-separated categorical covariates")
    p.add_argument("--model", choices=["ACE", "AE"], default="ACE")
    p.add_argument("--keep", default=None, help="Optional headerless PLINK FID/IID keep file")
    p.add_argument("--standard-errors", action="store_true",
                   help="Compute observed profile-Hessian/delta-method standard errors")
    p.add_argument("--hessian-step", type=float, default=2e-3,
                   help="Finite-difference step on the log-ratio scale")
    p.add_argument("--out", required=True)
    args = p.parse_args()
    factor_dir = Path(args.factor_dir)
    ids = pd.read_csv(factor_dir / "sample_ids.tsv", sep="\t", dtype=str)
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    if meta[["FID", "IID"]].duplicated().any():
        raise ValueError("Duplicate FID/IID in analysis metadata.")
    data = ids.merge(meta, on=["FID", "IID"], how="left", validate="one_to_one", indicator=True)
    if (data["_merge"] != "both").any():
        raise ValueError(f"Metadata absent for {(data['_merge'] != 'both').sum()} GRM samples.")
    data = data.drop(columns="_merge")
    analysis_sample = "full_factor_sample"
    if args.keep:
        keep = pd.read_csv(args.keep, sep=r"\s+", header=None, dtype=str)
        if keep.shape[1] == 1:
            keep.columns = ["IID"]
            mask = data["IID"].isin(set(keep["IID"]))
        else:
            keep = keep.iloc[:, :2]
            keep.columns = ["FID", "IID"]
            wanted = set(map(tuple, keep[["FID", "IID"]].to_numpy()))
            mask = np.fromiter(
                ((fid, iid) in wanted for fid, iid in data[["FID", "IID"]].to_numpy()),
                dtype=bool, count=len(data)
            )
        if int(np.sum(mask)) != len(keep):
            raise ValueError(
                f"Keep file has {len(keep)} rows but {int(np.sum(mask))} matched factor samples."
            )
        data = data.loc[mask].reset_index(drop=True)
        analysis_sample = "keep_subset"
    else:
        mask = np.ones(len(data), dtype=bool)
    continuous = [x for x in args.continuous.split(",") if x]
    categorical = [x for x in args.categorical.split(",") if x]
    needed = [args.phenotype, args.family_column] + continuous + categorical
    missing = [x for x in needed if x not in data]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if data[needed].isna().any().any():
        counts = data[needed].isna().sum()
        raise ValueError(f"Missing analysis values; rebuild analysis.keep first:\n{counts[counts > 0]}")
    y = pd.to_numeric(data[args.phenotype], errors="raise").to_numpy(float)
    design = pd.DataFrame({"intercept": np.ones(len(data))})
    for column in continuous:
        values = pd.to_numeric(data[column], errors="raise").to_numpy(float)
        sd = values.std(ddof=0)
        if sd == 0:
            raise ValueError(f"Continuous covariate {column} has zero variance.")
        design[column] = (values - values.mean()) / sd
    if categorical:
        dummies = pd.get_dummies(data[categorical].astype(str), drop_first=True, dtype=float)
        design = pd.concat([design, dummies], axis=1)
    x = design.to_numpy(float)
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("Covariate design is rank deficient; revise categorical/reference levels.")
    labels = data[args.family_column].astype(str).to_numpy()
    sizes = data.groupby(args.family_column).size()
    if args.model == "ACE" and int((sizes > 1).sum()) == 0:
        raise ValueError("ACE is not identifiable: every family is a singleton. Fit --model AE.")
    factor_all = np.load(factor_dir / "nystrom_factor.npy", mmap_mode="r")
    delta_all = np.load(factor_dir / "diagonal_correction.npy")
    factor = factor_all[mask]
    delta = delta_all[mask]
    covariance = LowRankACE(factor, delta, labels)
    result = ExactProfileML(covariance).fit(
        y, x, fix_c2_zero=args.model == "AE",
        compute_standard_errors=args.standard_errors,
        hessian_step=args.hessian_step,
    )
    payload = {
        "model": args.model, "n": len(data), "rank": int(factor.shape[1]),
        "analysis_sample": analysis_sample,
        "keep_file": args.keep,
        "phenotype": args.phenotype, "family_column": args.family_column,
        "n_families": int(sizes.size), "n_multi_participant_families": int((sizes > 1).sum()),
        "max_family_size": int(sizes.max()), "design_columns": list(design.columns),
        "genetic_variance": result.genetic_variance,
        "shared_environment_variance": result.shared_environment_variance,
        "residual_variance": result.residual_variance,
        "h2": result.h2, "c2": result.c2, "beta": result.beta.tolist(),
        "converged": result.converged, "iterations": result.iterations,
        "profile_objective": result.objective, "message": result.message,
        "se_h2": result.se_h2, "se_c2": result.se_c2,
        "se_h2_plus_c2": result.se_h2_plus_c2,
        "covariance_proportions": (
            None if result.covariance_proportions is None
            else result.covariance_proportions.tolist()
        ),
        "hessian_logratio": (
            None if result.hessian_logratio is None else result.hessian_logratio.tolist()
        ),
        "hessian_condition_number": result.hessian_condition_number,
        "hessian_min_eigenvalue": result.hessian_min_eigenvalue,
        "se_method": result.se_message,
        "hessian_step": args.hessian_step if args.standard_errors else None,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n")
    Path(args.out).with_suffix(".design.tsv").write_text("\n".join(design.columns) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
