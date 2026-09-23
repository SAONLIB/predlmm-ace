#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)

    factors = pd.read_csv(root / "manifests/factors.tsv", sep="\t")
    rows = []
    for _, factor in factors.iterrows():
        fid = int(factor.factor_id)
        base = root / f"results/predlmm_{fid:06d}.json"
        se_path = root / f"results/predlmm_{fid:06d}_se.json"
        path = se_path if se_path.exists() else base
        if not path.exists():
            continue
        result = json.loads(path.read_text())
        rows.append(
            {
                **factor.to_dict(),
                "method": "PredLMM-ACE",
                "estimate_h2": result.get("h2"),
                "estimate_c2": result.get("c2"),
                "se_h2": result.get("se_h2"),
                "se_c2": result.get("se_c2"),
                "fit_runtime_sec": result.get("fit_runtime_sec"),
                "runtime_sec": result.get("runtime_sec"),
                "converged": result.get("converged", False),
                "se_available": se_path.exists() and result.get("se_h2") is not None,
            }
        )

    replicates = pd.read_csv(root / "manifests/replicates.tsv", sep="\t")
    for _, replicate in replicates.iterrows():
        rid = int(replicate.task_id)
        external_dir = root / f"external/task_{rid:06d}"
        for filename in ("gcta_results.csv", "fastgwa_results.csv"):
            result_path = external_dir / filename
            if not result_path.exists():
                continue
            external = pd.read_csv(result_path)
            for _, result in external.iterrows():
                external_row = {**replicate.to_dict(), **result.to_dict()}
                external_row["converged"] = bool(result.get("success", False))
                rows.append(external_row)

    data = pd.DataFrame(rows)
    data.to_csv(root / "all_results.csv", index=False)
    if data.empty:
        raise RuntimeError("No PredLMM result files were found")

    for column in ("knot_method", "rank", "knot_fraction", "sparse_cutoff"):
        if column not in data:
            data[column] = np.nan

    group_keys = [
        "simulation", "N", "h2", "c2", "method",
        "knot_method", "rank", "knot_fraction", "sparse_cutoff",
    ]
    summaries = []
    for key, frame in data.groupby(group_keys, dropna=False):
        estimate_h2 = numeric(frame.estimate_h2)
        estimate_c2 = numeric(frame["estimate_c2"]) if "estimate_c2" in frame else pd.Series(np.nan, index=frame.index)
        se_h2 = numeric(frame["se_h2"]) if "se_h2" in frame else pd.Series(np.nan, index=frame.index)
        se_c2 = numeric(frame["se_c2"]) if "se_c2" in frame else pd.Series(np.nan, index=frame.index)
        truth_h2 = float(frame.h2.iloc[0])
        truth_c2 = float(frame.c2.iloc[0])
        valid_h2_se = estimate_h2.notna() & se_h2.notna() & (se_h2 > 0)
        valid_c2_se = estimate_c2.notna() & se_c2.notna() & (se_c2 > 0)

        coverage_h2 = np.nan
        if valid_h2_se.any():
            coverage_h2 = (
                (estimate_h2[valid_h2_se] - 1.96 * se_h2[valid_h2_se] <= truth_h2)
                & (truth_h2 <= estimate_h2[valid_h2_se] + 1.96 * se_h2[valid_h2_se])
            ).mean()

        coverage_c2 = np.nan
        # Regular log-ratio Wald intervals are not valid when true c2 is on
        # the boundary; report boundary frequency instead.
        if truth_c2 > 0 and valid_c2_se.any():
            coverage_c2 = (
                (estimate_c2[valid_c2_se] - 1.96 * se_c2[valid_c2_se] <= truth_c2)
                & (truth_c2 <= estimate_c2[valid_c2_se] + 1.96 * se_c2[valid_c2_se])
            ).mean()

        summaries.append(
            {
                **dict(zip(group_keys, key)),
                "n_expected": len(frame),
                "n_h2": int(estimate_h2.notna().sum()),
                "n_c2": int(estimate_c2.notna().sum()),
                "convergence_rate": (
                    pd.Series(frame["converged"]).fillna(False).astype(bool).mean()
                    if "converged" in frame else np.nan
                ),
                "mean_h2": estimate_h2.mean(),
                "bias_h2": (estimate_h2 - truth_h2).mean(),
                "rmse_h2": np.sqrt(((estimate_h2 - truth_h2) ** 2).mean()),
                "empirical_sd_h2": estimate_h2.std(ddof=1),
                "mean_c2": estimate_c2.mean(),
                "bias_c2": (estimate_c2 - truth_c2).mean(),
                "rmse_c2": np.sqrt(((estimate_c2 - truth_c2) ** 2).mean()),
                "empirical_sd_c2": estimate_c2.std(ddof=1),
                "n_se_h2": int(valid_h2_se.sum()),
                "mean_se_h2": se_h2[valid_h2_se].mean(),
                "coverage_h2": coverage_h2,
                "n_se_c2": int(valid_c2_se.sum()),
                "mean_se_c2": se_c2[valid_c2_se].mean(),
                "coverage_c2": coverage_c2,
                "boundary_rate_h2": float((estimate_h2 <= 1e-6).mean()),
                "boundary_rate_c2": float((estimate_c2 <= 1e-6).mean()),
                "median_runtime_sec": (
                    numeric(frame["runtime_sec"]).median()
                    if "runtime_sec" in frame else np.nan
                ),
            }
        )

    summary = pd.DataFrame(summaries)
    summary.to_csv(root / "summary.csv", index=False)
    print(f"COMPLETE PredLMM rows={len(data)} groups={len(summary)}", flush=True)


if __name__ == "__main__":
    main()
