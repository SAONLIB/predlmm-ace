#!/usr/bin/env python3
"""Select family- and stratum-aware knots aligned to a PLINK GRM."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def balanced_take(frame, n, strata, rng):
    if n <= 0 or frame.empty:
        return frame.iloc[0:0]
    if not strata:
        return frame.iloc[rng.permutation(len(frame))[:n]]
    work = frame.copy()
    work["_stratum"] = work[strata].astype(str).agg("|".join, axis=1)
    groups = {k: g.iloc[rng.permutation(len(g))].copy() for k, g in work.groupby("_stratum", sort=True)}
    selected = []
    while len(selected) < n and groups:
        for key in list(groups):
            if len(selected) >= n:
                break
            group = groups[key]
            if group.empty:
                groups.pop(key)
            else:
                selected.append(group.iloc[[0]])
                groups[key] = group.iloc[1:]
    return pd.concat(selected, ignore_index=False).drop(columns="_stratum") if selected else work.iloc[0:0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True, help="TSV containing FID, IID, family and strata")
    p.add_argument("--grm-id", required=True)
    p.add_argument("--family-column", required=True)
    p.add_argument("--strata", default="", help="Comma-separated columns, e.g. ancestry_group,site")
    p.add_argument("--rank", type=int, required=True)
    p.add_argument("--seed", type=int, default=20260918)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    required = {"FID", "IID", args.family_column}
    strata = [x for x in args.strata.split(",") if x]
    required.update(strata)
    missing_columns = required - set(meta)
    if missing_columns:
        raise ValueError(f"Missing metadata columns: {sorted(missing_columns)}")
    if meta[["FID", "IID"]].duplicated().any():
        raise ValueError("Duplicate FID/IID in metadata.")
    ids = pd.read_csv(args.grm_id, sep=r"\s+", names=["FID", "IID"], dtype=str)
    ids["_grm_order"] = np.arange(len(ids))
    data = ids.merge(meta, on=["FID", "IID"], how="left", validate="one_to_one", indicator=True)
    if (data["_merge"] != "both").any():
        raise ValueError(f"Metadata absent for {(data['_merge'] != 'both').sum()} GRM samples.")
    data = data.drop(columns="_merge")
    if data[list(required)].isna().any().any():
        raise ValueError("Missing family or stratum values are not allowed.")
    if not 0 < args.rank <= len(data):
        raise ValueError("rank must be between 1 and N.")

    family_size = data.groupby(args.family_column)["IID"].transform("size")
    related = data[family_size > 1]
    # One randomized representative per multi-participant family, then balance those representatives.
    reps = related.groupby(args.family_column, sort=False, group_keys=False).sample(n=1, random_state=args.seed)
    first = balanced_take(reps, min(args.rank, len(reps)), strata, rng)
    selected_orders = set(first["_grm_order"])
    remaining_n = args.rank - len(first)
    remaining = data[~data["_grm_order"].isin(selected_orders)]
    second = balanced_take(remaining, remaining_n, strata, rng)
    selected = pd.concat([first, second]).sort_values("_grm_order")
    selected[["FID", "IID"]].to_csv(args.out, sep="\t", index=False)
    represented = selected[args.family_column].nunique()
    represented_multi = selected[selected[args.family_column].isin(related[args.family_column])][args.family_column].nunique()
    print(f"N={len(data)} rank={len(selected)} represented_families={represented}")
    print(f"multi_participant_families={related[args.family_column].nunique()} represented={represented_multi}")
    print(f"seed={args.seed} output={args.out}")


if __name__ == "__main__":
    main()
