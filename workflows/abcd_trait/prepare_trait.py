#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--metadata", required=True)
p.add_argument("--grm-id", required=True)
p.add_argument("--trait", required=True)
p.add_argument("--family", required=True)
p.add_argument("--continuous", required=True)
p.add_argument("--categorical", required=True)
p.add_argument("--out-dir", required=True)
a = p.parse_args()

out = Path(a.out_dir)
(out / "data").mkdir(parents=True, exist_ok=True)
meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
ids = pd.read_csv(a.grm_id, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
continuous = [x for x in a.continuous.split(",") if x]
categorical = [x for x in a.categorical.split(",") if x]
required = ["FID", "IID", a.trait, a.family] + continuous + categorical
missing_columns = [x for x in required if x not in meta.columns]
if missing_columns:
    raise SystemExit(f"Metadata is missing columns: {missing_columns}")
if meta.IID.duplicated().any():
    raise SystemExit("Metadata contains duplicate IID values")

d = ids.merge(meta, on=["FID", "IID"], how="left", validate="one_to_one", indicator=True)
if (d._merge != "both").any():
    print(f"GRM IDs without metadata: {(d._merge != 'both').sum()}")
d = d.drop(columns="_merge")
for col in [a.trait] + continuous:
    d[col] = pd.to_numeric(d[col].replace({"n/a": None, "NA": None, "": None}), errors="coerce")
for col in categorical + [a.family]:
    d[col] = d[col].replace({"n/a": None, "NA": None, "": None})

complete = d[required].notna().all(axis=1)
d = d.loc[complete].copy()
if len(d) == 0:
    raise SystemExit("No complete cases")
if d[a.family].duplicated(keep=False).sum() == 0:
    print("WARNING: complete-case sample contains no multi-participant families")

# Exact GRM order is retained. First GRM-ordered member represents each family.
one = d.drop_duplicates(a.family, keep="first").copy()
d.to_csv(out / "data" / f"{a.trait}_metadata.tsv", sep="\t", index=False)
d[["FID", "IID"]].to_csv(out / "data" / f"{a.trait}_full.keep", sep="\t", index=False, header=False)
one[["FID", "IID"]].to_csv(out / "data" / f"{a.trait}_onefamily.keep", sep="\t", index=False, header=False)
d[["FID", "IID", a.trait]].to_csv(out / "data" / f"{a.trait}.pheno", sep="\t", index=False, header=False)
d[["FID", "IID"] + continuous].to_csv(out / "data" / f"{a.trait}.qcovar", sep="\t", index=False, header=False)

# GCTA discrete covariates are emitted as deterministic integer codes.
cov = d[["FID", "IID"]].copy()
for col in categorical:
    levels = sorted(d[col].astype(str).unique())
    mapping = {v: i + 1 for i, v in enumerate(levels)}
    cov[col] = d[col].astype(str).map(mapping)
cov.to_csv(out / "data" / f"{a.trait}.covar", sep="\t", index=False, header=False)

fam_size = d.groupby(a.family).size()
qc = pd.DataFrame({"metric": [
    "full_N", "onefamily_N", "families", "multi_participant_families",
    "participants_in_multi_participant_families", "maximum_family_size"
], "value": [
    len(d), len(one), d[a.family].nunique(), int((fam_size > 1).sum()),
    int(fam_size[fam_size > 1].sum()), int(fam_size.max())
]})
qc.to_csv(out / "data" / f"{a.trait}_qc.tsv", sep="\t", index=False)
print(qc.to_string(index=False))

