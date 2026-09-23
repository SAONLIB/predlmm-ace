#!/usr/bin/env python3
import argparse
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd


NUMBER = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"


def parse_h2(path):
    if not path.exists():
        return np.nan
    matches = re.findall(
        rf"Heritability\s*=\s*{NUMBER}",
        path.read_text(errors="replace"),
        flags=re.IGNORECASE,
    )
    return float(matches[-1]) if matches else np.nan


def run(command, log_path):
    start = time.time()
    with log_path.open("w") as log:
        return_code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    return return_code, time.time() - start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--grm-prefix", required=True)
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--cutoff", type=float, default=0.05)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    label = f"{args.cutoff:g}"
    sparse = out / f"sparse_{label}"

    rc_sparse, sparse_seconds = run(
        [
            args.gcta,
            "--grm", args.grm_prefix,
            "--make-bK-sparse", label,
            "--threads", str(args.threads),
            "--out", str(sparse),
        ],
        out / f"sparse_{label}.log",
    )

    h2 = np.nan
    rc_fastgwa = -1
    fastgwa_seconds = 0.0
    if rc_sparse == 0:
        prefix = out / f"fastgwa_{label}"
        rc_fastgwa, fastgwa_seconds = run(
            [
                args.gcta,
                "--bfile", str(data / "geno"),
                "--grm-sparse", str(sparse),
                "--fastGWA-mlm",
                "--model-only",
                "--est-vg", "REML",
                "--pheno", str(data / "phenotype.txt"),
                "--qcovar", str(data / "pcs.qcovar"),
                "--threads", str(args.threads),
                "--out", str(prefix),
            ],
            out / f"fastgwa_{label}.driver.log",
        )
        h2 = parse_h2(Path(str(prefix) + ".log"))

    result = pd.DataFrame([
        {
            "method": "fastGWA-REML",
            "sparse_cutoff": args.cutoff,
            "estimate_h2": h2,
            "runtime_sec": sparse_seconds + fastgwa_seconds,
            "success": rc_sparse == 0 and rc_fastgwa == 0,
        }
    ])
    result.to_csv(out / "fastgwa_results.csv", index=False)
    if rc_sparse != 0:
        raise SystemExit(rc_sparse)
    if rc_fastgwa != 0:
        raise SystemExit(rc_fastgwa)
    print("COMPLETE fastGWA", flush=True)


if __name__ == "__main__":
    main()
