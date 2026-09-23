#!/usr/bin/env python3
import argparse
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd


def parse_hsq(path):
    h2 = se = np.nan
    if path.exists():
        for line in path.read_text(errors="replace").splitlines():
            fields = line.split()
            if len(fields) >= 3 and fields[0] in {"V(G)/Vp", "V(1)/Vp"}:
                h2, se = map(float, fields[1:3])
    return h2, se


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--grm-prefix", required=True)
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    keep = data / "unrelated.keep"
    if not keep.exists():
        raise FileNotFoundError(f"Missing unrelated keep file: {keep}")

    prefix = out / "gcta_unrelated"
    command = [
        args.gcta,
        "--reml",
        "--grm", args.grm_prefix,
        "--keep", str(keep),
        "--pheno", str(data / "phenotype.txt"),
        "--qcovar", str(data / "pcs.qcovar"),
        "--thread-num", str(args.threads),
        "--out", str(prefix),
    ]
    start = time.time()
    with (out / "gcta_unrelated.driver.log").open("w") as log:
        return_code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    runtime = time.time() - start
    h2, se = parse_hsq(Path(str(prefix) + ".hsq"))
    result = pd.DataFrame([
        {
            "method": "GCTA-GREML-unrelated",
            "estimate_h2": h2,
            "se_h2": se,
            "runtime_sec": runtime,
            "success": return_code == 0,
        }
    ])
    result.to_csv(out / "gcta_results.csv", index=False)
    if return_code != 0:
        raise SystemExit(return_code)
    print("COMPLETE GCTA unrelated", flush=True)


if __name__ == "__main__":
    main()
