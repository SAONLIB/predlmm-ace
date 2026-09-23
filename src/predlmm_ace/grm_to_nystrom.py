#!/usr/bin/env python3
"""Extract a Nystrom factor from a PLINK/GCTA triangular binary GRM.

Expected inputs are PREFIX.grm.bin and PREFIX.grm.id produced by PLINK 1.9
--make-grm-bin. The file is streamed once; no N x N array is constructed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def read_ids(path: Path) -> list[tuple[str, str]]:
    rows = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected FID IID")
            rows.append((fields[0], fields[1]))
    if len(set(rows)) != len(rows):
        raise ValueError("Duplicate (FID, IID) rows in GRM ID file.")
    return rows


def read_knot_ids(path: Path) -> list[tuple[str, str]]:
    rows = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if line_number == 1 and len(fields) >= 2 and fields[0].upper() == "FID" and fields[1].upper() == "IID":
                continue
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected FID IID")
            rows.append((fields[0], fields[1]))
    if not rows or len(set(rows)) != len(rows):
        raise ValueError("Knot file is empty or contains duplicate IDs.")
    return rows


def extract_cross_grm(grm_path: Path, n: int, knot_indices: np.ndarray, output: Path):
    expected = n * (n + 1) // 2 * np.dtype("<f4").itemsize
    actual = grm_path.stat().st_size
    if actual != expected:
        raise ValueError(
            f"Unexpected GRM size: {actual:,} bytes; expected {expected:,} for N={n}. "
            "Confirm this is a PLINK/GCTA float32 triangular .grm.bin file."
        )
    rank = knot_indices.size
    knot_position = np.full(n, -1, dtype=np.int64)
    knot_position[knot_indices] = np.arange(rank)
    is_knot = knot_position >= 0
    cross = np.lib.format.open_memmap(output, mode="w+", dtype="float32", shape=(n, rank))
    diagonal = np.empty(n, dtype=np.float64)
    prior_knots = []
    with grm_path.open("rb") as stream:
        for i in range(n):
            row = np.fromfile(stream, dtype="<f4", count=i + 1)
            if row.size != i + 1:
                raise EOFError(f"GRM ended in row {i}.")
            diagonal[i] = float(row[i])
            if prior_knots:
                idx = np.asarray(prior_knots, dtype=np.int64)
                cross[i, knot_position[idx]] = row[idx]
            if is_knot[i]:
                pos = knot_position[i]
                cross[: i + 1, pos] = row
                prior_knots.append(i)
    cross.flush()
    return cross, diagonal


def build_from_plink_grm(
    grm_prefix: str,
    knot_file: str,
    output_directory: str,
    row_chunk_size: int = 4096,
    jitter: float = 1e-8,
):
    prefix = Path(grm_prefix)
    grm_path = Path(str(prefix) + ".grm.bin")
    id_path = Path(str(prefix) + ".grm.id")
    ids = read_ids(id_path)
    knot_ids = read_knot_ids(Path(knot_file))
    lookup = {key: i for i, key in enumerate(ids)}
    missing = [key for key in knot_ids if key not in lookup]
    if missing:
        raise ValueError(f"{len(missing)} knot IDs are absent from GRM; first: {missing[:5]}")
    knot_indices = np.asarray([lookup[key] for key in knot_ids], dtype=np.int64)
    outdir = Path(output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    cross_path = outdir / "cross_grm.npy"
    cross, grm_diagonal = extract_cross_grm(grm_path, len(ids), knot_indices, cross_path)
    knot_grm = np.asarray(cross[knot_indices, :], dtype=np.float64)
    knot_grm = (knot_grm + knot_grm.T) / 2.0
    scale = max(float(np.trace(knot_grm) / max(len(knot_ids), 1)), 1.0)
    ridge = jitter * scale
    chol = np.linalg.cholesky(knot_grm + ridge * np.eye(len(knot_ids)))
    factor_path = outdir / "nystrom_factor.npy"
    factor = np.lib.format.open_memmap(
        factor_path, mode="w+", dtype="float32", shape=(len(ids), len(knot_ids))
    )
    for start in range(0, len(ids), row_chunk_size):
        sl = slice(start, min(len(ids), start + row_chunk_size))
        factor[sl] = np.linalg.solve(chol, np.asarray(cross[sl], dtype=np.float64).T).T.astype("float32")
    factor.flush()
    explained = np.empty(len(ids), dtype=np.float64)
    for start in range(0, len(ids), row_chunk_size):
        sl = slice(start, min(len(ids), start + row_chunk_size))
        block = np.asarray(factor[sl], dtype=np.float64)
        explained[sl] = np.einsum("ij,ij->i", block, block)
    raw_delta = grm_diagonal - explained
    materially_negative = raw_delta < -1e-5
    if materially_negative.any():
        raise ValueError(
            f"Diagonal correction is below -1e-5 for {materially_negative.sum()} samples; "
            "increase --jitter only after checking knot-GRM conditioning."
        )
    delta = np.maximum(raw_delta, 0.0)
    np.save(outdir / "diagonal_correction.npy", delta)
    np.save(outdir / "grm_diagonal.npy", grm_diagonal)
    np.save(outdir / "knot_grm.npy", knot_grm)
    np.save(outdir / "knot_indices.npy", knot_indices)
    with (outdir / "sample_ids.tsv").open("w") as stream:
        stream.write("FID\tIID\n")
        for fid, iid in ids:
            stream.write(f"{fid}\t{iid}\n")
    with (outdir / "knot_ids.tsv").open("w") as stream:
        stream.write("FID\tIID\n")
        for fid, iid in knot_ids:
            stream.write(f"{fid}\t{iid}\n")
    condition = float(np.linalg.cond(knot_grm + ridge * np.eye(len(knot_ids))))
    metadata = {
        "grm_prefix": str(prefix), "n_samples": len(ids), "rank": len(knot_ids),
        "jitter_added": ridge, "regularized_knot_grm_condition_number": condition,
        "diagonal_correction_min": float(delta.min()),
        "diagonal_correction_median": float(np.median(delta)),
        "diagonal_correction_max": float(delta.max()),
        "sample_order": "sample_ids.tsv; inherited exactly from .grm.id",
    }
    (outdir / "factor_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--grm-prefix", required=True)
    p.add_argument("--knots", required=True, help="FID IID file")
    p.add_argument("--out", required=True)
    p.add_argument("--row-chunk-size", type=int, default=4096)
    p.add_argument("--jitter", type=float, default=1e-8)
    args = p.parse_args()
    metadata = build_from_plink_grm(
        args.grm_prefix, args.knots, args.out, args.row_chunk_size, args.jitter
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
