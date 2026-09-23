#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows generated for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_specs(text):
    specs = []
    for item in text.split(","):
        method, rank = item.split(":", 1)
        method = method.strip()
        rank = int(rank)
        if method not in {"structure", "random", "family"}:
            raise ValueError(f"Unsupported knot method: {method}")
        if rank < 2:
            raise ValueError(f"Rank must be at least 2: {rank}")
        specs.append((method, rank))
    if len(specs) != len(set(specs)):
        raise ValueError("FACTOR_SPECS contains duplicates")
    return specs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--N", type=int, default=12000)
    parser.add_argument(
        "--factor-specs",
        default="structure:3000,structure:4000,structure:5000,random:5000",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.reps < 1:
        raise ValueError("--reps must be positive")

    specs = parse_specs(args.factor_specs)
    replicate_rows = []
    task_id = 0
    for simulation in (1, 2):
        c2_values = (0.10, 0.25) if simulation == 1 else (0.0,)
        for h2 in (0.25, 0.70):
            for c2 in c2_values:
                for replicate in range(1, args.reps + 1):
                    replicate_rows.append(
                        dict(
                            task_id=task_id,
                            simulation=simulation,
                            N=args.N,
                            h2=h2,
                            c2=c2,
                            replicate=replicate,
                        )
                    )
                    task_id += 1

    factor_rows = []
    factor_id = 0
    for row in replicate_rows:
        for method, rank in specs:
            if rank > args.N:
                raise ValueError(f"Rank {rank} exceeds N={args.N}")
            factor_rows.append(
                dict(
                    factor_id=factor_id,
                    replicate_task_id=row["task_id"],
                    simulation=row["simulation"],
                    N=row["N"],
                    h2=row["h2"],
                    c2=row["c2"],
                    replicate=row["replicate"],
                    knot_fraction=rank / args.N,
                    rank=rank,
                    knot_method=method,
                )
            )
            factor_id += 1

    root = Path(args.out)
    write_tsv(root / "replicates.tsv", replicate_rows)
    write_tsv(root / "factors.tsv", factor_rows)
    print(
        f"N={args.N} scenarios=6 replicates={len(replicate_rows)} "
        f"factors={len(factor_rows)} configurations={len(specs)}"
    )


if __name__ == "__main__":
    main()
