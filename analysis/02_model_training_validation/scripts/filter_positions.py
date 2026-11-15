#!/usr/bin/env python3

import argparse
import bisect
import re

import pandas as pd


def parse_size(value):
    """Examples: 1000, 100kb, 100mb, 100mbp."""
    match = re.fullmatch(
        r"(\d+(?:\.\d+)?)\s*(bp|kb|kbp|mb|mbp|gb|gbp)?",
        value.strip().lower(),
    )
    if not match:
        raise argparse.ArgumentTypeError(
            "Use a size such as 1000, 100kb, or 100mb."
        )

    units = {
        None: 1, "bp": 1,
        "kb": 1_000, "kbp": 1_000,
        "mb": 1_000_000, "mbp": 1_000_000,
        "gb": 1_000_000_000, "gbp": 1_000_000_000,
    }
    return int(float(match.group(1)) * units[match.group(2)])


def main():
    parser = argparse.ArgumentParser(
        description="Remove position columns overlapping padded BED intervals."
    )
    parser.add_argument("--table", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--chrom", required=True, help="For example: chrX")
    parser.add_argument("--padding", type=parse_size, default=0)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--position-base",
        type=int,
        choices=[0, 1],
        default=0,
        help="Coordinate base of table positions [default: 0]",
    )
    args = parser.parse_args()

    # BED coordinates are always 0-based, half-open: [start, end).
    intervals = []
    matched_chrom = False

    with open(args.bed) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()
            if fields[0] != args.chrom:
                continue

            matched_chrom = True
            start, end = int(fields[1]), int(fields[2])

            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval: {line}")

            intervals.append((
                max(0, start - args.padding),
                end + args.padding,
            ))

    if not matched_chrom:
        raise ValueError(f"No BED intervals found for {args.chrom}")

    # Merge overlapping padded intervals.
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    starts = [start for start, end in merged]

    # Read whitespace-separated input; preserve metadata as strings.
    df = pd.read_csv(args.table, sep=r"\s+", dtype=str)

    keep_indices = [0, 1, 2]

    for index, column in enumerate(df.columns[3:], start=3):
        position = int(column) - args.position_base

        if position < 0:
            raise ValueError(f"Invalid position: {column}")

        interval_index = bisect.bisect_right(starts, position) - 1
        overlaps = (
            interval_index >= 0
            and position < merged[interval_index][1]
        )

        if not overlaps:
            keep_indices.append(index)

    df.iloc[:, keep_indices].to_csv(args.out, sep="\t", index=False)

    total = len(df.columns) - 3
    kept = len(keep_indices) - 3
    print(f"Position columns: {total}")
    print(f"Removed: {total - kept}")
    print(f"Retained: {kept}")


if __name__ == "__main__":
    main()