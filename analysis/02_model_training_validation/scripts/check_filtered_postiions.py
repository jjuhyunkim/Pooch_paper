#!/usr/bin/env python3

import argparse
import bisect
from collections import Counter

from filter_positions import parse_size


def read_positions(path, position_base):
    # Read only the header.
    with open(path) as handle:
        header = handle.readline().split()

    if header[:3] != ["Sample", "LeftOpHapLabel", "RightOpHapLabel"]:
        raise ValueError(f"Unexpected metadata columns in {path}")

    positions = [int(column) - position_base for column in header[3:]]

    if any(position < 0 for position in positions):
        raise ValueError(f"Negative position after coordinate conversion: {path}")

    if len(positions) != len(set(positions)):
        raise ValueError(f"Duplicate position columns found in {path}")

    return positions


def read_padded_intervals(path, chrom, padding):
    intervals = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()
            if fields[0] != chrom:
                continue

            start, end = int(fields[1]), int(fields[2])

            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval: {line}")

            intervals.append((
                max(0, start - padding),
                end + padding,
            ))

    if not intervals:
        raise ValueError(f"No BED intervals found for {chrom}")

    # Merge intervals to count each CpG only once.
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Verify CpG filtering against padded iDMR intervals."
    )
    parser.add_argument("--original", required=True)
    parser.add_argument("--filtered", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--padding", type=parse_size, default=0)
    parser.add_argument(
        "--position-base", type=int, choices=[0, 1], default=0
    )
    args = parser.parse_args()

    original = read_positions(args.original, args.position_base)
    filtered = read_positions(args.filtered, args.position_base)

    intervals = read_padded_intervals(
        args.bed, args.chrom, args.padding
    )
    starts = [start for start, end in intervals]

    def is_inside(position):
        index = bisect.bisect_right(starts, position) - 1
        return index >= 0 and position < intervals[index][1]

    original_inside = [p for p in original if is_inside(p)]
    original_outside = [p for p in original if not is_inside(p)]

    filtered_inside = [p for p in filtered if is_inside(p)]
    filtered_outside = [p for p in filtered if not is_inside(p)]

    expected = set(original_outside)
    observed = set(filtered)

    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)

    print(f"Chromosome: {args.chrom}")
    print(f"Padding on each side: {args.padding:,} bp")
    print()
    print("Dataset\tTotal_CpGs\tWithin_iDMR_plus_padding\tOutside_iDMR_plus_padding")
    print(
        f"Original\t{len(original)}\t"
        f"{len(original_inside)}\t{len(original_outside)}"
    )
    print(
        f"Filtered\t{len(filtered)}\t"
        f"{len(filtered_inside)}\t{len(filtered_outside)}"
    )

    print()
    print(f"Excluded CpGs remaining in output: {len(filtered_inside)}")
    print(f"Expected CpGs missing from output: {len(missing)}")
    print(f"Unexpected CpGs in output: {len(unexpected)}")

    # Compare every expected position, including multiplicity.
    passed = Counter(filtered) == Counter(original_outside)

    if passed:
        print("\nPASS: Output contains exactly the CpGs outside iDMR + padding.")
    else:
        # Show positions using the original table's coordinate convention.
        def preview(values):
            return [p + args.position_base for p in values[:10]]

        if filtered_inside:
            print("Excluded positions remaining (first 10):",
                  preview(filtered_inside))
        if missing:
            print("Missing positions (first 10):", preview(missing))
        if unexpected:
            print("Unexpected positions (first 10):", preview(unexpected))

        print("\nFAIL: Output does not match the expected CpG positions.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()