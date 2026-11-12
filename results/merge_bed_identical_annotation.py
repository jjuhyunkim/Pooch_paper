#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge overlapping or adjacent BED intervals when all "
            "columns after column 3 are identical."
        )
    )
    parser.add_argument("input", help="Input BED file")
    parser.add_argument("output", help="Output BED file")
    parser.add_argument(
        "--max-gap",
        type=int,
        default=0,
        help="Maximum gap allowed between intervals [0]",
    )
    args = parser.parse_args()

    if args.max_gap < 0:
        parser.error("--max-gap must be >= 0")

    rows = []
    with open(args.input) as fin:
        for line_number, line in enumerate(fin, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"Line {line_number}: expected at least 4 columns"
                )

            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(
                    f"Line {line_number}: invalid interval"
                )

            rows.append(fields)

    # Group by chromosome and every annotation column.
    rows.sort(
        key=lambda r: (r[0], tuple(r[3:]), int(r[1]), int(r[2]))
    )

    merged = []
    for row in rows:
        if merged:
            previous = merged[-1]

            same_annotations = (
                row[0] == previous[0]
                and row[3:] == previous[3:]
            )
            close_enough = (
                int(row[1]) <= int(previous[2]) + args.max_gap
            )

            if same_annotations and close_enough:
                previous[2] = str(
                    max(int(previous[2]), int(row[2]))
                )
                continue

        merged.append(row.copy())

    merged.sort(
        key=lambda r: (r[0], int(r[1]), int(r[2]), tuple(r[3:]))
    )

    with open(args.output, "w") as fout:
        for row in merged:
            fout.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()