#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Sum BED column 4 within non-overlapping genomic bins."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-o", "--output", required=True, help="Output bedGraph file")
    parser.add_argument(
        "--bin-size", type=int, default=50, help="Bin size in bp [50]"
    )
    args = parser.parse_args()

    if args.bin_size <= 0:
        parser.error("--bin-size must be greater than 0")

    df = pd.read_csv(
        args.input,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "value"],
        comment="#",
    )

    df["bin_start"] = (df["start"] // args.bin_size) * args.bin_size

    result = (
        df.groupby(["chrom", "bin_start"], as_index=False, sort=True)["value"]
        .sum()
        .rename(columns={"bin_start": "start"})
    )

    result["end"] = result["start"] + args.bin_size
    result[["chrom", "start", "end", "value"]].to_csv(
        args.output, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    main()