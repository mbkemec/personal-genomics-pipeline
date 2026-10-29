
#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vmiss", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rates = []
    total_missing = 0
    total_observations = 0

    with open(args.vmiss, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        required = {"#CHROM", "ID", "MISSING_CT", "OBS_CT", "F_MISS"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected columns: {reader.fieldnames}")

        for row in reader:
            missing = int(row["MISSING_CT"])
            observed = int(row["OBS_CT"])
            rate = float(row["F_MISS"])

            if observed <= 0 or not 0 <= missing <= observed:
                raise ValueError(f"Invalid counts for {row['ID']}")

            if abs(rate - missing / observed) > 1e-5:
                raise ValueError(f"Inconsistent F_MISS for {row['ID']}")

            rates.append(rate)
            total_missing += missing
            total_observations += observed

    if not rates:
        raise ValueError("No variants found")

    thresholds = [0.01, 0.02, 0.05, 0.10]

    report = {
        "input": args.vmiss,
        "total_snps": len(rates),
        "overall_missingness": total_missing / total_observations,
        "mean_snp_missingness": mean(rates),
        "median_snp_missingness": median(rates),
        "max_snp_missingness": max(rates),
        "zero_missing_snps": sum(r == 0 for r in rates),
        "thresholds": {
            f"{threshold:.0%}": {
                "retained": sum(r <= threshold for r in rates),
                "excluded": sum(r > threshold for r in rates),
            }
            for threshold in thresholds
        },
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
        f.write("\n")

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
