
#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def read_ids(path):
    with open(path, encoding="utf-8") as f:
        ids = [line.strip() for line in f if line.strip()]

    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate SNP IDs in {path}")

    return ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matched", required=True)
    parser.add_argument("--bim", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    matched = read_ids(args.matched)
    matched_set = set(matched)

    shared = []
    reference_ids = set()
    reference_count = 0

    with open(args.bim, encoding="utf-8") as f:
        for number, line in enumerate(f, 1):
            fields = line.split()

            if len(fields) != 6:
                raise ValueError(f"Invalid BIM row: {number}")

            snp_id = fields[1]

            if snp_id in reference_ids:
                raise ValueError(f"Duplicate reference SNP: {snp_id}")

            reference_ids.add(snp_id)
            reference_count += 1

            if snp_id in matched_set:
                shared.append(snp_id)

    report = {
        "matched_sample_snps": len(matched),
        "qc_reference_snps": reference_count,
        "shared_snps": len(shared),
        "matched_snps_excluded_by_reference_qc": (
            len(matched) - len(shared)
        ),
        "shared_fraction_of_matched": (
            len(shared) / len(matched)
        ),
        "selection": {
            "sample": "successfully harmonized HO SNPs",
            "reference": "HO SNPs passing 5% missingness filter",
            "ld_pruning_applied": False,
        },
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as f:
        for snp_id in shared:
            f.write(f"{snp_id}\n")

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
        f.write("\n")

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
