
#!/usr/bin/env python3

import argparse
import json
from collections import Counter
from pathlib import Path


def read_ind(path):
    records = []

    with open(path, encoding="utf-8") as f:
        for number, line in enumerate(f, 1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{number}: invalid .ind row")
            records.append(tuple(fields))

    ids = [row[0] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate IDs in {path}")

    return records


def read_selected(path):
    with open(path, encoding="utf-8") as f:
        ids = [line.strip() for line in f if line.strip()]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate selected IDs")

    return set(ids)


def count_lines(path):
    with open(path, "rb") as f:
        return sum(1 for _ in f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-ind", required=True)
    parser.add_argument("--selected", required=True)
    parser.add_argument("--extracted-ind", required=True)
    parser.add_argument("--original-snp", required=True)
    parser.add_argument("--extracted-snp", required=True)
    parser.add_argument("--extracted-geno", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    original = read_ind(args.original_ind)
    extracted = read_ind(args.extracted_ind)
    selected = read_selected(args.selected)

    expected = [row for row in original if row[0] in selected]

    checks = {
        "selected_ids_match": {r[0] for r in extracted} == selected,
        "individual_order_and_metadata_match": extracted == expected,
        "individual_count_match": len(extracted) == len(selected),
        "snp_count_match": (
            count_lines(args.original_snp)
            == count_lines(args.extracted_snp)
        ),
        "snp_content_match": (
            Path(args.original_snp).read_bytes()
            == Path(args.extracted_snp).read_bytes()
        ),
        "genotype_file_nonempty": (
            Path(args.extracted_geno).stat().st_size > 0
        ),
    }

    report = {
        "reference_individuals": len(original),
        "selected_individuals": len(selected),
        "extracted_individuals": len(extracted),
        "original_snps": count_lines(args.original_snp),
        "extracted_snps": count_lines(args.extracted_snp),
        "extracted_genotype_bytes": (
            Path(args.extracted_geno).stat().st_size
        ),
        "population_counts": dict(sorted(Counter(
            row[2] for row in extracted
        ).items())),
        "checks": checks,
        "status": "PASS" if all(checks.values()) else "FAIL",
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(json.dumps({
        "status": report["status"],
        "individuals": len(extracted),
        "snps": report["extracted_snps"],
        "checks": checks,
    }, indent=2))

    if report["status"] != "PASS":
        raise SystemExit("Reference validation failed")


if __name__ == "__main__":
    main()
