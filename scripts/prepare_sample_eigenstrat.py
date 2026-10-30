#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp", required=True)
    parser.add_argument("--codes", required=True)
    parser.add_argument("--selected", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sex", choices=["M", "F", "U"], required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    selected = set(Path(args.selected).read_text().splitlines())

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    written = set()

    with (
        open(args.snp, encoding="utf-8") as snps,
        open(args.codes, encoding="utf-8") as codes,
        open(f"{prefix}.snp", "w", encoding="utf-8") as snp_out,
        open(f"{prefix}.geno", "w", encoding="utf-8") as geno_out,
    ):
        for snp_line, code_line in zip(snps, codes, strict=True):
            snp_id = snp_line.split()[0]

            if snp_id not in selected:
                continue

            code = code_line.strip()

            if code not in {"0", "1", "2", "9"}:
                raise ValueError(f"Invalid genotype code: {snp_id}={code}")

            if snp_id in written:
                raise ValueError(f"Duplicate SNP: {snp_id}")

            written.add(snp_id)
            snp_out.write(snp_line)
            geno_out.write(code + "\n")

    if written != selected:
        missing = selected - written
        raise ValueError(
            f"{len(missing)} selected SNPs not found in reference. "
            f"Examples: {sorted(missing)[:5]}"
        )

    with open(f"{prefix}.ind", "w", encoding="utf-8") as ind_out:
        ind_out.write(f"{args.sample_id} {args.sex} TARGET\n")

    print(f"Prepared {len(written)} SNPs for {args.sample_id}")


if __name__ == "__main__":
    main()
