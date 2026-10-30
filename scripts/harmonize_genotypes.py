
#!/usr/bin/env python3

import json
import sys
from collections import Counter
from pathlib import Path

if len(sys.argv) != 6:
    sys.exit(
        "Usage: harmonize_genotypes.py "
        "<tped> <snp> <codes> <qc.json> <matched_snps.txt>"
    )

tped_file = sys.argv[1]
snp_file = sys.argv[2]
output_file = Path(sys.argv[3])
qc_file = Path(sys.argv[4])
matched_file = Path(sys.argv[5])

complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
valid_alleles = set(complement)

sample = {}
duplicate_positions = set()

with open(tped_file, encoding="utf-8") as f:
    for line in f:
        chrom, rsid, dist, pos, a1, a2 = line.split()
        key = (chrom, pos)

        if key in sample:
            duplicate_positions.add(key)
        else:
            sample[key] = (rsid, a1, a2)

counts = Counter()
counts["duplicate_sample_positions"] = len(duplicate_positions)

for path in (output_file, qc_file, matched_file):
    path.parent.mkdir(parents=True, exist_ok=True)

with (
    open(snp_file, encoding="utf-8") as snps,
    open(output_file, "w", encoding="utf-8") as out,
    open(matched_file, "w", encoding="utf-8") as matched,
):
    for line in snps:
        rsid, chrom, dist, pos, ref1, ref2 = line.split()

        key = (chrom, pos)
        code = 9

        if key in duplicate_positions:
            counts["duplicate_position_excluded"] += 1

        elif key not in sample:
            counts["not_in_sample"] += 1

        else:
            sample_rsid, a1, a2 = sample[key]

            if sample_rsid != rsid:
                counts["rsid_mismatch"] += 1

            elif "0" in (a1, a2):
                counts["missing"] += 1

            elif not {a1, a2, ref1, ref2} <= valid_alleles:
                counts["invalid_alleles"] += 1

            elif {ref1, ref2} in ({"A", "T"}, {"C", "G"}):
                counts["ambiguous"] += 1

            elif {a1, a2} <= {ref1, ref2}:
                code = int(a1 == ref1) + int(a2 == ref1)
                counts["encoded_direct"] += 1

            else:
                flipped = (
                    complement[a1],
                    complement[a2],
                )

                if set(flipped) <= {ref1, ref2}:
                    code = (
                        int(flipped[0] == ref1)
                        + int(flipped[1] == ref1)
                    )
                    counts["encoded_complement"] += 1
                else:
                    counts["allele_mismatch"] += 1

        out.write(f"{code}\n")

        if code != 9:
            matched.write(f"{rsid}\n")

        counts["total_reference_snps"] += 1

counts["encoded"] = (
    counts["encoded_direct"]
    + counts["encoded_complement"]
)

with open(qc_file, "w", encoding="utf-8") as f:
    json.dump(dict(counts), f, indent=2)
    f.write("\n")

print(dict(counts))
