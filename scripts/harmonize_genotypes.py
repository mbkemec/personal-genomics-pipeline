
import json
import sys
from collections import Counter
from pathlib import Path

tped_file = sys.argv[1]
snp_file = sys.argv[2]
output_file = Path(sys.argv[3])
qc_file = Path(sys.argv[4])

sample = {}

with open(tped_file) as f:
    for line in f:
        chrom, rsid, dist, pos, a1, a2 = line.split()
        sample[(chrom, pos)] = (rsid, a1, a2)

counts = Counter()
complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

output_file.parent.mkdir(parents=True, exist_ok=True)
qc_file.parent.mkdir(parents=True, exist_ok=True)

with open(snp_file) as snps, open(output_file, "w") as out:
    for line in snps:
        rsid, chrom, dist, pos, ref1, ref2 = line.split()

        key = (chrom, pos)
        code = 9

        if key not in sample:
            counts["not_in_sample"] += 1
        else:
            sample_rsid, a1, a2 = sample[key]

            if sample_rsid != rsid:
                counts["rsid_mismatch"] += 1
            elif "0" in (a1, a2):
                counts["missing"] += 1
            elif {ref1, ref2} in ({"A", "T"}, {"C", "G"}):
                counts["ambiguous"] += 1
            elif {a1, a2} <= {ref1, ref2}:
                code = int(a1 == ref1) + int(a2 == ref1)
                counts["encoded"] += 1
            else:
                counts["allele_mismatch"] += 1

        out.write(str(code) + "\n")

counts["total_reference_snps"] = sum(counts.values())

with open(qc_file, "w") as f:
    json.dump(dict(counts), f, indent=4)

print(dict(counts))
