
import json
import sys
from collections import Counter
from pathlib import Path

tped_file = sys.argv[1]
snp_file = sys.argv[2]
output_file = Path(sys.argv[3])

# Read the cleaned autosomal genotypes
sample = {}

with open(tped_file) as f:
    for line in f:
        chrom, rsid, genetic_pos, pos, a1, a2 = line.split()

        if a1 == "0" or a2 == "0":
            gt = None
        else:
            gt = {a1, a2}

        sample[(chrom, pos)] = (rsid, gt)

complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

counts = Counter()

# Compare with the AADR reference panel
with open(snp_file) as f:
    for line in f:
        rsid, chrom, genetic_pos, pos, a1, a2 = line.split()

        key = (chrom, pos)

        if key not in sample:
            continue

        sample_rsid, gt = sample[key]

        counts["shared_positions"] += 1

        if sample_rsid != rsid:
            counts["rsid_mismatch"] += 1
            continue

        counts["shared_snps"] += 1

        if gt is None:
            counts["missing"] += 1
            continue

        reference = {a1, a2}

        if reference in ({"A", "T"}, {"C", "G"}):
            counts["ambiguous"] += 1
        elif gt <= reference:
            counts["matching"] += 1
        elif {complement[a] for a in gt} <= reference:
            counts["strand_flip"] += 1
        else:
            counts["mismatch"] += 1

output_file.parent.mkdir(parents=True, exist_ok=True)

with open(output_file, "w") as f:
    json.dump(dict(counts), f, indent=4)

print(dict(counts))
