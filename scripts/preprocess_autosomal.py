import sys
import csv
import json
from pathlib import Path

input_file = sys.argv[1]
output_dir = Path(sys.argv[2])
sample_id = sys.argv[3]
sex = sys.argv[4]
output_dir.mkdir(parents=True, exist_ok=True)

variants = {}
conflict_positions = set()

total = 0
excluded = 0
duplicates = 0
conflicts = 0

with open(input_file) as f:
    reader = csv.DictReader(f)

    for row in reader:
        total += 1

        chrom = row["CHROMOSOME"]
        pos = row["POSITION"]
        gt = row["RESULT"]

        if chrom not in [str(i) for i in range(1, 23)]:
            excluded += 1
            continue

        if not pos.isdigit() or int(pos) <= 0:
            excluded += 1
            continue

        if gt != "--" and (len(gt) != 2 or any(b not in "ACGT" for b in gt)):
            excluded += 1
            continue

        key = (chrom, int(pos))

        if key in conflict_positions:
            continue

        if key not in variants:
            variants[key] = {
                "rsid": row["RSID"],
                "gt": gt
            }
        else:
            duplicates += 1
            old_gt = variants[key]["gt"]

            if old_gt == "--" and gt != "--":
                variants[key] = {
                    "rsid": row["RSID"],
                    "gt": gt
                }

            elif old_gt != "--" and gt != "--":
                if sorted(old_gt) != sorted(gt):
                    conflicts += 1
                    conflict_positions.add(key)
                    del variants[key]

tped_file = output_dir / "autosomal.tped"
tfam_file = output_dir / "autosomal.tfam"

with open(tped_file, "w") as f:
    for (chrom, pos), variant in sorted(
        variants.items(), key=lambda x: (int(x[0][0]), x[0][1])
    ):
        gt = variant["gt"]

        if gt == "--":
            a1, a2 = "0", "0"
        else:
            a1, a2 = gt[0], gt[1]

        rsid = variant["rsid"]

        f.write(f"{chrom}\t{rsid}\t0\t{pos}\t{a1}\t{a2}\n")

with open(tfam_file, "w") as f:
    f.write(f"FTDNA {sample_id} 0 0 {sex} -9\n")

summary = {
    "total_records": total,
    "excluded_records": excluded,
    "duplicate_records": duplicates,
    "conflicting_positions": conflicts,
    "output_variants": len(variants)
}

with open(output_dir / "autosomal_qc.json", "w") as f:
    json.dump(summary, f, indent=4)

print(summary)
