
#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Prepare an AADR .ind file for reference filtering."
    )
    parser.add_argument("--ind", required=True)
    parser.add_argument("--selected", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.selected, encoding="utf-8") as f:
        selected = [line.strip() for line in f if line.strip()]

    if len(selected) != len(set(selected)):
        raise ValueError("Duplicate IDs in selected sample list")

    selected_set = set(selected)
    found = set()
    total = 0

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.ind, encoding="utf-8") as source, \
         open(args.output, "w", encoding="utf-8") as target:

        for line_number, line in enumerate(source, start=1):
            fields = line.split()

            if len(fields) != 3:
                raise ValueError(
                    f"Invalid .ind record at line {line_number}"
                )

            sample_id, sex, population = fields

            if sample_id in selected_set:
                if sample_id in found:
                    raise ValueError(f"Duplicate reference ID: {sample_id}")
                found.add(sample_id)
            else:
                population = "Ignore"

            target.write(f"{sample_id}\t{sex}\t{population}\n")
            total += 1

    missing = selected_set - found

    if missing:
        raise ValueError(
            f"{len(missing)} selected IDs missing from reference"
        )

    print(f"Reference individuals: {total}")
    print(f"Selected individuals: {len(found)}")
    print(f"Marked as Ignore: {total - len(found)}")


if __name__ == "__main__":
    main()
