
#!/usr/bin/env python3

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import yaml


def load_config(path):
    with open(path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    population_regions = {}

    for region, populations in config["regions"].items():
        for population in populations:
            if population in population_regions:
                raise ValueError(
                    f"Population assigned to multiple regions: {population}"
                )
            population_regions[population] = region

    return config, population_regions


def read_annotations(path):
    annotations = {}

    with open(path, encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        if len(header) < 48:
            raise ValueError("Unexpected AADR annotation format")

        for row in reader:
            if len(row) < 48:
                continue

            genetic_id = row[0].strip()

            if genetic_id in annotations:
                raise ValueError(f"Duplicate annotation ID: {genetic_id}")

            annotations[genetic_id] = {
                "date_bp": row[10].strip(),
                "population": row[14].strip(),
                "data_type": row[21].strip(),
                "assessment": row[47].strip(),
            }

    return annotations


def select_samples(config, population_regions, annotations, ind_path):
    settings = config["selection"]

    max_per_population = settings.get("max_individuals_per_population")

    if max_per_population is not None:
        if type(max_per_population) is not int or max_per_population < 1:
            raise ValueError("max_individuals_per_population must be positive")

    candidates = defaultdict(list)
    exclusions = defaultdict(Counter)
    observed = Counter()
    seen_ids = set()

    with open(ind_path, encoding="utf-8") as f:
        for line in f:
            fields = line.split()

            if len(fields) < 3:
                continue

            genetic_id, sex, population = fields[:3]

            if genetic_id in seen_ids:
                raise ValueError(f"Duplicate .ind ID: {genetic_id}")

            seen_ids.add(genetic_id)

            if population not in population_regions:
                continue

            observed[population] += 1
            info = annotations.get(genetic_id)

            reason = None

            if info is None:
                reason = "missing_annotation"
            elif info["population"] != population:
                reason = "population_mismatch"
            elif not genetic_id.endswith(settings["genetic_id_suffix"]):
                reason = "invalid_id_suffix"
            elif info["date_bp"] != str(settings["date_bp"]):
                reason = "invalid_date"
            elif info["data_type"] != settings["data_type"]:
                reason = "invalid_data_type"
            elif info["assessment"] != settings["assessment"]:
                reason = "failed_assessment"
            elif any(
                pattern.lower() in population.lower()
                for pattern in settings["exclude_label_patterns"]
            ):
                reason = "excluded_population_label"

            if reason:
                exclusions[population][reason] += 1
            else:
                candidates[population].append(genetic_id)

    selected = {}

    for population in population_regions:
        ids = candidates[population]

        # Preserve original reference order.
        if max_per_population is None or len(ids) <= max_per_population:
            selected[population] = ids
        else:
            raise ValueError(
                "Population downsampling is not implemented. "
                "Use max_individuals_per_population: null."
            )

    selected_ids = {
        genetic_id
        for ids in selected.values()
        for genetic_id in ids
    }

    # Preserve the original .ind order for downstream genotype extraction.
    ordered_ids = []

    with open(ind_path, encoding="utf-8") as f:
        for line in f:
            fields = line.split()
            if fields and fields[0] in selected_ids:
                ordered_ids.append(fields[0])

    report = {
        "reference_panel": config["reference_panel"],
        "selection_criteria": settings,
        "total_selected": len(ordered_ids),
        "total_candidates_in_configured_populations": sum(observed.values()),
        "total_eligible_before_downsampling": sum(
            len(ids) for ids in candidates.values()
        ),
        "missing_populations": [
            population
            for population in population_regions
            if observed[population] == 0
        ],
        "regions": {},
        "exclusions_by_reason": {},
    }

    total_exclusions = Counter()

    for region, populations in config["regions"].items():
        population_report = {}

        for population in populations:
            excluded = dict(exclusions[population])
            total_exclusions.update(excluded)

            population_report[population] = {
                "region": region,
                "observed": observed[population],
                "eligible": len(candidates[population]),
                "selected": len(selected[population]),
                "excluded": excluded,
            }

        report["regions"][region] = {
            "observed": sum(observed[p] for p in populations),
            "selected": sum(len(selected[p]) for p in populations),
            "populations": population_report,
        }

    report["exclusions_by_reason"] = dict(total_exclusions)

    assert len(ordered_ids) == len(selected_ids)
    assert (
        report["total_candidates_in_configured_populations"]
        == report["total_selected"] + sum(total_exclusions.values())
    )

    return ordered_ids, report


def main():
    parser = argparse.ArgumentParser(
        description="Select modern AADR reference individuals."
    )

    parser.add_argument("--config", required=True)
    parser.add_argument("--ind", required=True)
    parser.add_argument("--anno", required=True)
    parser.add_argument("--output-ids", required=True)
    parser.add_argument("--output-report", required=True)

    args = parser.parse_args()

    config, population_regions = load_config(args.config)
    annotations = read_annotations(args.anno)

    selected_ids, report = select_samples(
        config, population_regions, annotations, args.ind
    )

    Path(args.output_ids).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_report).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output_ids, "w", encoding="utf-8") as f:
        for genetic_id in selected_ids:
            f.write(genetic_id + "\n")

    with open(args.output_report, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(f"Selected individuals: {report['total_selected']}")
    print(f"Regions: {len(report['regions'])}")
    print(
        "Populations with selected individuals:",
        sum(
            population["selected"] > 0
            for region in report["regions"].values()
            for population in region["populations"].values()
        ),
    )
    print(f"Missing populations: {len(report['missing_populations'])}")


if __name__ == "__main__":
    main()
