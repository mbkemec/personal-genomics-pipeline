
#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_evec(path):
    rows = []

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            fields = line.split()
            coordinates = [float(value) for value in fields[1:-1]]

            row = {
                "sample_id": fields[0],
                "population": fields[-1],
            }

            for i, value in enumerate(coordinates, start=1):
                row[f"PC{i}"] = value

            rows.append(row)

    df = pd.DataFrame(rows)

    if df.empty:
        raise ValueError("No PCA samples found")

    if df["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample IDs in PCA output")

    return df


def read_metadata(path):
    metadata = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # AADR v66.p1 compatibility_HO annotation columns.
    columns = metadata.columns

    selected = metadata.iloc[:, [0, 6, 7, 15, 16]].copy()

    selected.columns = [
        "sample_id",
        "publication",
        "doi",
        "locality",
        "country",
    ]

    selected["sample_id"] = selected["sample_id"].str.strip()

    if selected["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample IDs in AADR metadata")

    return selected


def calculate_distances(reference, target, n_components):
    columns = [
        f"PC{i}" for i in range(1, n_components + 1)
    ]

    if not set(columns).issubset(reference.columns):
        raise ValueError(
            f"Missing components for PC1-PC{n_components}"
        )

    differences = (
        reference[columns].to_numpy()
        - target[columns].to_numpy()[0]
    )

    distances = np.linalg.norm(differences, axis=1)

    result = reference[
        ["sample_id", "population"]
    ].copy()

    result["distance"] = distances

    summary = (
        result.groupby("population")
        .agg(
            n=("sample_id", "size"),
            median_distance=("distance", "median"),
            mean_distance=("distance", "mean"),
        )
        .reset_index()
        .sort_values(
            ["median_distance", "population"]
        )
        .reset_index(drop=True)
    )

    return summary


def calculate_individual_distances(reference, target, n_components):
    columns = [
        f"PC{i}" for i in range(1, n_components + 1)
    ]

    if not set(columns).issubset(reference.columns):
        raise ValueError(
            f"Missing components for PC1-PC{n_components}"
        )

    differences = (
        reference[columns].to_numpy()
        - target[columns].to_numpy()[0]
    )

    distances = np.linalg.norm(differences, axis=1)

    result = reference[
        ["sample_id", "population"]
    ].copy()

    result["distance"] = distances

    result = (
        result.sort_values(
            ["distance", "sample_id"]
        )
        .reset_index(drop=True)
    )

    result.insert(
        0,
        "rank",
        np.arange(1, len(result) + 1)
    )

    return result


def plot_scatter(reference, target, selected, output, zoom=False):
    fig, ax = plt.subplots(figsize=(13, 9))

    colors = plt.get_cmap("tab20", len(selected))

    if not zoom:
        background = reference[
            ~reference["population"].isin(selected)
        ]

        ax.scatter(
            background["PC1"],
            background["PC2"],
            s=14,
            alpha=0.2,
            color="gray",
            label="Other reference populations",
            rasterized=True,
        )

    for i, population in enumerate(selected):
        subset = reference[
            reference["population"] == population
        ]

        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            s=35 if zoom else 25,
            alpha=0.75,
            color=colors(i),
            label=f"{population} (n={len(subset)})",
            rasterized=True,
        )

    target_x = target["PC1"].iloc[0]
    target_y = target["PC2"].iloc[0]

    ax.scatter(
        target_x,
        target_y,
        s=300 if zoom else 260,
        marker="*",
        color="red",
        edgecolors="black",
        linewidths=0.8,
        label=target["sample_id"].iloc[0],
        zorder=10,
    )

    if zoom:
        ax.set_xlim(target_x - 0.018, target_x + 0.018)
        ax.set_ylim(target_y - 0.018, target_y + 0.018)
        ax.set_title("AADR PCA - Target Neighborhood")
    else:
        ax.set_title(
            "AADR PCA - Nearest Reference Populations"
        )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.grid(alpha=0.15)

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_distance_ranking(summary, n_components, output):
    # Only sufficiently represented populations are plotted.
    plotted = (
        summary[summary["n"] >= 5]
        .head(15)
        .iloc[::-1]
    )

    if plotted.empty:
        raise ValueError(
            "No populations meet the minimum sample size"
        )

    labels = [
        f"{row.population} (n={row.n})"
        for row in plotted.itertuples()
    ]

    fig, ax = plt.subplots(figsize=(12, 8))

    bars = ax.barh(
        labels,
        plotted["median_distance"],
        color="steelblue",
        alpha=0.85,
    )

    ax.bar_label(
        bars,
        fmt="%.4f",
        padding=4,
        fontsize=9,
    )

    ax.set_xlim(
        0,
        plotted["median_distance"].max() * 1.22,
    )

    ax.set_xlabel(
        f"Median Euclidean distance (PC1-PC{n_components})"
    )

    ax.set_title(
        f"Nearest Reference Populations - First {n_components} PCs"
    )

    ax.grid(axis="x", alpha=0.2)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evec", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = read_evec(args.evec)

    target = df[df["sample_id"] == args.target]
    reference = df[df["sample_id"] != args.target].copy()

    if len(target) != 1:
        raise ValueError(
            "Expected exactly one target sample"
        )

    metadata = read_metadata(args.metadata)
    summaries = {}

    for n_components in (2, 5, 10):
        summary = calculate_distances(
            reference,
            target,
            n_components,
        )

        summaries[n_components] = summary

    # Existing PC1-PC2 outputs.
    summaries[2].to_csv(
        output_dir / "pca_population_distances.csv",
        index=False,
    )

    selected = (
        summaries[2]
        .head(12)["population"]
        .tolist()
    )

    plot_scatter(
        reference,
        target,
        selected,
        output_dir / "pca_colored.png",
    )

    plot_scatter(
        reference,
        target,
        selected,
        output_dir / "pca_zoom.png",
        zoom=True,
    )

    # Additional PC1-PC5 and PC1-PC10 outputs.
    for n_components in (5, 10):
        prefix = f"pca_pc1_pc{n_components}"

        summaries[n_components].to_csv(
            output_dir / f"{prefix}.csv",
            index=False,
        )

        plot_distance_ranking(
            summaries[n_components],
            n_components,
            output_dir / f"{prefix}.png",
        )


    # Individual-level PCA distances with AADR metadata.
    for n_components in (2, 5, 10):
        individuals = calculate_individual_distances(
            reference,
            target,
            n_components,
        )

        individuals = individuals.merge(
            metadata,
            on="sample_id",
            how="left",
            validate="one_to_one",
            indicator=True,
        )

        unmatched = individuals["_merge"] != "both"

        if unmatched.any():
            missing_ids = individuals.loc[
                unmatched, "sample_id"
            ].tolist()

            raise ValueError(
                f"Missing AADR metadata for {len(missing_ids)} "
                f"reference samples. Examples: {missing_ids[:5]}"
            )

        individuals = individuals.drop(columns="_merge")

        individuals = individuals[
            [
                "rank",
                "sample_id",
                "population",
                "locality",
                "country",
                "publication",
                "doi",
                "distance",
            ]
        ]

        individuals.to_csv(
            output_dir
            / f"pca_nearest_individuals_pc{n_components}.csv",
            index=False,
        )


    print("PCA visualization completed.")
    print("Generated 10 output files.")


if __name__ == "__main__":
    main()
