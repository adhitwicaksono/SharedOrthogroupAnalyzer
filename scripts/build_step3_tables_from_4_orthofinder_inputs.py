#!/usr/bin/env python3
"""
Generate Step 3.3 orthogroup tables and a preliminary Step 3.4 QC summary
from four corrected OrthoFinder PEP outputs.

Required input files:
1. OrthoFinder_PEP_2026-06-23_corrected_Orthogroups.tsv
2. OrthoFinder_PEP_2026-06-23_corrected_UnassignedGenes.tsv
3. OrthoFinder_PEP_2026-06-23_corrected_Statistics_PerSpecies.tsv
4. OrthoFinder_PEP_2026-06-23_corrected_Statistics_Overall.tsv

Outputs:
- orthogroup_presence_absence.tsv
- orthogroup_gene_counts.tsv
- orthogroup_classification.tsv
- core_orthogroups.tsv
- softcore_orthogroups.tsv
- accessory_orthogroups.tsv
- private_orthogroups.tsv
- unassigned_singleton_genes.tsv
- private_orthogroups_by_cultivar.tsv
- private_orthogroup_summary_by_cultivar.tsv
- quality_control_summary.tsv
- orthogroup_class_summary.tsv
- cultivar_pair_column_similarity.tsv
- reference_only_orthogroups.tsv
- README_processing_report.txt
"""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------

CULTIVARS = ["BSM", "CI", "GJ", "IN", "KM", "MPE", "PP"]
REFERENCE = "AGIS1.0"
ALL_DATASETS = [REFERENCE] + CULTIVARS

ORTHOGROUPS_FILE = Path(
    "OrthoFinder_PEP_2026-07-20_corrected_Orthogroups.tsv"
)
UNASSIGNED_FILE = Path(
    "OrthoFinder_PEP_2026-07-20_corrected_UnassignedGenes.tsv"
)
PER_SPECIES_FILE = Path(
    "OrthoFinder_PEP_2026-07-20_corrected_Statistics_PerSpecies.tsv"
)
OVERALL_FILE = Path(
    "OrthoFinder_PEP_2026-07-20_corrected_Statistics_Overall.tsv"
)

OUTPUT_DIR = Path("OrthoFinder_PEP_2026-07-20_Step3_tables")
OUTPUT_ZIP = Path("OrthoFinder_PEP_2026-07-20_Step3_tables.zip")


# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

def clean_species_name(name: str) -> str:
    """Convert Galaxy-style OrthoFinder labels to short dataset names."""
    if name == "Orthogroup":
        return name

    match = re.search(r"Augustus on (.+?)_ Protein sequence", str(name))
    if match:
        return match.group(1).strip()

    return str(name).strip()


def split_gene_ids(value) -> list[str]:
    """Split a comma-separated OrthoFinder cell into individual gene IDs."""
    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text:
        return []

    return [item.strip() for item in text.split(",") if item.strip()]


def classify_occupancy(n: int) -> str:
    """Classify orthogroups using occupancy among the all cultivars."""
    if n == 7:
        return "Core"
    if n == 6:
        return "Soft core"
    if 2 <= n <= 5:
        return "Accessory"
    if n == 1:
        return "Private"
    if n == 0:
        return "Reference-only"
    return "Unclassified"


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing input file(s):\n- " + "\n- ".join(missing)
        )


def parse_orthofinder_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    ).rename(columns=clean_species_name)

    required = ["Orthogroup"] + ALL_DATASETS
    missing = [column for column in required if column not in df.columns]

    if missing:
        raise ValueError(
            f"{path.name} is missing required columns: {missing}"
        )

    return df[required].copy()


def parse_per_species_statistics(path: Path) -> pd.DataFrame:
    """
    Parse the initial summary block of Statistics_PerSpecies.tsv.
    """
    lines = path.read_text(encoding="utf-8-sig").splitlines()

    if len(lines) < 10:
        raise ValueError(
            f"{path.name} is unexpectedly short or malformed."
        )

    header = lines[0].split("\t")
    species_names = [clean_species_name(x) for x in header[1:]]

    metric_map = {
        "Number of genes": "total_predicted_proteins",
        "Number of genes in orthogroups": "genes_assigned_to_orthogroups",
        "Number of unassigned genes": "unassigned_genes",
        "Percentage of genes in orthogroups": "percentage_assigned",
        "Percentage of unassigned genes": "percentage_unassigned",
        "Number of orthogroups containing species":
            "orthogroups_containing_cultivar",
        "Percentage of orthogroups containing species":
            "percentage_orthogroups_containing_cultivar",
        "Number of species-specific orthogroups":
            "orthofinder_species_specific_orthogroups_nine_dataset_run",
        "Number of genes in species-specific orthogroups":
            "genes_in_orthofinder_species_specific_orthogroups",
    }

    records = {
        species: {"cultivar": species}
        for species in species_names
    }

    for line in lines[1:12]:
        fields = line.split("\t")
        metric_name = fields[0]

        if metric_name not in metric_map:
            continue

        output_name = metric_map[metric_name]

        for species, value in zip(species_names, fields[1:]):
            records[species][output_name] = value

    stats = pd.DataFrame(records.values())

    for column in stats.columns:
        if column != "cultivar":
            stats[column] = pd.to_numeric(
                stats[column],
                errors="coerce",
            )

    return stats


def parse_overall_statistics(path: Path) -> dict[str, str]:
    metrics = {}

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            metrics[parts[0]] = parts[1]

    return metrics


# ---------------------------------------------------------------------
# MAIN PIPELINE
# ---------------------------------------------------------------------

def main() -> None:
    require_files([
        ORTHOGROUPS_FILE,
        UNASSIGNED_FILE,
        PER_SPECIES_FILE,
        OVERALL_FILE,
    ])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    membership = parse_orthofinder_matrix(ORTHOGROUPS_FILE)
    unassigned = parse_orthofinder_matrix(UNASSIGNED_FILE)

    # -------------------------------------------------------------
    # DUPLICATE-CULTIVAR GUARD
    # -------------------------------------------------------------

    pair_checks = []

    for i, cultivar_a in enumerate(CULTIVARS):
        for cultivar_b in CULTIVARS[i + 1:]:
            identical_membership = int(
                (
                    membership[cultivar_a]
                    == membership[cultivar_b]
                ).sum()
            )

            counts_a = membership[cultivar_a].map(
                lambda value: len(split_gene_ids(value))
            )
            counts_b = membership[cultivar_b].map(
                lambda value: len(split_gene_ids(value))
            )

            identical_counts = int((counts_a == counts_b).sum())

            membership_fraction = (
                identical_membership / len(membership)
            )
            count_fraction = identical_counts / len(membership)

            pair_checks.append({
                "cultivar_a": cultivar_a,
                "cultivar_b": cultivar_b,
                "identical_membership_cells": identical_membership,
                "total_orthogroups": len(membership),
                "identical_membership_fraction": membership_fraction,
                "identical_gene_count_cells": identical_counts,
                "identical_gene_count_fraction": count_fraction,
                "critical_duplicate_flag": (
                    membership_fraction >= 0.999
                    or count_fraction >= 0.9999
                ),
            })

    pair_checks_df = pd.DataFrame(pair_checks)

    pair_checks_df.to_csv(
        OUTPUT_DIR / "cultivar_pair_column_similarity.tsv",
        sep="\t",
        index=False,
    )

    duplicate_pairs = pair_checks_df[
        pair_checks_df["critical_duplicate_flag"]
    ]

    if not duplicate_pairs.empty:
        pairs = ", ".join(
            duplicate_pairs["cultivar_a"]
            + "-"
            + duplicate_pairs["cultivar_b"]
        )

        raise RuntimeError(
            "Critical duplicated or near-duplicated cultivar "
            f"columns detected: {pairs}. Pipeline aborted."
        )

    # -------------------------------------------------------------
    # GENE-COUNT MATRIX
    # -------------------------------------------------------------

    gene_counts = membership.copy()

    for dataset in ALL_DATASETS:
        gene_counts[dataset] = gene_counts[dataset].map(
            lambda value: len(split_gene_ids(value))
        )

    # -------------------------------------------------------------
    # PRESENCE/ABSENCE MATRIX
    # -------------------------------------------------------------

    presence_absence = gene_counts.copy()

    for dataset in ALL_DATASETS:
        presence_absence[dataset] = (
            presence_absence[dataset] > 0
        ).astype(int)

    cultivar_occupancy = presence_absence[
        CULTIVARS
    ].sum(axis=1)

    classes = cultivar_occupancy.map(classify_occupancy)

    private_cultivar = []

    for _, row in presence_absence[CULTIVARS].iterrows():
        present = [
            cultivar
            for cultivar in CULTIVARS
            if row[cultivar] == 1
        ]

        private_cultivar.append(
            present[0] if len(present) == 1 else ""
        )

    presence_output = presence_absence.copy()
    presence_output["cultivar_occupancy"] = cultivar_occupancy
    presence_output["orthogroup_class"] = classes
    presence_output["private_cultivar"] = private_cultivar

    presence_output.to_csv(
        OUTPUT_DIR / "orthogroup_presence_absence.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # GENE-COUNT OUTPUT
    # -------------------------------------------------------------

    gene_count_output = gene_counts.copy()

    gene_count_output["cultivar_occupancy"] = (
        cultivar_occupancy
    )
    gene_count_output["total_genes_all_cultivars"] = (
        gene_counts[CULTIVARS].sum(axis=1)
    )
    gene_count_output[
        "total_genes_including_reference"
    ] = gene_counts[ALL_DATASETS].sum(axis=1)
    gene_count_output["orthogroup_class"] = classes
    gene_count_output["private_cultivar"] = private_cultivar

    gene_count_output.to_csv(
        OUTPUT_DIR / "orthogroup_gene_counts.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # CLASSIFICATION OUTPUT
    # -------------------------------------------------------------

    classification = pd.DataFrame({
        "Orthogroup": membership["Orthogroup"],
        "cultivar_occupancy": cultivar_occupancy,
        "orthogroup_class": classes,
        "private_cultivar": private_cultivar,
        "reference_present":
            presence_absence[REFERENCE].astype(int),
        "total_genes_all_cultivars":
            gene_counts[CULTIVARS].sum(axis=1),
        "total_genes_including_reference":
            gene_counts[ALL_DATASETS].sum(axis=1),
    })

    classification.to_csv(
        OUTPUT_DIR / "orthogroup_classification.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # CATEGORY TABLES
    # -------------------------------------------------------------

    category_files = {
        "Core": "core_orthogroups.tsv",
        "Soft core": "softcore_orthogroups.tsv",
        "Accessory": "accessory_orthogroups.tsv",
        "Private": "private_orthogroups.tsv",
    }

    for category, filename in category_files.items():
        category_ids = classification.loc[
            classification["orthogroup_class"] == category,
            [
                "Orthogroup",
                "cultivar_occupancy",
                "private_cultivar",
                "reference_present",
                "total_genes_all_cultivars",
            ],
        ]

        category_membership = membership.merge(
            category_ids,
            on="Orthogroup",
            how="inner",
        )

        category_membership.to_csv(
            OUTPUT_DIR / filename,
            sep="\t",
            index=False,
        )

    # Reference-only orthogroups
    reference_only_ids = classification.loc[
        classification["orthogroup_class"]
        == "Reference-only",
        [
            "Orthogroup",
            "cultivar_occupancy",
            "reference_present",
            "total_genes_including_reference",
        ],
    ]

    membership.merge(
        reference_only_ids,
        on="Orthogroup",
        how="inner",
    ).to_csv(
        OUTPUT_DIR / "reference_only_orthogroups.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # PRIVATE ORTHOGROUPS BY CULTIVAR
    # -------------------------------------------------------------

    private_rows = []

    private_classification = classification[
        classification["orthogroup_class"] == "Private"
    ]

    for row in private_classification.itertuples(index=False):
        cultivar = row.private_cultivar

        gene_text = membership.loc[
            membership["Orthogroup"] == row.Orthogroup,
            cultivar,
        ].iloc[0]

        gene_ids = split_gene_ids(gene_text)

        private_rows.append({
            "cultivar": cultivar,
            "orthogroup": row.Orthogroup,
            "gene_count": len(gene_ids),
            "gene_ids": ", ".join(gene_ids),
            "reference_present": row.reference_present,
        })

    private_long = pd.DataFrame(private_rows)

    if not private_long.empty:
        private_long = private_long.sort_values(
            ["cultivar", "orthogroup"]
        )

    private_long.to_csv(
        OUTPUT_DIR / "private_orthogroups_by_cultivar.tsv",
        sep="\t",
        index=False,
    )

    private_summary = (
        private_long.groupby(
            "cultivar",
            as_index=False,
        )
        .agg(
            private_orthogroups=(
                "orthogroup",
                "nunique",
            ),
            genes_in_private_orthogroups=(
                "gene_count",
                "sum",
            ),
        )
        .set_index("cultivar")
        .reindex(CULTIVARS)
        .reset_index()
    )

    private_summary.to_csv(
        OUTPUT_DIR
        / "private_orthogroup_summary_by_cultivar.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # UNASSIGNED SINGLETON GENES
    # -------------------------------------------------------------

    unassigned_rows = []

    for row in unassigned.itertuples(index=False):
        row_dict = row._asdict()

        for cultivar in CULTIVARS:
            for gene_id in split_gene_ids(
                row_dict[cultivar]
            ):
                unassigned_rows.append({
                    "cultivar": cultivar,
                    "orthogroup_label":
                        row_dict["Orthogroup"],
                    "gene_id": gene_id,
                })

    unassigned_long = pd.DataFrame(unassigned_rows)

    if not unassigned_long.empty:
        unassigned_long = unassigned_long.sort_values(
            [
                "cultivar",
                "orthogroup_label",
                "gene_id",
            ]
        )

    unassigned_long.to_csv(
        OUTPUT_DIR / "unassigned_singleton_genes.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # ORTHOFINDER PER-SPECIES QC
    # -------------------------------------------------------------

    species_stats = parse_per_species_statistics(
        PER_SPECIES_FILE
    )

    species_stats = species_stats[
        species_stats["cultivar"].isin(CULTIVARS)
    ].copy()

    unassigned_summary = (
        unassigned_long.groupby(
            "cultivar",
            as_index=False,
        )
        .agg(
            unassigned_singleton_genes_from_matrix=(
                "gene_id",
                "count",
            )
        )
        .set_index("cultivar")
        .reindex(CULTIVARS)
        .fillna(0)
        .reset_index()
    )

    qc = (
        species_stats
        .merge(
            private_summary,
            on="cultivar",
            how="left",
        )
        .merge(
            unassigned_summary,
            on="cultivar",
            how="left",
        )
    )

    # Blank columns to be filled later using BUSCO and assembly reports
    pending_columns = [
        "busco_complete_pct",
        "busco_single_copy_pct",
        "busco_duplicated_pct",
        "busco_fragmented_pct",
        "busco_missing_pct",
        "assembly_size_bp",
        "contig_count",
        "contig_n50_bp",
        "contig_l50",
        "longest_contig_bp",
        "gc_pct",
    ]

    for column in pending_columns:
        qc[column] = pd.NA

    # Preliminary IQR-based QC flags
    available_metrics = [
        "total_predicted_proteins",
        "genes_assigned_to_orthogroups",
        "unassigned_genes",
        "private_orthogroups",
        "genes_in_private_orthogroups",
        "orthogroups_containing_cultivar",
    ]

    for metric in available_metrics:
        values = pd.to_numeric(
            qc[metric],
            errors="coerce",
        )

        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr

        qc[f"{metric}_iqr_outlier"] = (
            (values < lower)
            | (values > upper)
        )

    def qc_comment(row: pd.Series) -> str:
        comments = []

        if row[
            "total_predicted_proteins_iqr_outlier"
        ]:
            comments.append(
                "unusual predicted-protein count"
            )

        if row["unassigned_genes_iqr_outlier"]:
            comments.append(
                "unusual unassigned-gene count"
            )

        if row["private_orthogroups_iqr_outlier"]:
            comments.append(
                "unusual private-orthogroup count"
            )

        if row[
            "genes_in_private_orthogroups_iqr_outlier"
        ]:
            comments.append(
                "unusual number of genes in private orthogroups"
            )

        if not comments:
            comments.append(
                "no IQR outlier among current OrthoFinder metrics"
            )

        comments.append(
            "BUSCO and assembly-contiguity metrics pending"
        )

        return "; ".join(comments)

    qc["current_qc_interpretation"] = qc.apply(
        qc_comment,
        axis=1,
    )

    qc.to_csv(
        OUTPUT_DIR / "quality_control_summary.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # CLASS SUMMARY
    # -------------------------------------------------------------

    class_order = [
        "Core",
        "Soft core",
        "Accessory",
        "Private",
        "Reference-only",
    ]

    class_summary = (
        classification["orthogroup_class"]
        .value_counts()
        .reindex(
            class_order,
            fill_value=0,
        )
        .rename_axis("orthogroup_class")
        .reset_index(name="orthogroup_count")
    )

    class_summary["percentage_of_all_orthogroups"] = (
        class_summary["orthogroup_count"]
        / len(classification)
        * 100
    )

    class_summary.to_csv(
        OUTPUT_DIR / "orthogroup_class_summary.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------
    # PROCESSING REPORT
    # -------------------------------------------------------------

    overall_metrics = parse_overall_statistics(
        OVERALL_FILE
    )

    largest_private = private_summary.loc[
        private_summary["private_orthogroups"].idxmax()
    ]

    largest_unassigned = qc.loc[
        qc["unassigned_genes"].idxmax()
    ]

    report = f"""OrthoFinder PEP Step 3.3–3.4 processing report
=================================================

Input validation
----------------
Orthogroups: {len(membership):,}
Total proteins: {int(overall_metrics.get("Number of genes", 0)):,}
Assigned proteins: {int(overall_metrics.get("Number of genes in orthogroups", 0)):,}
Unassigned proteins: {int(overall_metrics.get("Number of unassigned genes", 0)):,}
Assignment rate: {overall_metrics.get("Percentage of genes in orthogroups", "NA")}%

Duplicate-column check
----------------------
No cultivar pair met the critical duplicate-column threshold.

All cultivar orthogroup classification
----------------------------------------
Core (7/7): {int((classes == "Core").sum()):,}
Soft core (6/7): {int((classes == "Soft core").sum()):,}
Accessory (2–5/7): {int((classes == "Accessory").sum()):,}
Private (1/7): {int((classes == "Private").sum()):,}
Reference-only: {int((classes == "Reference-only").sum()):,}

Preliminary QC observations
---------------------------
Cultivar with most private orthogroups:
{largest_private["cultivar"]} ({int(largest_private["private_orthogroups"])})

Cultivar with most unassigned proteins:
{largest_unassigned["cultivar"]} ({int(largest_unassigned["unassigned_genes"])})

Interpret these only as screening signals.
BUSCO and assembly-contiguity metrics remain to be added.
"""

    (
        OUTPUT_DIR / "README_processing_report.txt"
    ).write_text(
        report,
        encoding="utf-8",
    )

    # -------------------------------------------------------------
    # ZIP ALL OUTPUTS
    # -------------------------------------------------------------

    with zipfile.ZipFile(
        OUTPUT_ZIP,
        "w",
        zipfile.ZIP_DEFLATED,
    ) as archive:
        for file in sorted(OUTPUT_DIR.iterdir()):
            archive.write(
                file,
                arcname=file.name,
            )

    print(report)
    print(f"Outputs written to: {OUTPUT_DIR}")
    print(f"ZIP archive written to: {OUTPUT_ZIP}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)
