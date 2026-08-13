#!/usr/bin/env python3
"""
Build core/soft-core/accessory/private orthogroup tables and cultivar-level QC.

Required inputs:
  1) OrthoFinder Orthogroups.tsv-like matrix
  2) OrthoFinder Orthogroups_UnassignedGenes.tsv-like matrix
  3) OrthoFinder Statistics_PerSpecies.tsv-like file

Optional:
  4) A tab-separated QC metrics file containing cultivar-level BUSCO and
     assembly-contiguity metrics.

The script intentionally aborts if it detects duplicated cultivar columns,
because such duplication would invalidate core/accessory/private inference.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


CULTIVARS = ["BSM", "CI", "GJ", "IN", "KM", "MPE", "PP"]
REFERENCE_LABEL = "AGIS1.0"


def clean_species_name(column: str) -> str:
    """Convert Galaxy/OrthoFinder column labels to short species names."""
    if column == "Orthogroup":
        return column
    match = re.search(r"Augustus on (.+?)_ Protein sequence", column)
    if match:
        return match.group(1).strip()
    return column.strip()


def split_gene_ids(value) -> List[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split(",") if item.strip()]


def parse_orthogroup_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df = df.rename(columns={col: clean_species_name(col) for col in df.columns})

    required = ["Orthogroup", REFERENCE_LABEL, *CULTIVARS]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            "Missing required orthogroup columns: " + ", ".join(missing)
        )

    return df[required].copy()


def duplicate_column_checks(membership: pd.DataFrame) -> pd.DataFrame:
    """Detect identical or near-identical cultivar membership columns."""
    rows = []
    n = len(membership)

    for i, a in enumerate(CULTIVARS):
        for b in CULTIVARS[i + 1:]:
            identical_cells = int((membership[a] == membership[b]).sum())
            identical_fraction = identical_cells / n if n else 0.0

            count_a = membership[a].map(lambda x: len(split_gene_ids(x)))
            count_b = membership[b].map(lambda x: len(split_gene_ids(x)))
            identical_count_cells = int((count_a == count_b).sum())
            identical_count_fraction = identical_count_cells / n if n else 0.0

            rows.append(
                {
                    "cultivar_a": a,
                    "cultivar_b": b,
                    "identical_membership_cells": identical_cells,
                    "total_orthogroups": n,
                    "identical_membership_fraction": identical_fraction,
                    "identical_gene_count_cells": identical_count_cells,
                    "identical_gene_count_fraction": identical_count_fraction,
                    "critical_duplicate_flag": (
                        identical_fraction >= 0.999
                        or identical_count_fraction >= 0.9999
                    ),
                }
            )

    return pd.DataFrame(rows)


def parse_per_species_stats(path: Path) -> pd.DataFrame:
    """
    Parse the first summary block of Statistics_PerSpecies.tsv.
    This file contains several tables; only rows 1-11 are needed.
    """
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 11:
        raise ValueError("Per-species statistics file is unexpectedly short.")

    header = lines[0].split("\t")
    species = [clean_species_name(x) for x in header[1:]]

    wanted = {
        "Number of genes": "total_predicted_proteins",
        "Number of genes in orthogroups": "genes_assigned_to_orthogroups",
        "Number of unassigned genes": "unassigned_genes",
        "Percentage of genes in orthogroups": "percentage_assigned",
        "Percentage of unassigned genes": "percentage_unassigned",
        "Number of orthogroups containing species": "orthogroups_containing_cultivar",
        "Number of species-specific orthogroups": "orthofinder_species_specific_orthogroups",
        "Number of genes in species-specific orthogroups": "genes_in_orthofinder_species_specific_orthogroups",
    }

    records: Dict[str, Dict[str, object]] = {
        sp: {"cultivar": sp} for sp in species
    }

    for line in lines[1:12]:
        fields = line.split("\t")
        metric = fields[0]
        if metric not in wanted:
            continue
        out_name = wanted[metric]
        for sp, value in zip(species, fields[1:]):
            records[sp][out_name] = value

    df = pd.DataFrame(records.values())
    numeric_cols = [c for c in df.columns if c != "cultivar"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_optional_qc(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    qc = pd.read_csv(path, sep="\t")
    if "cultivar" not in qc.columns:
        raise ValueError("Optional QC metrics file must contain a 'cultivar' column.")
    return qc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--unassigned", required=True, type=Path)
    parser.add_argument("--per-species-stats", required=True, type=Path)
    parser.add_argument("--qc-metrics", type=Path, default=None)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--allow-duplicate-columns",
        action="store_true",
        help="Not recommended. Continue despite duplicate-column detection.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    membership = parse_orthogroup_matrix(args.orthogroups)

    # ---------- Hard QC gate: duplicated cultivar inputs ----------
    duplicate_report = duplicate_column_checks(membership)
    duplicate_report.to_csv(
        args.output_dir / "cultivar_duplicate_column_check.tsv",
        sep="\t",
        index=False,
    )

    critical = duplicate_report[
        duplicate_report["critical_duplicate_flag"] == True
    ]

    if not critical.empty and not args.allow_duplicate_columns:
        pairs = ", ".join(
            f"{r.cultivar_a}-{r.cultivar_b}"
            for r in critical.itertuples(index=False)
        )
        raise RuntimeError(
            "CRITICAL: duplicated or near-duplicated cultivar columns detected: "
            f"{pairs}. No biological output tables were generated. "
            "Correct the OrthoFinder input and rerun."
        )

    # ---------- Gene-count matrix ----------
    gene_counts = membership.copy()
    for species in [REFERENCE_LABEL, *CULTIVARS]:
        gene_counts[species] = gene_counts[species].map(
            lambda x: len(split_gene_ids(x))
        )

    gene_counts.to_csv(
        args.output_dir / "orthogroup_gene_counts.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Presence/absence matrix ----------
    presence = gene_counts.copy()
    for species in [REFERENCE_LABEL, *CULTIVARS]:
        presence[species] = (presence[species] > 0).astype(int)

    presence.to_csv(
        args.output_dir / "orthogroup_presence_absence.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Classification among eight cultivars ----------
    classification = gene_counts.copy()
    classification["cultivar_occupancy"] = presence[CULTIVARS].sum(axis=1)
    classification["reference_present"] = presence[REFERENCE_LABEL]

    classification["orthogroup_class"] = classification[
        "cultivar_occupancy"
    ].map(
        lambda n:
            "Core" if n == 7 else
            "Soft core" if n == 6 else
            "Accessory" if 2 <= n <= 5 else
            "Private" if n == 1 else
            "Reference-only" if n == 0 else
            "Unclassified"
    )

    private_owner = []
    for _, row in presence.iterrows():
        present_cultivars = [c for c in CULTIVARS if row[c] == 1]
        private_owner.append(
            present_cultivars[0] if len(present_cultivars) == 1 else ""
        )
    classification["private_cultivar"] = private_owner

    classification.to_csv(
        args.output_dir / "orthogroup_classification.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Category-specific tables ----------
    category_to_filename = {
        "Core": "core_orthogroups.tsv",
        "Soft core": "softcore_orthogroups.tsv",
        "Accessory": "accessory_orthogroups.tsv",
        "Private": "private_orthogroups.tsv",
    }

    for category, filename in category_to_filename.items():
        ids = classification.loc[
            classification["orthogroup_class"] == category,
            ["Orthogroup", "cultivar_occupancy", "reference_present",
             "private_cultivar"],
        ]
        category_membership = membership.merge(ids, on="Orthogroup", how="inner")
        category_membership.to_csv(
            args.output_dir / filename,
            sep="\t",
            index=False,
        )

    # ---------- Private orthogroups by cultivar, long format ----------
    private_ids = classification.loc[
        classification["orthogroup_class"] == "Private",
        ["Orthogroup", "private_cultivar", "reference_present"],
    ]

    private_rows = []
    for row in private_ids.itertuples(index=False):
        gene_text = membership.loc[
            membership["Orthogroup"] == row.Orthogroup,
            row.private_cultivar,
        ].iloc[0]
        gene_ids = split_gene_ids(gene_text)
        private_rows.append(
            {
                "cultivar": row.private_cultivar,
                "orthogroup": row.Orthogroup,
                "gene_count": len(gene_ids),
                "gene_ids": ", ".join(gene_ids),
                "reference_present": row.reference_present,
            }
        )

    private_long = pd.DataFrame(private_rows).sort_values(
        ["cultivar", "orthogroup"]
    )
    private_long.to_csv(
        args.output_dir / "private_orthogroups_by_cultivar.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Unassigned singleton genes ----------
    unassigned = parse_orthogroup_matrix(args.unassigned)
    unassigned_rows = []

    for row in unassigned.itertuples(index=False):
        row_dict = row._asdict()
        orthogroup = row_dict["Orthogroup"]
        for cultivar in CULTIVARS:
            for gene_id in split_gene_ids(row_dict[cultivar]):
                unassigned_rows.append(
                    {
                        "cultivar": cultivar,
                        "orthogroup_label": orthogroup,
                        "gene_id": gene_id,
                    }
                )

    unassigned_long = pd.DataFrame(unassigned_rows)
    if not unassigned_long.empty:
        unassigned_long = unassigned_long.sort_values(
            ["cultivar", "orthogroup_label", "gene_id"]
        )
    unassigned_long.to_csv(
        args.output_dir / "unassigned_singleton_genes.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Summary counts ----------
    class_summary = (
        classification["orthogroup_class"]
        .value_counts()
        .rename_axis("orthogroup_class")
        .reset_index(name="orthogroup_count")
    )
    class_summary.to_csv(
        args.output_dir / "orthogroup_class_summary.tsv",
        sep="\t",
        index=False,
    )

    # ---------- Cultivar-level QC ----------
    species_stats = parse_per_species_stats(args.per_species_stats)
    species_stats = species_stats[
        species_stats["cultivar"].isin(CULTIVARS)
    ].copy()

    private_summary = (
        private_long.groupby("cultivar")
        .agg(
            private_orthogroups=("orthogroup", "nunique"),
            genes_in_private_orthogroups=("gene_count", "sum"),
        )
        .reset_index()
    )

    unassigned_summary = (
        unassigned_long.groupby("cultivar")
        .agg(unassigned_singleton_rows=("gene_id", "count"))
        .reset_index()
        if not unassigned_long.empty
        else pd.DataFrame({"cultivar": CULTIVARS,
                           "unassigned_singleton_rows": 0})
    )

    qc = (
        species_stats.merge(private_summary, on="cultivar", how="left")
        .merge(unassigned_summary, on="cultivar", how="left")
    )

    optional_qc = load_optional_qc(args.qc_metrics)
    if optional_qc is not None:
        qc = qc.merge(optional_qc, on="cultivar", how="left")
    else:
        for col in [
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
        ]:
            qc[col] = pd.NA

    # Simple robust outlier flags based on IQR where possible.
    for metric in [
        "total_predicted_proteins",
        "unassigned_genes",
        "private_orthogroups",
        "genes_in_private_orthogroups",
        "busco_fragmented_pct",
        "contig_count",
        "contig_n50_bp",
    ]:
        if metric not in qc.columns:
            continue
        values = pd.to_numeric(qc[metric], errors="coerce")
        if values.notna().sum() >= 4:
            q1 = values.quantile(0.25)
            q3 = values.quantile(0.75)
            iqr = q3 - q1
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            qc[f"{metric}_outlier_flag"] = (
                (values < lower) | (values > upper)
            )
        else:
            qc[f"{metric}_outlier_flag"] = pd.NA

    def combine_flags(row) -> str:
        concerns = []
        if row.get("unassigned_genes_outlier_flag") is True:
            concerns.append("unusually many/few unassigned genes")
        if row.get("private_orthogroups_outlier_flag") is True:
            concerns.append("unusual private-orthogroup count")
        if row.get("total_predicted_proteins_outlier_flag") is True:
            concerns.append("unusual predicted-protein count")
        if row.get("busco_fragmented_pct_outlier_flag") is True:
            concerns.append("unusual BUSCO fragmentation")
        if row.get("contig_count_outlier_flag") is True:
            concerns.append("unusual contig count")
        if row.get("contig_n50_bp_outlier_flag") is True:
            concerns.append("unusual contig N50")
        return "; ".join(concerns) if concerns else "No automatic outlier flag"

    qc["automatic_qc_comment"] = qc.apply(combine_flags, axis=1)

    qc.to_csv(
        args.output_dir / "quality_control_summary.tsv",
        sep="\t",
        index=False,
    )

    print("Completed successfully.")
    print(f"Outputs written to: {args.output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)