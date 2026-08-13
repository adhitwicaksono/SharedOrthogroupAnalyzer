# Shared Orthogroup Analyzer

**Repository owner:** Adhityo Wicaksono  
**Primary language:** Python  
**Original application:** comparative predicted-proteome analysis of seven Indonesian rice cultivars

Shared Orthogroup Analyzer is a small collection of Python scripts for converting OrthoFinder tabular outputs into orthogroup occupancy tables, cultivar-level quality-control summaries, and publication-ready visualizations.

These scripts were made for the Indonesian seven-cultivar whole-genome sequencing (WGS) project involving Boawae Seratus Malam (BSM), Cempo Ireng (CI), Gogo Jak (GJ), Inpago (IN), Kisol Manggarai (KM), Merah Pari Eja (MPE), and Putih Payo (PP), with the complete Nipponbare assembly AGIS1.0 used as the reference dataset. The current implementation therefore uses these dataset labels and a seven-cultivar occupancy scheme by default.

## What the repository does

The workflow takes OrthoFinder protein-orthogroup outputs and:

1. standardizes Galaxy/OrthoFinder dataset labels;
2. checks for duplicated or near-duplicated cultivar columns;
3. converts orthogroup membership into gene-count and presence/absence matrices;
4. classifies orthogroups by occupancy among the seven cultivars;
5. extracts core, soft-core, accessory, private, and reference-only orthogroups;
6. summarizes private orthogroups and unassigned proteins by cultivar;
7. combines these results with OrthoFinder per-species statistics and optional assembly/BUSCO metrics;
8. applies preliminary interquartile-range-based QC flags; and
9. generates a stacked-bar summary and a selected-intersection UpSet-style plot.

## Orthogroup definitions

AGIS1.0 is retained as a reference but is not counted when assigning cultivar occupancy classes.

| Class | Definition among BSM, CI, GJ, IN, KM, MPE, and PP |
|---|---|
| Core | Present in all seven cultivars |
| Soft core | Present in exactly six cultivars |
| Accessory | Present in two to five cultivars |
| Private | Present in exactly one cultivar |
| Reference-only | Absent from all seven cultivars but present in AGIS1.0 |

These classes describe **predicted orthogroup representation**. They do not, by themselves, demonstrate biological gene presence, absence, gain, loss, or cultivar specificity. Gene prediction errors, fragmented models, sequence divergence, and reference-guided genome construction can all affect orthogroup recovery.

The table-building scripts currently write the six-of-seven category as `Soft core`, whereas the figure labels use `Semi-core`. These terms refer to the same occupancy class in this repository.

## Scripts

| Script | Purpose |
|---|---|
| `build_orthogroup_tables_and_qc.py` | Reusable command-line workflow. Reads three required OrthoFinder tables and an optional cultivar-level QC table, performs duplicate-column checks, constructs orthogroup tables, and produces a QC summary. |
| `build_step3_tables_from_4_orthofinder_inputs.py` | Reproducibility snapshot for the original Indonesian seven-WGS analysis. Uses four hard-coded, corrected OrthoFinder PEP filenames from the 20 July 2026 run, generates the complete Step 3 table set and processing report, and packages the outputs as a ZIP archive. |
| `fig_ortho_stackbar.py` | Builds a stacked bar chart showing core, soft-core, accessory, private, and unassigned groups for AGIS1.0 and the seven cultivars. Also exports the plotted source data. |
| `fig_UpSet_orthogroups(1).py` | Builds an UpSet-style plot of selected exact orthogroup intersections among the seven cultivars, deliberately excluding AGIS1.0 from intersection definitions. Also exports the selected intersections. |

For new analyses, use `build_orthogroup_tables_and_qc.py`. The longer `build_step3_tables_from_4_orthofinder_inputs.py` is retained to reproduce the original project-specific run.

## Requirements

- Python 3.10 or later
- pandas
- NumPy
- Matplotlib

Install the required Python packages with:

```bash
python -m pip install pandas numpy matplotlib
```

## Expected inputs

### Reusable workflow

`build_orthogroup_tables_and_qc.py` requires:

- an OrthoFinder `Orthogroups.tsv`-like membership matrix;
- an OrthoFinder `Orthogroups_UnassignedGenes.tsv`-like matrix; and
- an OrthoFinder `Statistics_PerSpecies.tsv`-like file.

An optional tab-separated QC file may also be supplied. It must contain a `cultivar` column and can contain the following fields:

```text
busco_complete_pct
busco_single_copy_pct
busco_duplicated_pct
busco_fragmented_pct
busco_missing_pct
assembly_size_bp
contig_count
contig_n50_bp
contig_l50
longest_contig_bp
gc_pct
```

The orthogroup matrices must resolve to the following columns after label cleaning:

```text
Orthogroup, AGIS1.0, BSM, CI, GJ, IN, KM, MPE, PP
```

Galaxy-style names such as `Augustus on BSM_ Protein sequence` are converted automatically to their short labels.

### Project-specific reproducibility workflow

`build_step3_tables_from_4_orthofinder_inputs.py` expects the following files in the working directory:

```text
OrthoFinder_PEP_2026-07-20_corrected_Orthogroups.tsv
OrthoFinder_PEP_2026-07-20_corrected_UnassignedGenes.tsv
OrthoFinder_PEP_2026-07-20_corrected_Statistics_PerSpecies.tsv
OrthoFinder_PEP_2026-07-20_corrected_Statistics_Overall.tsv
```

## Usage

### 1. Run the reusable analysis

```bash
python build_orthogroup_tables_and_qc.py \
  --orthogroups Orthogroups.tsv \
  --unassigned Orthogroups_UnassignedGenes.tsv \
  --per-species-stats Statistics_PerSpecies.tsv \
  --output-dir results
```

With optional BUSCO and assembly metrics:

```bash
python build_orthogroup_tables_and_qc.py \
  --orthogroups Orthogroups.tsv \
  --unassigned Orthogroups_UnassignedGenes.tsv \
  --per-species-stats Statistics_PerSpecies.tsv \
  --qc-metrics cultivar_qc_metrics.tsv \
  --output-dir results
```

The script stops if duplicated or near-duplicated cultivar columns are detected. A bypass option exists for exceptional diagnostic use:

```bash
--allow-duplicate-columns
```

Bypassing this safeguard is not recommended for biological analysis.

### 2. Generate the stacked-bar figure

Run the figure script from the directory containing the generated category tables:

```bash
(cd results && python ../fig_ortho_stackbar.py)
```

### 3. Generate the UpSet-style figure

Run the intersection script from the directory containing `orthogroup_presence_absence.tsv`:

```bash
(cd results && python '../fig_UpSet_orthogroups(1).py')
```

### 4. Reproduce the original project-specific Step 3 run

Place the four corrected OrthoFinder files listed above beside the script, then run:

```bash
python build_step3_tables_from_4_orthofinder_inputs.py
```

This creates:

```text
OrthoFinder_PEP_2026-07-20_Step3_tables/
OrthoFinder_PEP_2026-07-20_Step3_tables.zip
```

## Main outputs

The table-building scripts generate all or most of the following files:

| Output | Contents |
|---|---|
| `orthogroup_presence_absence.tsv` | Binary presence/absence matrix across AGIS1.0 and the seven cultivars |
| `orthogroup_gene_counts.tsv` | Number of predicted proteins assigned to each orthogroup per dataset |
| `orthogroup_classification.tsv` | Occupancy class, reference status, and private-cultivar assignment for every orthogroup |
| `core_orthogroups.tsv` | Orthogroups represented in all seven cultivars |
| `softcore_orthogroups.tsv` | Orthogroups represented in six cultivars |
| `accessory_orthogroups.tsv` | Orthogroups represented in two to five cultivars |
| `private_orthogroups.tsv` | Orthogroups represented in one cultivar |
| `reference_only_orthogroups.tsv` | Orthogroups represented only in AGIS1.0; generated by the project-specific workflow |
| `private_orthogroups_by_cultivar.tsv` | Long-format private-orthogroup membership table |
| `private_orthogroup_summary_by_cultivar.tsv` | Private orthogroup and gene counts per cultivar; generated by the project-specific workflow |
| `unassigned_singleton_genes.tsv` | Long-format list of unassigned cultivar proteins |
| `orthogroup_class_summary.tsv` | Total number of orthogroups in each occupancy class |
| `cultivar_duplicate_column_check.tsv` or `cultivar_pair_column_similarity.tsv` | Pairwise duplicate/near-duplicate input check |
| `quality_control_summary.tsv` | OrthoFinder statistics, optional external QC metrics, and preliminary outlier flags |
| `README_processing_report.txt` | Plain-text processing summary from the project-specific workflow |

The visualization scripts additionally generate:

```text
Fig_A_stacked_bar_source_data_bottom_legend.tsv
Fig_A_stacked_bar_orthogroup.png
Fig_A_stacked_bar_orthogroup.pdf
Fig_B_UpSet_selected_intersections_v2.tsv
Fig_B_UpSet_orthogroups_v2.png
Fig_B_UpSet_orthogroups_v2.pdf
```

## Duplicate-column quality control

Duplicated sample columns can produce invalid core, accessory, and private classifications. The reusable workflow therefore compares all cultivar pairs and raises a critical error when either:

- at least 99.9% of orthogroup membership cells are identical; or
- at least 99.99% of per-orthogroup gene-count cells are identical.

The duplicate-column report is written before the analysis stops, allowing the problematic pair to be identified and the upstream OrthoFinder input to be corrected.

## Adapting the scripts to another dataset

The current release is customized for seven Indonesian rice cultivars and AGIS1.0. To apply it to another study, edit the dataset constants near the top of the scripts:

```python
CULTIVARS = ["BSM", "CI", "GJ", "IN", "KM", "MPE", "PP"]
REFERENCE_LABEL = "AGIS1.0"
```

The occupancy-classification rules must also be adjusted if the number of focal samples changes. Replacing only the names without updating the numerical thresholds will produce incorrect class assignments.

## Interpretation and limitations

- The scripts analyze predicted proteins and OrthoFinder assignments; they do not validate gene function.
- Private and reference-only orthogroups are candidates for follow-up, not confirmed biological presence/absence variants.
- Reference-guided consensus genomes can underrepresent non-reference sequence and structural variation.
- Fragmented, missed, or overpredicted gene models can alter apparent orthogroup occupancy.
- IQR-based flags are screening indicators, not formal evidence that a sample or biological result is invalid.
- Candidate private or missing loci should be checked using genomic read support, nucleotide-level synteny, local assembly, curated gene models, and transcript evidence where available.

## Artificial-intelligence use declaration

OpenAI ChatGPT 5.6 Sol (Thinking mode, Extra High reasoning effort, Fast speed) was used to assist with generation and refinement of the Python scripts, workflow logic, quality-control checks, visualization code, and repository documentation. The model worked from OrthoFinder tables, genome-analysis outputs, biological definitions, and analytical instructions supplied by Adhityo Wicaksono. It did not generate the WGS data, predicted proteomes, or OrthoFinder results and did not provide independent experimental validation. All scripts, thresholds, outputs, visualizations, interpretations, and documentation were reviewed and approved by Adhityo Wicaksono, who retains responsibility for the repository and its scientific use.

## Attribution

Developed and maintained by **Adhityo Wicaksono** for the comparative WGS analysis of seven Indonesian rice cultivars.

If this repository is used in another project, please cite the repository and the corresponding OrthoFinder publication, and clearly report any changes made to cultivar definitions, occupancy thresholds, or QC criteria.
