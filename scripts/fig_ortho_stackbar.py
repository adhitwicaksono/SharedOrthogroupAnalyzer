
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

files = {
    'Core': 'core_orthogroups.tsv',
    'Semi-core': 'softcore_orthogroups.tsv',
    'Accessory': 'accessory_orthogroups.tsv',
    'Private': 'private_orthogroups.tsv',
    'Unassigned singletons': 'unassigned_singleton_genes.tsv',
}

cultivars = ['AGIS1.0', 'BSM', 'CI', 'GJ', 'IN', 'KM', 'MPE', 'PP']


def read_tsv(path):
    return pd.read_csv(path, sep='	')


def count_presence_by_cultivar(df, cultivar_columns):
    counts = {}
    for col in cultivar_columns:
        if col not in df.columns:
            counts[col] = 0
            continue
        present = df[col].notna() & (df[col].astype(str).str.strip() != '')
        counts[col] = int(present.sum())
    return counts

core_df = read_tsv(files['Core'])
softcore_df = read_tsv(files['Semi-core'])
accessory_df = read_tsv(files['Accessory'])
private_df = read_tsv(files['Private'])
unassigned_df = read_tsv(files['Unassigned singletons'])

summary = pd.DataFrame(index=cultivars)
summary['Core'] = pd.Series(count_presence_by_cultivar(core_df, cultivars))
summary['Semi-core'] = pd.Series(count_presence_by_cultivar(softcore_df, cultivars))
summary['Accessory'] = pd.Series(count_presence_by_cultivar(accessory_df, cultivars))
summary['Private'] = pd.Series(count_presence_by_cultivar(private_df, cultivars))

unassigned_counts = (
    unassigned_df.groupby('cultivar')['gene_id']
    .nunique()
    .reindex(cultivars, fill_value=0)
)
summary['Unassigned singletons'] = unassigned_counts
summary = summary.fillna(0).astype(int)
summary.insert(0, 'Cultivar', summary.index)

source_path = 'Fig_A_stacked_bar_source_data_bottom_legend.tsv'
summary.to_csv(source_path, sep='	', index=False)

fig, ax = plt.subplots(figsize=(12.5, 8))
categories = ['Core', 'Semi-core', 'Accessory', 'Private', 'Unassigned singletons']
bottom = pd.Series([0] * len(summary), index=summary.index)

for cat in categories:
    ax.bar(summary['Cultivar'], summary[cat], bottom=bottom.values, label=cat, edgecolor='black', linewidth=0.4)
    bottom += summary[cat]

ax.set_xlabel('Cultivar / reference', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{int(x):,}'))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for tick in ax.get_xticklabels():
    tick.set_rotation(25)
    tick.set_ha('right')

# Legend at the bottom, centered, no overlap
ax.legend(
    frameon=False,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.18),
    ncol=3,
    title='Category'
)

fig.text(
    0.5, 0.02,
    "Note: 'Unassigned singletons' in the supplied table were available only for the eight local cultivars; AGIS1.0 is shown as 0 in that category.",
    ha='center', fontsize=9
)

fig.tight_layout(rect=[0, 0.12, 1, 1])

out_png = 'Fig_A_stacked_bar_orthogroup.png'
out_pdf = 'Fig_A_stacked_bar_orthogroup.pdf'
fig.savefig(out_png, dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
