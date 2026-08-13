import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np

PRESENCE_ABSENCE_FILE = 'orthogroup_presence_absence.tsv'
OUT_PNG = 'Fig_B_UpSet_orthogroups_v2.png'
OUT_PDF = 'Fig_B_UpSet_orthogroups_v2.pdf'
OUT_SELECTED = 'Fig_B_UpSet_selected_intersections_v2.tsv'

# Seven Indonesian cultivars (AGIS1.0 excluded here on purpose)
CULTIVARS = ['BSM', 'CI', 'GJ', 'IN', 'KM', 'MPE', 'PP']
N_ACCESSORY = 10  # number of larger accessory intersections to keep in the main figure

pa = pd.read_csv(PRESENCE_ABSENCE_FILE, sep='\t')

# Build exact intersection counts
agg = pa.groupby(CULTIVARS, dropna=False).size().reset_index(name='count')
agg['members'] = agg.apply(lambda r: ';'.join([c for c in CULTIVARS if int(r[c]) == 1]), axis=1)
agg['occupancy'] = agg[CULTIVARS].sum(axis=1)
agg['category'] = np.where(
    agg['occupancy'] == 7, 'Core (7/7)',
    np.where(
        agg['occupancy'] == 6, 'Semi-core',
        np.where(
            agg['occupancy'] == 1, 'Private',
            np.where(agg['occupancy'].between(2, 5), 'Accessory', 'Absent')
        )
    )
)
full_df = agg[['members', 'count', 'occupancy', 'category']].sort_values(['occupancy', 'count'], ascending=[False, False])

# Select intersections for plotting
core_df = full_df[full_df['occupancy'] == 8].copy()
seven_df = full_df[full_df['occupancy'] == 7].sort_values('count', ascending=False).copy()
acc_df = full_df[full_df['occupancy'].between(2, 6)].sort_values(['occupancy', 'count'], ascending=[False, False]).head(N_ACCESSORY).copy()
private_df = full_df[full_df['occupancy'] == 1].copy()
private_df['cultivar'] = private_df['members']
private_df = private_df.set_index('cultivar').reindex(CULTIVARS).reset_index().rename(columns={'index': 'cultivar'})

selected_rows = []
selected_rows.extend(core_df.to_dict('records'))
selected_rows.extend(seven_df.to_dict('records'))
selected_rows.extend(acc_df.to_dict('records'))
selected_rows.extend(private_df.to_dict('records'))
selected_df = pd.DataFrame(selected_rows).drop_duplicates(subset=['members']).reset_index(drop=True)
selected_df['plot_group'] = np.select(
    [
        selected_df['occupancy'] == 7,
        selected_df['occupancy'] == 6,
        selected_df['occupancy'].between(2, 5),
        selected_df['occupancy'] == 1,
    ],
    ['Core (7/7)', 'Semi-core', 'Accessory', 'Private'],
    default='Other'
)
plot_order = core_df['members'].tolist() + seven_df['members'].tolist() + acc_df['members'].tolist() + private_df['members'].tolist()
selected_df['plot_order'] = selected_df['members'].map({m: i for i, m in enumerate(plot_order)})
selected_df = selected_df.sort_values('plot_order').reset_index(drop=True)

# Short, clean labels for the x-axis
labels = []
acc_n = 0
for _, row in selected_df.iterrows():
    members = row['members'].split(';') if row['members'] else []
    if row['occupancy'] == 7:
        lbl = 'All 7'
    elif row['occupancy'] == 6:
        missing = [c for c in CULTIVARS if c not in members][0]
        lbl = f'no {missing}'
    elif row['occupancy'] == 1:
        lbl = members[0]
    else:
        acc_n += 1
        lbl = f'Acc{acc_n}'
    labels.append(lbl)
selected_df['plot_label'] = labels
selected_df.to_csv(OUT_SELECTED, sep='\t', index=False)

# Set sizes per cultivar
set_sizes = {c: int(pa[c].sum()) for c in CULTIVARS}

# Plotting
plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(16, 8.5))
gs = GridSpec(2, 2, width_ratios=[1.8, 6.2], height_ratios=[3.0, 2.0], hspace=0.08, wspace=0.08)
ax_top = fig.add_subplot(gs[0, 1])
ax_mat = fig.add_subplot(gs[1, 1])
ax_set = fig.add_subplot(gs[1, 0], sharey=ax_mat)
fig.add_subplot(gs[0, 0]).axis('off')

color_map = {
    'Core (7/7)': '#1f77b4',
    'Semi-core': '#ff7f0e',
    'Accessory': '#2ca02c',
    'Private': '#d62728'
}

x = np.arange(len(selected_df))
bar_colors = [color_map[g] for g in selected_df['plot_group']]
bars = ax_top.bar(x, selected_df['count'], color=bar_colors, width=0.8)
ax_top.set_ylabel('Orthogroup count')
ax_top.set_xticks([])
ax_top.spines[['top', 'right']].set_visible(False)

for grp in ['Core (7/7)', 'Semi-core', 'Accessory', 'Private']:
    idx = selected_df.index[selected_df['plot_group'] == grp].tolist()
    if idx:
        left = min(idx) - 0.5
        right = max(idx) + 0.5
        mid = (left + right) / 2
        ax_top.axvline(left, color='lightgray', lw=1, ls='--', zorder=0)
        ax_mat.axvline(left, color='lightgray', lw=1, ls='--', zorder=0)
        ax_top.text(mid, ax_top.get_ylim()[1] * 1.01, grp.replace(' (8/8)', ''), ha='center', va='bottom', fontsize=10)
ax_top.axvline(len(selected_df) - 0.5, color='lightgray', lw=1, ls='--', zorder=0)
ax_mat.axvline(len(selected_df) - 0.5, color='lightgray', lw=1, ls='--', zorder=0)

for rect, val in zip(bars, selected_df['count']):
    if val >= 0:
        ax_top.text(rect.get_x() + rect.get_width() / 2,
                    rect.get_height() + max(selected_df['count']) * 0.008,
                    f'{int(val)}', ha='center', va='bottom', rotation=90, fontsize=7)

for i in range(len(CULTIVARS)):
    if i % 2 == 0:
        ax_mat.axhspan(i - 0.5, i + 0.5, color='#f5f5f5', zorder=0)

for xi, row in selected_df.iterrows():
    members = row['members'].split(';') if row['members'] else []
    y_present = [CULTIVARS.index(c) for c in members]
    y_absent = [i for i, c in enumerate(CULTIVARS) if c not in members]
    if y_absent:
        ax_mat.scatter([xi] * len(y_absent), y_absent, s=22, color='lightgray', zorder=2)
    if y_present:
        ax_mat.scatter([xi] * len(y_present), y_present, s=55, color='black', zorder=3)
        if len(y_present) > 1:
            ax_mat.plot([xi, xi], [min(y_present), max(y_present)], color='black', lw=1.2, zorder=1)

ax_mat.set_xlim(-0.6, len(selected_df) - 0.4)
ax_mat.set_ylim(-0.5, len(CULTIVARS) - 0.5)
ax_mat.set_yticks(range(len(CULTIVARS)))
ax_mat.set_yticklabels(CULTIVARS)
ax_mat.set_xticks(x)
ax_mat.set_xticklabels(selected_df['plot_label'], rotation=90)
ax_mat.spines[['top', 'right', 'bottom']].set_visible(False)
ax_mat.set_xlabel('Selected exact orthogroup intersections')

sizes = [set_sizes[c] for c in CULTIVARS]
y = np.arange(len(CULTIVARS))
ax_set.barh(y, sizes, color='#4c78a8')
ax_set.invert_xaxis()
ax_set.set_xlabel('Total orthogroups\nper cultivar')
ax_set.tick_params(axis='y', left=False, labelleft=False)
ax_set.spines[['top', 'right', 'left']].set_visible(False)
for yy, val in zip(y, sizes):
    ax_set.text(val * 0.98, yy, f'{val}', va='center', ha='right', fontsize=8)

legend = [
    Line2D([0], [0], color=color_map['Core (7/7)'], lw=8, label='Core (shared by all 7)'),
    Line2D([0], [0], color=color_map['Semi-core'], lw=8, label='Semi-core (shared by 6)'),
    Line2D([0], [0], color=color_map['Accessory'], lw=8, label='Accessory (between 2 and 5)'),
    Line2D([0], [0], color=color_map['Private'], lw=8, label='Private (just 1)')
]
fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(0.5, -0.005), ncol=2, frameon=False)
fig.subplots_adjust(top=0.92, bottom=0.15, left=0.07, right=0.98)
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
fig.savefig(OUT_PDF, bbox_inches='tight')
print('Saved:', OUT_PNG)
print('Saved:', OUT_PDF)
print('Saved:', OUT_SELECTED)
