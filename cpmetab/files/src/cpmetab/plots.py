"""Figures, drawn with the same code and styling as the original scripts.
Every function returns a matplotlib Figure (nothing is shown or saved here)."""
from __future__ import annotations

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

FONT = {"fontname": "Arial"}


def mpl_defaults():
    """Style of the scripts that never called sns.set (Figures 3B, 6, 7b)."""
    import matplotlib
    matplotlib.rcdefaults()


def seaborn_white():
    """Style after ``sns.set(style="white")`` (cell_count_version_4.py, Figure 7a)."""
    sns.set(style="white")


# --------------------------------------------------------- cell-count heatmap
def toxicity_heatmap(matrix: pd.DataFrame, cell_line_order) -> Figure:
    """``plot_clustermap`` (cell_count_version_4.py)."""
    masked = matrix.apply(lambda c: c.map(lambda x: min(x, 100)))
    masked = masked.reindex([c for c in cell_line_order if c in masked.index])
    figsize = (len(masked.columns) * 1, len(masked.index) * 1.5)
    seaborn_white()
    g = sns.clustermap(masked, cmap="coolwarm", linewidths=0.5, figsize=figsize,
                       col_cluster=False, row_cluster=False, dendrogram_ratio=(.1, .2),
                       cbar_kws={"label": "Cell Count % (DMSO=100%)"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), rotation=90, ha='right',
                                 fontsize=50)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), rotation=0, fontsize=50)
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    return g.figure


# ---------------------------------------------------------------- ED heatmaps
def _cbar(ax):
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Euclidean Distance', fontsize=18, **FONT)
    cbar.ax.yaxis.labelpad = 15


def ed_heatmap(matrix: pd.DataFrame, variant: str, vmax: float = 6.9) -> Figure:
    """``create_heatmap`` ("heatmap1"), ``create_heatmap_2`` ("heatmap2") or
    ``create_heatmap_with_order`` ("heatmap_ordered")."""
    size = {"heatmap1": (6, 18), "heatmap2": (6, 10), "heatmap_ordered": (8, 10)}[variant]
    seaborn_white()
    fig = plt.figure(figsize=size)
    ax = sns.heatmap(matrix, cmap='BuPu', linewidths=.5, vmin=0, vmax=vmax,
                     cbar_kws={'label': 'Euclidean Distance'})
    _cbar(ax)
    ysize = 16 if variant == "heatmap1" else 18
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=ysize, rotation=0, **FONT)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=18, rotation=90, ha='center', **FONT)
    ax.set_xlabel('')
    ax.set_ylabel('')
    if variant == "heatmap_ordered":
        plt.tight_layout()
    return fig


# ------------------------------------------------------- correlation boxplot
def correlation_boxplot(corr: pd.DataFrame, cell_line_order) -> Figure:
    """``plot_combined_boxplot`` (cell_count_version_4.py)."""
    d = corr.rename(columns={"r": "correlation"})
    present = [c for c in cell_line_order if c in set(d['cell-line'])]
    d = d[d['cell-line'].isin(present)].copy()
    d['cell-line'] = pd.Categorical(d['cell-line'], categories=present, ordered=True)
    seaborn_white()
    fig = plt.figure(figsize=(6, 4))
    plt.gca().set_facecolor('#f7f7f7')
    sns.boxplot(x='cell-line', y='correlation', data=d, hue='cell-line', palette='Set2',
                legend=False)
    plt.xlabel('Cell Line', fontsize=16)
    plt.ylabel('Pearson Correlation Coefficient', fontsize=16)
    plt.ylim(0, 1.1)
    plt.yticks(fontsize=14)
    plt.grid(axis='y', linestyle='--', linewidth=0.7, color='gray', alpha=0.7)
    plt.xticks(rotation=45, ha='right', fontsize=14)
    return fig


# ------------------------------- Controls: technical/biological replicate validation
def control_ed_plot(ctrl: pd.DataFrame, control: str, cell_line_order) -> Figure:
    """Negative/positive controls: Euclidean distance of each well's EMD
    (Wasserstein) profile, per cell line and biological replicate
    (figure3b_control_ed.py). One figure per control, broken y-axis."""
    from brokenaxes import brokenaxes
    mpl_defaults()
    d = ctrl[ctrl["control"] == control].copy()
    lines = [c for c in cell_line_order if c in set(d["cell-line"])]
    reps = sorted(ctrl["replicate"].unique(), key=lambda r: int(r) if str(r).isdigit() else r)
    fig = plt.figure(figsize=(max(18 / 7 * len(lines), 3.5), 4))
    gs = GridSpec(1, len(lines), width_ratios=[1] * len(lines), figure=fig)
    axes = [brokenaxes(ylims=((0, 6), (6, 8)), height_ratios=[0.25, 0.75], hspace=0, d=0,
                       subplot_spec=gs[i]) for i in range(len(lines))]
    d["biological_replicate"] = d["replicate"].map({r: i for i, r in enumerate(reps)})
    for i, cl in enumerate(lines):
        e = d[d["cell-line"] == cl]
        if e.empty:
            continue
        stats = e.groupby('biological_replicate')['ed'].agg(
            median_emd='median', q25_emd=lambda x: x.quantile(0.25),
            q75_emd=lambda x: x.quantile(0.75)).reset_index()
        for replicate in range(len(reps)):
            row = stats[stats['biological_replicate'] == replicate]
            if row.empty:
                continue
            row = row.iloc[0]
            axes[i].errorbar(x=replicate, y=row['median_emd'],
                             yerr=[[row['median_emd'] - row['q25_emd']],
                                   [row['q75_emd'] - row['median_emd']]],
                             fmt='', ecolor='black', capsize=5, elinewidth=2)
            axes[i].plot([replicate - 0.3, replicate + 0.3], [row['median_emd']] * 2,
                         color='black', linestyle='-', linewidth=2)
        well_col = "well-id" if "well-id" in e else None
        groups = e.groupby(well_col) if well_col else [(None, e)]
        for well, w in groups:
            jitter = 0.05 * np.random.randn(len(w))
            axes[i].scatter(w['biological_replicate'] + jitter, w['ed'], edgecolors='blue',
                            facecolors='none', label=well, alpha=0.7)
        for ax in axes[i].axs:      # directly on the sub-axes: works for any number of replicates
            ax.set_xlim(-0.5, len(reps) - 0.5)
            ax.set_xticks(list(range(len(reps))))
            ax.set_xticklabels([f'Biological Rep {r}' for r in reps],
                               rotation=45, ha='right', fontsize=12, **FONT)
        axes[i].set_xlabel(f'{cl}', labelpad=80, fontsize=14, **FONT)
        axes[i].axs[1].set_yticks([0, 1, 2, 3, 4, 5])
        axes[i].axs[0].set_yticks([6, 7])
        axes[i].set_ylabel('Euclidean Distance' if i == 0 else '', fontsize=12)
    fig.suptitle(f'{control}', fontsize=16, y=0.9, **FONT)
    fig.subplots_adjust(hspace=0.05, wspace=0.5)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        if lines:   # the original plt.legend() call: an empty legend box on the last panel
            axes[len(lines) - 1].axs[0].legend(handles=[])
    return fig


# --------------------- EMD feature deviation counts between biological replicates
def outlier_bars(counts: pd.DataFrame, replicate, cell_line: str, treatments,
                 show_ylabel: bool = False, show_legend: bool = False) -> Figure:
    """Number of features whose signed Wasserstein distance to replicate 1
    deviates (figure6_step4_plot.py); one figure per cell line and replicate."""
    from brokenaxes import brokenaxes
    mpl_defaults()
    total_features = 405
    cell_df = counts[(counts['cell-line'] == cell_line)
                     & (counts['replicate'].astype(str) == str(replicate))]
    low, up = [], []
    for t in treatments:
        td = cell_df[cell_df['treatment'] == t]
        low.append((td['low_count'].sum() / total_features) * total_features)
        up.append((td['up_count'].sum() / total_features) * total_features)
    fig = plt.figure(figsize=(4, 4))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = brokenaxes(ylims=((0, 180), (350, 405)), height_ratios=[0.2, 0.8], hspace=0.2,
                    d=0.01, subplot_spec=gs[0])
    x = np.arange(len(treatments))
    ax.bar(x, low, 0.6, color='#4472C4', edgecolor='black', linewidth=1.5,
           label='Threshold -0.15')
    ax.bar(x, up, 0.6, bottom=low, color='#C00000', edgecolor='black', linewidth=1.5,
           label='Threshold 0.15')
    ax.set_xticks(x, minor=True)
    ax.set_xticklabels(list(treatments), minor=True, rotation=45, ha='right', fontsize=18, **FONT)
    ax.set_xticks(x)
    ax.set_xticklabels(np.append('', list(treatments)), rotation=45, ha='right', fontsize=18,
                       **FONT)
    ax.set_xlabel(f'{cell_line}', labelpad=100, fontsize=20, **FONT)
    if show_ylabel:
        ax.set_ylabel('Number of outliers', fontsize=20, labelpad=40, **FONT)
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    ax.axs[1].set_yticks([0, 40, 80, 120, 160, 200])
    ax.axs[0].set_yticks([350, 405])
    ax.axs[1].set_yticklabels([0, 40, 80, 120, 160, 200], fontsize=18, **FONT)
    ax.axs[0].set_yticklabels([350, 405], fontsize=18, **FONT)
    if show_legend:
        legend = ax.legend(handles=[
            Patch(facecolor='#4472C4', edgecolor='black', linewidth=2, label='IQR threshold -0.15'),
            Patch(facecolor='#C00000', edgecolor='black', linewidth=2, label='IQR threshold 0.15')],
            loc='upper right', fontsize=14, frameon=True)
        for text in legend.get_texts():
            text.set_fontname('Arial')
    return fig
