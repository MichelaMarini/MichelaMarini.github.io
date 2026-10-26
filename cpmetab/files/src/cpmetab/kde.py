"""Figure 7, as in original_scripts/figure7a_kde_replicates.py and figure7b_ridgeplot.py.

Values are min-max scaled to [0, 1] per curve group (each replicate in 7a,
each treatment in 7b), as in the original scripts.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from .config import canon

FIG7A = {"cell_line": "U2OS", "feature": "Shape_Area_cell", "treatment": "DMSO",
         "color": "#74B6A1", "xlabel": "Area cell"}
FIG7B = {"feature": "Shape_Area_cyto", "replicate": "1", "xlabel": "Cell area",
         "treatments": ("DMSO", "Toyocamycin", "17-dmag", "Tms", "Artemisinin", "Fccp")}


def _minmax(s: pd.Series) -> pd.Series:
    """sklearn MinMaxScaler on one column."""
    lo, hi = s.min(), s.max()
    return (s - lo) / ((hi - lo) or 1.0)


def _style(ax, xlabel):
    ax.set_xlabel(xlabel, fontsize=18, fontname="Arial")
    ax.set_ylabel("Density curves (per well, FOV)", fontsize=18, fontname="Arial")
    ax.set_ylim(0, 14)
    ax.set_yticks(range(0, 15, 2))
    ax.set_xticks([0, 0.5, 1])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def replicate_density_figures(cells: pd.DataFrame, cell_line: str = FIG7A["cell_line"],
                              feature: str = FIG7A["feature"],
                              treatment: str = FIG7A["treatment"], label: str = "DMSO",
                              color: str = FIG7A["color"], xlabel: str = FIG7A["xlabel"]
                              ) -> dict[str, Figure]:
    """Figure 7a: per replicate (wells/FOV + global), overlay, and globals only."""
    d = cells[(cells["cell-line"] == cell_line)
              & (cells["treatment"].map(canon) == canon(treatment))]
    if d.empty:
        raise ValueError(f"No cells of {treatment} in {cell_line}")
    reps = sorted(d["replicate"].unique(), key=lambda r: int(r) if str(r).isdigit() else r)
    reps_df = {}
    for r in reps:
        f = d[d["replicate"] == r].copy()
        f[feature] = _minmax(f[feature].astype(float))
        reps_df[r] = f
    figs = {}
    from .plots import seaborn_white
    seaborn_white()
    for idx, r in enumerate(reps):
        f = reps_df[r]
        fig = plt.figure(figsize=(7, 6))
        for _, g in f.groupby(['well-id', 'fov']):
            sns.kdeplot(g[feature], lw=1, alpha=0.9, color=color)
        sns.kdeplot(f[feature], color="black", lw=2)
        ax = plt.gca()
        ax.legend(handles=[Line2D([0], [0], color=color, lw=2, label=f"{label} wells"),
                           Line2D([0], [0], color="black", lw=2, label=f"Global {label}")],
                  loc="best", fontsize=16)
        ax.set_title(f"Biological replicate {r}", fontsize=22, fontname="Arial")
        _style(ax, xlabel)
        figs[f"rep{r}"] = fig
    if len(reps) >= 2:
        rep_colors = ["#74C6AC", "#E99675", "#95A3C3", "#DB96C0", "#A2C865", "#E5C949", "#DBC29E"]
        styles = ["-", "--", ":", "-.", (0, (5, 1)), (0, (1, 3)), (0, (3, 1, 1, 1))]
        dfs = [reps_df[r] for r in reps]
        fig = plt.figure(figsize=(7, 6))
        for groups in zip(*(f.groupby(['well-id', 'fov']) for f in dfs)):
            for i, (_, g) in enumerate(groups):
                sns.kdeplot(g[feature], lw=1, alpha=0.5, color=rep_colors[i % 7])
        for i, f in enumerate(dfs):
            sns.kdeplot(f[feature], color="black", lw=2, linestyle=styles[i % 7])
        handles = []
        for i, r in enumerate(reps):
            handles += [Line2D([0], [0], color=rep_colors[i % 7], lw=2,
                               label=f"{label} wells replicate {r}"),
                        Line2D([0], [0], color="black", lw=2, linestyle=styles[i % 7],
                               label=f"Global {label} replicate {r}")]
        ax = plt.gca()
        ax.legend(handles=handles, loc="best", fontsize=16)
        ax.set_title(cell_line, fontsize=22, fontname="Arial")
        _style(ax, xlabel)
        figs["overlay"] = fig

        fig = plt.figure(figsize=(7, 6))
        for i, f in enumerate(dfs):
            sns.kdeplot(f[feature], color="black", lw=2, linestyle=styles[i % 7])
        ax = plt.gca()
        ax.legend(handles=[Line2D([0], [0], color="black", lw=2, linestyle=styles[i % 7],
                                  label=f"Global {label} replicate {r}")
                           for i, r in enumerate(reps)], loc="best", fontsize=16)
        ax.set_title(cell_line, fontsize=22, fontname="Arial")
        _style(ax, xlabel)
        figs["global_only"] = fig
    return figs


def ridge_plot(cells: pd.DataFrame, cell_line: str, feature: str = FIG7B["feature"],
               treatments=FIG7B["treatments"], replicate: str = FIG7B["replicate"],
               xlabel: str = FIG7B["xlabel"]) -> Figure:
    """Figure 7b: joyplot of the treatments (figure7b_ridgeplot.py)."""
    import importlib
    import joypy
    _jj = importlib.import_module("joypy.joyplot")
    if not getattr(_jj, "_cpmetab_patched", False):   # joypy 0.2.6 + pandas >= 3
        _orig = _jj._flatten
        _jj._flatten = lambda axes: list(_orig(axes))
        _jj._cpmetab_patched = True
    from .plots import mpl_defaults
    mpl_defaults()
    df = cells[(cells["cell-line"] == cell_line) & (cells["replicate"] == str(replicate))]
    palette = sns.color_palette("Set2", n_colors=len(treatments))
    plot_data, colors = [], []
    treatments = [t for t in treatments if (df['treatment'].map(canon) == canon(t)).any()]
    if not treatments:
        raise ValueError(f"No cells for {cell_line} replicate {replicate}")
    palette = sns.color_palette("Set2", n_colors=len(treatments))
    for i, t in enumerate(treatments):
        f = df[df['treatment'].map(canon) == canon(t)].copy()
        f[feature] = _minmax(f[feature].astype(float))
        plot_data.append(pd.DataFrame({feature: f[feature], 'Treatment': t}))
        colors.append(palette[i])
    jdf = pd.concat(plot_data)
    jdf['Treatment'] = pd.Categorical(jdf['Treatment'], categories=list(treatments), ordered=True)
    fig, axes = plt.subplots(figsize=(6, 5))
    joypy.joyplot(data=jdf, by="Treatment", column=feature, overlap=0.6, fade=True, ax=axes,
                  linewidth=2, color=colors)
    plt.xticks([0, 0.5, 1], ['0', '0.5', '1'], fontsize=20)
    plt.xlabel(xlabel, fontsize=22, fontname="Arial")
    plt.title(f"{cell_line}", fontsize=24, fontname="Arial")
    plt.tight_layout()
    return plt.gcf()
