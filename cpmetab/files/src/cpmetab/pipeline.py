"""Runs the analysis of cell_count_version_4.py and figure3b_control_ed.py."""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from . import distance as D
from . import reproducibility as R
from . import toxicity as T
from .config import (CONTROL_FEATURE_START_COL, DMSO, FEATURE_START_COL,
                     HEATMAP1_EXCLUDE, HEATMAP1_LAST, HEATMAP2_CELL_ORDER, Settings, canon,
                     ordered)
from .io import Dataset


@dataclass
class Results:
    settings: Settings
    cell_lines: list[str]                     #: cell lines in the data, in figure order
    replicates: list[str]                     #: biological replicates in the data
    cell_counts: pd.DataFrame                 #: % of DMSO per replicate/cell line/treatment
    toxicity: dict[str, pd.DataFrame]         #: per replicate: cell lines x treatments
    ed: pd.DataFrame                          #: ED per replicate/cell line/treatment (all)
    ed_selected: pd.DataFrame                 #: rows of treatments with ED >= threshold
    heatmaps: dict[str, dict[str, pd.DataFrame]]  #: variant -> replicate -> matrix
    correlations: pd.DataFrame
    correlation_summary: pd.DataFrame
    control_ed: pd.DataFrame                  #: Figure 3B, per well
    diagnostics: dict = field(default_factory=dict)


def _heatmap1(ed_sel: pd.DataFrame, rep: str, cell_order) -> pd.DataFrame:
    """``create_heatmap``: exclusion list, custom order, all responsive treatments."""
    d = ed_sel[(ed_sel["replicate"] == rep) & ~ed_sel["treatment"].isin(HEATMAP1_EXCLUDE)]
    present = set(d["treatment"])
    order = ([DMSO] + sorted(present - {DMSO, *HEATMAP1_LAST}) + list(HEATMAP1_LAST))
    return D.ed_matrix(d, rep, [t for t in order if t in present], cell_order)


def _heatmap2(ed_sel: pd.DataFrame, rep: str, compounds, cell_lines) -> pd.DataFrame:
    """``create_heatmap_2``: compound list in list order, other cell-line order."""
    rename, _ = D.match(sorted(ed_sel["treatment"].unique()), compounds)
    rows = [c for c in compounds if c in set(rename.values())]
    return D.ed_matrix(ed_sel, rep, rows, ordered(cell_lines, HEATMAP2_CELL_ORDER), rename)


def run_analysis(data: Dataset, settings: Settings | None = None) -> Results:
    s = settings or Settings()
    wells, features = data.wells, data.features
    lines = [c for c in (s.cell_line_order or ordered(data.cell_lines)) if c in data.cell_lines]
    diag = {"n_files": len(data.file_order), "n_features": len(features),
            "replicates": data.replicates, "cell_lines": data.cell_lines}

    # --- cell counts -------------------------------------------------------
    pct = T.percent_of_dmso(T.mean_cell_counts(wells, s.min_dmso_cells))
    tox = T.toxicity_matrices(pct, s.cell_loss_threshold, lines)
    diag["treatments_below_threshold"] = T.flagged_treatments(pct, s.cell_loss_threshold)

    # --- Euclidean distance ------------------------------------------------
    med, mad = D.dmso_reference(wells, features)
    means = D.treatment_means(wells, features)
    ref = {}
    for rep in data.replicates:
        first = [cl for _, cl, r in data.file_order if r == rep][0]
        given = s.reference_cell_line
        ref[rep] = given.get(rep, first) if isinstance(given, dict) else (given or first)
    diag["reference_cell_line"] = ref
    means["ed"] = D.euclidean_distance(means, med, mad, features, ref, s.mad_tolerance)
    ed = means[["replicate", "cell-line", "treatment", "ed"]].copy()
    keep = D.responsive_treatments(ed, s.ed_threshold)
    ed_sel = ed[ed["treatment"].isin(keep)]
    _, missing = D.match(sorted(ed["treatment"].unique()), s.compounds)
    diag["compounds_not_in_data"] = missing
    order = D.global_order(ed_sel, s.compounds, lines)
    rename, _ = D.match(sorted(ed_sel["treatment"].unique()), s.compounds)
    heatmaps = {
        "heatmap1": {r: _heatmap1(ed_sel, r, lines) for r in data.replicates},
        "heatmap2": {r: _heatmap2(ed_sel, r, s.compounds, lines) for r in data.replicates},
        "heatmap_ordered": {r: D.ed_matrix(ed_sel, r, order, lines, rename)
                            for r in data.replicates},
    }

    # --- reproducibility (features + ED column, |r|) -------------------------
    profiles = means.rename(columns={"ed": "euclidean_distance"})
    corr = R.replicate_correlations(profiles, features + ["euclidean_distance"])
    summary = R.correlation_summary(corr, lines)

    # --- Figure 3B: control wells --------------------------------------------
    cfeat = D.control_features(features, FEATURE_START_COL, CONTROL_FEATURE_START_COL)
    label = {canon(c): c for c in s.controls}
    cw = wells[wells["treatment"].map(canon).isin(label)]
    cw = cw[~((cw["treatment"] == DMSO) & (cw["cell-count"] < s.min_dmso_cells))].copy()
    cw["ed"] = D.control_ed(cw, cfeat)
    cw["control"] = cw["treatment"].map(lambda t: label[canon(t)])
    cols = ["replicate", "cell-line", "control"] + (["well-id"] if "well-id" in cw else []) + ["ed"]
    control = cw[cols].reset_index(drop=True)
    diag["controls_not_in_data"] = [c for c in s.controls if c not in set(control["control"])]
    diag["control_features"] = f"{cfeat[0]} .. {cfeat[-1]} ({len(cfeat)})" if cfeat else "none"

    return Results(s, lines, data.replicates, pct, tox, ed, ed_sel, heatmaps, corr, summary, control, diag)
