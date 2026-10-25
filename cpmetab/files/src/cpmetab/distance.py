"""Euclidean distance (ED) from DMSO, as in cell_count_version_4.py and
figure3b_control_ed.py."""
from __future__ import annotations

import numpy as np
import pandas as pd
from statistics import NormalDist

from .config import DMSO, canon

KEYS = ["replicate", "cell-line"]
_MAD_C = NormalDist().inv_cdf(0.75)   # statsmodels.robust.mad constant


def dmso_reference(wells: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per replicate and cell line: DMSO median and MAD (statsmodels.robust.mad)."""
    d = wells[wells["treatment"] == DMSO]
    med = d.groupby(KEYS)[features].median()
    mad = d.groupby(KEYS)[features].agg(
        lambda x: np.median(np.abs(x.to_numpy() - np.median(x.to_numpy()))) / _MAD_C)
    return med, mad


def treatment_means(wells: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    """Mean feature values per replicate, cell line and treatment."""
    return wells.groupby(KEYS + ["treatment"], as_index=False)[features].mean()


def euclidean_distance(means: pd.DataFrame, med: pd.DataFrame, mad: pd.DataFrame,
                       features: list[str], reference_cell_line: dict[str, str],
                       tol: float) -> pd.Series:
    """ED as computed by the original script.

    A row is z-normalised with the DMSO median/MAD unless a DMSO MAD of its
    cell line is below ``tol``, in which case raw values are used. ED is the
    Euclidean distance to the DMSO medians of ``reference_cell_line`` (the
    first file listed, for every row).
    """
    zero = (mad[features].abs() < tol).any(axis=1)
    if not zero.all():
        raise RuntimeError(
            "The original script only runs when every cell line/replicate has a feature "
            "with zero DMSO MAD; these data do not satisfy that.")
    out = pd.Series(np.nan, index=means.index)
    for rep, ref_cl in reference_cell_line.items():
        rows = means["replicate"] == rep
        ref = med.loc[(rep, ref_cl), features].to_numpy(dtype=float)
        diff = means.loc[rows, features].to_numpy(dtype=float) - ref
        out[rows] = np.sqrt((diff ** 2).sum(axis=1))
    return out


def responsive_treatments(ed: pd.DataFrame, threshold: float) -> list[str]:
    """Treatments with ED >= threshold in at least one cell line and replicate, plus DMSO."""
    return sorted(set(ed.loc[ed["ed"] >= threshold, "treatment"]) | {DMSO})


def match(available, wanted) -> tuple[dict[str, str], list[str]]:
    """Map data names to the spelling in ``wanted`` (case/hyphen-insensitive)."""
    lookup = {canon(w): w for w in wanted}
    mapping = {a: lookup[canon(a)] for a in available if canon(a) in lookup}
    missing = [w for w in wanted if canon(w) not in {canon(a) for a in available}]
    return mapping, missing


def global_order(ed: pd.DataFrame, compounds, cell_line_order) -> list[str]:
    """``compute_global_treatment_order``: DMSO first, then by mean ED (descending)."""
    mapping, _ = match(sorted(ed["treatment"].unique()), compounds)
    d = ed[ed["treatment"].isin(mapping) & ed["cell-line"].isin(cell_line_order)].copy()
    d["treatment"] = d["treatment"].map(mapping)
    agg = d.pivot_table(index="treatment", columns="cell-line", values="ed", aggfunc="mean")
    score = agg.reindex(columns=list(cell_line_order)).mean(axis=1)
    top = [t for t in score.index if canon(t) == canon(DMSO)]
    rest = sorted([t for t in score.index if canon(t) != canon(DMSO)], key=lambda t: -score[t])
    return top + rest


def ed_matrix(ed: pd.DataFrame, replicate: str, rows, columns, rename=None) -> pd.DataFrame:
    """Treatments x cell lines for one replicate; non-finite values set to 0."""
    d = ed[ed["replicate"] == replicate].copy()
    if rename:
        d = d[d["treatment"].isin(rename)]
        d["treatment"] = d["treatment"].map(rename)
    m = d.pivot(index="treatment", columns="cell-line", values="ed")
    m = m.reindex(index=list(rows), columns=list(columns))
    return m.where(np.isfinite(m), 0.0)


def control_features(features: list[str], feature_start: int, control_start: int) -> list[str]:
    """``df.columns[16:-2]`` after one appended column = columns 16 .. second-to-last."""
    return features[control_start - feature_start:-1]


def control_ed(wells: pd.DataFrame, features: list[str]) -> pd.Series:
    """Figure 3B: norm of the EMD vector of each well."""
    x = wells[features].to_numpy(dtype=float)
    return pd.Series(np.sqrt((x ** 2).sum(axis=1)), index=wells.index)
