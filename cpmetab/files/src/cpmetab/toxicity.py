"""Cell-count heatmaps, as in cell_count_version_4.py."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import DMSO

KEYS = ["replicate", "cell-line"]


def mean_cell_counts(wells: pd.DataFrame, min_dmso_cells: int = 100) -> pd.DataFrame:
    """Mean cell count per treatment (rounded), after dropping DMSO wells < min_dmso_cells."""
    bad = (wells["treatment"] == DMSO) & (wells["cell-count"] < min_dmso_cells)
    out = wells.loc[~bad].groupby(KEYS + ["treatment"], as_index=False)["cell-count"].mean()
    out["cell-count"] = out["cell-count"].round(0).astype(int)
    return out


def percent_of_dmso(counts: pd.DataFrame) -> pd.DataFrame:
    """Cell count as % of DMSO of the same cell line and replicate, capped at 100."""
    dmso = counts.loc[counts["treatment"] == DMSO].set_index(KEYS)["cell-count"]
    missing = counts.set_index(KEYS).index.unique().difference(dmso.index)
    if len(missing):
        raise ValueError(f"No DMSO row found for {list(missing)}")
    ref = dmso.reindex(pd.MultiIndex.from_frame(counts[KEYS])).to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = counts["cell-count"].to_numpy(dtype=float) / ref * 100.0
    out = counts.copy()
    out["pct_dmso"] = np.minimum(pct, 100.0)
    out.loc[out["treatment"] == DMSO, "pct_dmso"] = 100.0
    return out


def flagged_treatments(pct: pd.DataFrame, threshold: float = 50.0) -> list[str]:
    hit = pct.loc[(pct["pct_dmso"] < threshold) & (pct["treatment"] != DMSO), "treatment"]
    return sorted(hit.unique())


def toxicity_matrices(pct: pd.DataFrame, threshold: float, cell_line_order) -> dict[str, pd.DataFrame]:
    """Per replicate: cell lines x (DMSO + treatments below threshold anywhere).

    Missing values are 0; columns ordered by the summed % of the first
    replicate (descending), DMSO first.
    """
    flagged = flagged_treatments(pct, threshold)
    reps = sorted(pct["replicate"].unique(), key=lambda r: (not str(r).isdigit(),
                  int(r) if str(r).isdigit() else r))
    mats = {}
    for rep in reps:
        sub = pct[(pct["replicate"] == rep) & pct["treatment"].isin(flagged + [DMSO])]
        m = sub.pivot(index="cell-line", columns="treatment", values="pct_dmso")
        mats[rep] = m.reindex(index=list(cell_line_order), columns=[DMSO] + flagged).fillna(0.0)
    first = mats[reps[0]][flagged].sum(axis=0)
    order = sorted(flagged, key=lambda t: (-first[t], t))
    return {r: m[[DMSO] + order] for r, m in mats.items()}
