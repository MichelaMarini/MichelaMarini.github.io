"""Replicate Pearson correlations, as in cell_count_version_4.py."""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def replicate_correlations(profiles: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """|Pearson r| between replicate pairs, per cell line and treatment,
    computed across ``columns`` (the original: all EMD features plus the
    ``euclidean_distance`` column). Only cell lines and treatments present in
    every replicate are used."""
    reps = sorted(profiles["replicate"].unique(), key=lambda r: (not str(r).isdigit(),
                  int(r) if str(r).isdigit() else r))
    by_rep = {r: profiles[profiles["replicate"] == r].set_index(["cell-line", "treatment"])
              for r in reps}
    common = set.intersection(*(set(df.index) for df in by_rep.values()))
    rows = []
    if len(reps) < 2:
        return pd.DataFrame(columns=["cell-line", "treatment", "pair", "r"])
    for cl, tr in sorted(common):
        v = {r: by_rep[r].loc[(cl, tr), columns].to_numpy(dtype=float) for r in reps}
        for a, b in combinations(reps, 2):
            with np.errstate(all="ignore"):
                r = pearsonr(v[a], v[b]).statistic
            rows.append({"cell-line": cl, "treatment": tr, "pair": f"{a}-{b}", "r": abs(r)})
    return pd.DataFrame(rows)


def correlation_summary(corr: pd.DataFrame, cell_line_order) -> pd.DataFrame:
    """Median and SD (np.std) of the correlations per cell line, and overall."""
    g = corr.groupby("cell-line")["r"]
    out = pd.DataFrame({"n": g.count(), "median_r": g.median(),
                        "sd_r": g.agg(lambda x: np.std(x.dropna()))})
    out = out.reindex([c for c in cell_line_order if c in out.index])
    out.loc["All"] = [corr["r"].count(), corr["r"].median(), np.std(corr["r"].dropna())]
    out["n"] = out["n"].astype("Int64")
    return out
