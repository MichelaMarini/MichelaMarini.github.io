"""Figure 6, as in original_scripts/figure6_step1..3.

Step 1: single cells of six treatments from the SPACe Step4_Features CSVs.
Step 2: per cell line, treatment and feature, each replicate is min-max
scaled; the Wasserstein distance of replicate 2 and 3 to replicate 1 is
signed by the median difference.
Step 3: features with value <= -0.15 (low) or 0.15 <= value <= 2 (up) are
counted per row of the step-2 table (all numeric columns of the row are
tested, including "Replicate").
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

from .config import canon
from .io import clean_treatment, parse_filename

#: Treatments analysed for Figure 6.
FIG6_TREATMENTS = ("DMSO", "Toyocamycin", "17-dmag", "Tms", "Artemisinin", "Fccp")
META_COLS = ["exp-id", "well-id", "fov", "treatment"]
FEATURE_TABLES = ("shape", "intensity", "texture")
#: Experiment folders: "<anything>_<CellLine>-<replicate>", e.g. "20250220-MariniPlates_5637-1".
FOLDER_RE = re.compile(r".*_([A-Za-z0-9]+)-(\d+)$")


@dataclass
class DistributionSettings:
    treatments: tuple[str, ...] = FIG6_TREATMENTS
    threshold: float = 0.15
    reference_replicate: str = "1"


@dataclass
class DistributionResults:
    settings: DistributionSettings
    distances: pd.DataFrame     #: cell-line, treatment, replicate, feature, signed_w
    counts: pd.DataFrame        #: per cell line, replicate, treatment: low/up counts
    n_features: int
    diagnostics: dict = field(default_factory=dict)


def feature_name(table: str, column: str) -> str:
    """Name of a single-cell feature column: the table name is added as a
    prefix only when the column does not already start with it, so SPACe's
    ``Shape_Area_cell`` stays ``Shape_Area_cell`` (not ``shape_Shape_Area_cell``)
    and ``mean_dapi`` from the intensity table becomes ``intensity_mean_dapi``."""
    column = str(column)
    return column if column.lower().startswith(table + "_") else f"{table}_{column}"


def clean_feature_name(name: str) -> str:
    """Remove the doubled prefix written by the original step 1 script
    (``shape_Shape_Area_cell`` -> ``Shape_Area_cell``)."""
    m = re.match(r"^(shape|intensity|texture)_(.*)$", str(name), flags=re.I)
    if m and m.group(2).lower().startswith(m.group(1).lower() + "_"):
        return m.group(2)
    return str(name)


def _clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    new = [clean_feature_name(c) for c in df.columns]
    if len(set(new)) == len(new):          # only if no two columns would get the same name
        df = df.set_axis(new, axis=1)
    return df


def _is_feature(df: pd.DataFrame, c) -> bool:
    return (str(c).split("_", 1)[0].lower() in FEATURE_TABLES
            and pd.api.types.is_numeric_dtype(df[c]))


def _find_csv(folder: Path, key: str) -> Path:
    hits = sorted(p for p in folder.glob("*.csv") if key.lower() in p.name.lower())
    if not hits:
        raise FileNotFoundError(f"No CSV containing '{key}' in {folder}")
    return hits[0]


def _n_rows(path: Path) -> int:
    with open(path, "rb") as f:
        return sum(1 for _ in f) - 1


def _read_rows(path: Path, rows: np.ndarray) -> pd.DataFrame:
    keep = set((rows + 1).tolist())  # +1: header line
    df = pd.read_csv(path, skiprows=lambda i: i > 0 and i not in keep)
    return df.loc[:, ~df.columns.str.contains(r"^Unnamed", case=False)]


def read_experiment(exp_dir: Path, treatments) -> pd.DataFrame | None:
    """Single cells of the requested treatments from one Step4_Features folder.

    The metadata, shape, intensity and texture tables are matched by row
    position (as SPACe writes them); their row counts must be equal.
    """
    step4 = exp_dir / "Step4_Features"
    if not step4.is_dir():
        return None
    try:
        cell_line, rep = parse_filename(exp_dir.name)
    except ValueError:
        return None
    paths = {"meta": _find_csv(step4, "metadata_features")}
    paths.update({k: _find_csv(step4, k) for k in FEATURE_TABLES})
    n = {k: _n_rows(p) for k, p in paths.items()}
    if len(set(n.values())) != 1:
        raise ValueError(f"{exp_dir.name}: tables have different row counts {n}; "
                         "rows cannot be matched by position")
    meta = pd.read_csv(paths["meta"])
    missing = [c for c in META_COLS if c not in meta.columns]
    if missing:
        raise KeyError(f"{exp_dir.name}: missing metadata columns {missing}")
    wanted = {canon(t) for t in treatments}
    mask = meta["treatment"].astype(str).str.strip().map(canon).isin(wanted).to_numpy()
    rows = np.flatnonzero(mask)
    if rows.size == 0:
        return None
    parts = [meta.loc[rows, META_COLS].reset_index(drop=True)]
    for k in FEATURE_TABLES:
        t = _read_rows(paths[k], rows)
        parts.append(t.set_axis([feature_name(k, c) for c in t.columns], axis=1)
                     .reset_index(drop=True))
    df = pd.concat(parts, axis=1)
    df["treatment"] = df["treatment"].map(clean_treatment)
    df.insert(0, "replicate", rep)
    df.insert(0, "cell-line", cell_line)
    return df


#: Files written by the original Figure 6 step 1: distribution_features_<CellLine>_<rep>.xlsx
EXTRACTED_RE = re.compile(r"_([A-Za-z0-9]+)[-_](\d+)$")


def read_extracted(path, treatments, name: str | None = None) -> pd.DataFrame | None:
    """One ``distribution_features_<CellLine>_<rep>.xlsx`` (or .csv) file,
    as written by the original step 1 script; ``path`` may be a file object
    when ``name`` is given."""
    name = name or Path(path).name
    try:
        cell_line, rep = parse_filename(name)
    except ValueError:
        return None
    df = pd.read_csv(path) if name.lower().endswith(".csv") else pd.read_excel(path)
    df = df.loc[:, ~df.columns.astype(str).str.contains(r"^Unnamed", case=False)]
    df = _clean_columns(df)
    if treatments is None:                       # keep every treatment in the file
        treatments = df["treatment"].astype(str).unique()
    missing = [c for c in META_COLS if c not in df.columns]
    if missing:
        raise KeyError(f"{path.name}: missing metadata columns {missing}")
    wanted = {canon(t) for t in treatments}
    df = df[df["treatment"].astype(str).str.strip().map(canon).isin(wanted)].copy()
    if df.empty:
        return None
    df["treatment"] = df["treatment"].map(clean_treatment)
    df.insert(0, "replicate", rep)
    df.insert(0, "cell-line", cell_line)
    return df


def load_single_cell_files(items, treatments=None) -> tuple[pd.DataFrame, list[str]]:
    """Single cells from uploaded distribution_features files: (name, file object) pairs."""
    frames = [df for name, obj in sorted(items, key=lambda x: x[0])
              if (df := read_extracted(obj, treatments, name)) is not None]
    if not frames:
        raise ValueError("No usable distribution_features_<CellLine>_<replicate> files")
    cells = pd.concat(frames, ignore_index=True)
    features = [c for c in cells.columns if _is_feature(cells, c)]
    return cells, features


def load_single_cell(root: str | Path, treatments=FIG6_TREATMENTS) -> tuple[pd.DataFrame, list[str]]:
    """Load single cells from ``root``. Returns (cells, feature names).

    ``root`` can be either the folder of SPACe experiment folders
    (``<plate>_<CellLine>-<rep>/Step4_Features``), or the "Distribution
    Features" folder of ``distribution_features_<CellLine>_<rep>.xlsx`` files
    written by the original Figure 6 step 1 script. The Excel files only
    contain the treatments that script kept.
    """
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"Folder not found: {root}")
    frames = [df for p in sorted(root.iterdir()) if p.is_dir()
              if (df := read_experiment(p, treatments)) is not None]
    if not frames:
        xlsx = sorted(p for p in root.glob("*distribution_features_*.xlsx")
                      if not p.name.startswith("~$"))
        frames = [df for p in xlsx if (df := read_extracted(p, treatments)) is not None]
    if not frames:
        raise FileNotFoundError(
            f"Nothing to read in {root}: expected '<plate>_<CellLine>-<rep>/Step4_Features' "
            "folders or 'distribution_features_<CellLine>_<rep>.xlsx' files")
    cells = pd.concat(frames, ignore_index=True)
    features = [c for c in cells.columns if _is_feature(cells, c)]
    return cells, features


def _minmax(x: np.ndarray) -> np.ndarray:
    """sklearn MinMaxScaler: NaN ignored, zero range -> scale 1."""
    lo, hi = np.nanmin(x), np.nanmax(x)
    rng = hi - lo
    return (x - lo) / (rng if rng != 0 else 1.0)


def signed_wasserstein(cells: pd.DataFrame, features: list[str],
                       s: DistributionSettings) -> pd.DataFrame:
    rows = []
    ref = s.reference_replicate
    for (cl, tr), g in cells.groupby(["cell-line", "treatment"], sort=True):
        reps = sorted(g["replicate"].unique(), key=lambda r: int(r) if r.isdigit() else r)
        if ref not in reps:
            continue
        others = [r for r in reps if r != ref]
        data = {r: g.loc[g["replicate"] == r, features].to_numpy(dtype=float) for r in reps}
        for j, f in enumerate(features):
            cols = {r: data[r][:, j] for r in reps}
            scaled = {r: _minmax(v) for r, v in cols.items()}   # each replicate alone
            for r in others:
                if not len(scaled[ref]) or not len(scaled[r]):
                    w = np.nan
                else:
                    a, b = scaled[ref], scaled[r]
                    w = np.sign(np.median(b) - np.median(a)) * wasserstein_distance(a, b)
                rows.append((cl, tr, r, f, w))
    return pd.DataFrame(rows, columns=["cell-line", "treatment", "replicate", "feature", "signed_w"])


def count_outliers(dist: pd.DataFrame, threshold: float = 0.15) -> pd.DataFrame:
    d = dist.assign(low=dist["signed_w"] <= -threshold,
                    up=(dist["signed_w"] >= threshold) & (dist["signed_w"] <= 2))
    g = d.groupby(["cell-line", "replicate", "treatment"], sort=False)
    out = g.agg(low_count=("low", "sum"), up_count=("up", "sum"),
                n_nan=("signed_w", lambda x: int(x.isna().sum()))).reset_index()
    out["low_features"] = g.apply(lambda x: ", ".join(x.loc[x["low"], "feature"]),
                                  include_groups=False).to_numpy()
    out["up_features"] = g.apply(lambda x: ", ".join(x.loc[x["up"], "feature"]),
                                 include_groups=False).to_numpy()
    # step 3 tests every numeric column of the row, including "Replicate"
    extra = out["replicate"].astype(int).between(threshold, 2)
    out.loc[extra, "up_count"] += 1
    out.loc[extra, "up_features"] = ("Replicate, " + out.loc[extra, "up_features"]
                                     ).str.rstrip(", ")
    return out


def run_distribution_analysis(root: str | Path, settings: DistributionSettings | None = None,
                              cells: pd.DataFrame | None = None,
                              features: list[str] | None = None) -> DistributionResults:
    s = settings or DistributionSettings()
    if cells is None:
        cells, features = load_single_cell(root, s.treatments)
    dist = signed_wasserstein(cells, features, s)
    counts = count_outliers(dist, s.threshold)
    found = set(cells["treatment"].map(canon))
    diag = {"n_cells": len(cells),
            "n_features": len(features),
            "cells_per_group": cells.groupby(["cell-line", "replicate", "treatment"]).size()
            .unstack("treatment").to_dict(orient="index"),
            "treatments_not_found": [t for t in s.treatments if canon(t) not in found],
            "n_nan_distances": int(dist["signed_w"].isna().sum())}
    diag["cells_per_group"] = {f"{k[0]} rep {k[1]}": v for k, v in diag["cells_per_group"].items()}
    return DistributionResults(s, dist, counts, len(features), diag)


# --------------------------------------------------------------------------
# Tables: reading back (csv/xlsx, this package's or the original format)
LONG_COLS = ["cell-line", "treatment", "replicate", "feature", "signed_w"]


def _from_frame(df: pd.DataFrame) -> pd.DataFrame | None:
    cols = {c.lower(): c for c in df.columns}
    if set(LONG_COLS) <= set(cols):
        out = df[[cols[c] for c in LONG_COLS]].copy()
        out.columns = LONG_COLS
    else:
        wide = [c for c in df.columns if str(c).endswith("_Wess_Distance")]
        if not wide or not {"cell-line", "replicate", "treatment"} <= set(cols):
            return None
        out = df.melt(id_vars=[cols["cell-line"], cols["treatment"], cols["replicate"]],
                      value_vars=wide, var_name="feature", value_name="signed_w")
        out.columns = LONG_COLS
        out["feature"] = out["feature"].str.removesuffix("_Wess_Distance")
    from .io import canonical_cell_line
    out["cell-line"] = out["cell-line"].astype(str).map(canonical_cell_line)
    out["feature"] = out["feature"].map(clean_feature_name)
    out["replicate"] = out["replicate"].astype(str).str.replace(r"\.0$", "", regex=True)
    out["signed_w"] = pd.to_numeric(out["signed_w"], errors="coerce")
    return out


def read_wasserstein_table(source, name: str) -> pd.DataFrame:
    """Signed Wasserstein distances from a .csv or .xlsx file.

    Accepts the long table written by this package (``wasserstein_distances.csv``
    / ``.xlsx``) and the wide tables of the original scripts
    (``wasserstein_distances.xlsx``, ``wasserstein_distances_with_counts.xlsx``).
    """
    if str(name).lower().endswith(".csv"):
        frames = [pd.read_csv(source)]
    else:
        frames = list(pd.read_excel(source, sheet_name=None).values())
    for df in frames:
        out = _from_frame(df)
        if out is not None:
            return out
    raise ValueError(f"`{name}` does not contain signed Wasserstein distances "
                     "(expected columns ending in `_Wess_Distance`, or "
                     "cell-line/treatment/replicate/feature/signed_w)")


def to_original_layout(dist: pd.DataFrame, counts: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """The two sheets of the original ``wasserstein_distances_with_counts.xlsx``."""
    wide = dist.pivot_table(index=["cell-line", "replicate", "treatment"], columns="feature",
                            values="signed_w", aggfunc="first", dropna=False)
    feats = list(dict.fromkeys(dist["feature"]))
    wide = wide.reindex(columns=feats)
    wide.columns = [f"{f}_Wess_Distance" for f in feats]
    wide = wide.reset_index()
    c = counts.rename(columns={"low_count": "Low_Count", "up_count": "Up_Count",
                               "low_features": "Low_Features", "up_features": "Up_Features"})
    keys = ["cell-line", "replicate", "treatment"]
    sheet1 = wide.merge(c[keys + ["Low_Count", "Up_Count", "Low_Features", "Up_Features"]],
                        on=keys, how="left")
    ren = {"cell-line": "Cell-line", "replicate": "Replicate", "treatment": "Treatment"}
    sheet1 = sheet1.rename(columns=ren)
    sheet2 = sheet1[["Cell-line", "Replicate", "Treatment", "Low_Count", "Up_Count",
                     "Low_Features", "Up_Features"]]
    return {"Sheet1": sheet1, "Sheet2": sheet2}


# --------------------------------------------------------------------------
# Per-feature summary: which features deviate between biological replicates
def feature_reproducibility(dist: pd.DataFrame, threshold: float = 0.15,
                            max_pct_deviating: float = 0.0) -> pd.DataFrame:
    """One row per feature: in how many replicate comparisons (cell line x
    treatment x replicate 2/3 vs 1) its signed Wasserstein distance is
    <= -threshold or between threshold and 2 (the Figure 6 rule).

    ``status`` is "consistent" when the percentage of deviating comparisons
    is <= ``max_pct_deviating`` (default 0: never deviates), otherwise
    "deviating".
    """
    w = dist["signed_w"]
    d = dist.assign(low=w <= -threshold, up=(w >= threshold) & (w <= 2),
                    valid=w.notna())
    d["dev"] = d["low"] | d["up"]
    g = d.groupby("feature", sort=False)
    out = pd.DataFrame({
        "n_comparisons": g["valid"].sum(),
        "n_low": g["low"].sum(),
        "n_up": g["up"].sum(),
        "n_deviating": g["dev"].sum(),
    })
    out["pct_deviating"] = (100 * out["n_deviating"] / out["n_comparisons"]).round(2)
    out["max_abs_w"] = g["signed_w"].apply(lambda x: x.abs().max())
    per_cl = d[d["dev"]].groupby(["feature", "cell-line"]).size().unstack(fill_value=0)
    out["cell_lines_affected"] = per_cl.gt(0).sum(axis=1).reindex(out.index).fillna(0).astype(int)
    out["deviating_in"] = d[d["dev"]].groupby("feature").apply(
        lambda x: "; ".join(f"{r['cell-line']} {r['treatment']} rep{r['replicate']}"
                            for _, r in x.iterrows()), include_groups=False
    ).reindex(out.index).fillna("")
    out["status"] = np.where(out["pct_deviating"] <= max_pct_deviating, "consistent", "deviating")
    out = out.reset_index().sort_values(["status", "pct_deviating", "feature"])
    per_cl = per_cl.add_prefix("n_deviating_").reset_index()
    out = out.merge(per_cl, on="feature", how="left").fillna(
        {c: 0 for c in per_cl.columns if c != "feature"})
    return out.reset_index(drop=True)
