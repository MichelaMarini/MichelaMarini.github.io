"""Reading SPACe well-wise distance-map spreadsheets."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable

import pandas as pd

from .config import DMSO, EXCLUDED_TREATMENTS, FEATURE_START_COL


@dataclass
class Dataset:
    wells: pd.DataFrame          #: one row per well: metadata + features
    features: list[str]          #: feature column names
    file_order: list[tuple[str, str, str]]  #: (file name, cell line, replicate)

    @property
    def replicates(self) -> list[str]:
        return sorted(self.wells["replicate"].unique(), key=_rep_key)

    @property
    def cell_lines(self) -> list[str]:
        return sorted(self.wells["cell-line"].unique())


def _rep_key(r: str):
    return (0, int(r)) if str(r).isdigit() else (1, str(r))


def canonical_cell_line(name: str) -> str:
    """Map spelling variants to the names used in the figures
    (u2os -> U2OS, PANC1 -> Panc1, MDA-MB-231 -> MDAMB231, ...)."""
    from .config import CELL_LINE_ORDER, canon
    key = canon(name)
    for known in CELL_LINE_ORDER:
        if key == canon(known):
            return known
    for known in CELL_LINE_ORDER:            # e.g. "MDAMB231-new"
        if key.startswith(canon(known)):
            return known
    return name


def parse_filename(name: str) -> tuple[str, str]:
    """File or folder name -> (cell line, replicate).

    Works for every cell line and for ``..._U2OS_1``, ``..._U2OS-1``,
    ``..._MDA-MB-231_1``, ``..._Panc1-2 (1)`` (Windows copy), etc.
    """
    stem = Path(name).stem if Path(name).suffix.lower() in (".xlsx", ".csv") else Path(name).name
    stem = re.sub(r"\s*\(\d+\)$", "", stem.strip())        # " (1)" added by Windows
    m = re.search(r"_([A-Za-z0-9][A-Za-z0-9-]*?)[-_](\d+)$", stem)
    if not m:
        raise ValueError(f"Cannot read cell line and replicate from `{Path(name).name}`: "
                         "expected a name ending in `_CellLine_1` or `_CellLine-1`")
    return canonical_cell_line(m.group(1)), m.group(2)


def clean_treatment(name) -> str:
    s = str(name).strip()
    return DMSO if s.lower() == "dmso" else s.lower().capitalize()


def read_plate(source: str | Path | BinaryIO, name: str | None = None,
               feature_start: int = FEATURE_START_COL) -> tuple[pd.DataFrame, list[str]]:
    """Read one distance-map spreadsheet (one cell line, one replicate)."""
    name = name or str(source)
    cell_line, replicate = parse_filename(name)
    df = pd.read_excel(source)
    for col in ("treatment", "cell-count"):
        if col not in df.columns[:feature_start]:
            raise ValueError(f"`{name}`: column `{col}` not found in the first "
                             f"{feature_start} (metadata) columns")
    features = [str(c) for c in df.columns[feature_start:]]
    df.columns = [str(c) for c in df.columns]
    non_numeric = [c for c in features if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        raise ValueError(f"'{name}': non-numeric feature columns {non_numeric[:5]}; "
                         f"check FEATURE_START_COL ({feature_start})")
    df["treatment_raw"] = df["treatment"]
    df["treatment"] = df["treatment"].map(clean_treatment)
    df = df[~df["treatment_raw"].astype(str).str.lower().isin(EXCLUDED_TREATMENTS)].copy()
    df["cell-line"] = cell_line
    df["replicate"] = replicate
    df["source_file"] = Path(name).name
    return df, features


def load_dataset(source: str | Path | Iterable[tuple[str, BinaryIO]],
                 feature_start: int = FEATURE_START_COL) -> Dataset:
    """Load every ``*.xlsx`` in a folder, or an iterable of (file name, file object)."""
    if isinstance(source, (str, Path)):
        paths = sorted(p for p in Path(source).glob("*.xlsx") if not p.name.startswith("~$"))
        if not paths:
            raise FileNotFoundError(f"No .xlsx files found in {source}")
        items = [(p.name, p) for p in paths]
    else:
        items = sorted(source, key=lambda x: x[0])
    frames, features, order = [], None, []
    for name, obj in items:
        df, feats = read_plate(obj, name, feature_start)
        if features is None:
            features = feats
        elif feats != features:
            raise ValueError(f"'{name}' has different feature columns from the first file")
        frames.append(df)
        order.append((Path(name).name, df["cell-line"].iat[0] if len(df) else "",
                      df["replicate"].iat[0] if len(df) else ""))
    wells = pd.concat(frames, ignore_index=True)
    dup = wells.groupby(["cell-line", "replicate"])["source_file"].nunique()
    if (dup > 1).any():
        raise ValueError(f"More than one file for the same cell line and replicate: "
                         f"{dup[dup > 1].index.tolist()}")
    return Dataset(wells=wells, features=features, file_order=order)
