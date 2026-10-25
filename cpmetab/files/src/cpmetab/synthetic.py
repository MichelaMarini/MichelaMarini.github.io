"""Generate small synthetic SPACe-like distance-map files (for tests and the
web-app demo). Values are random and carry no biological meaning.

Real input: one Excel file per cell line and biological replicate, named
``<prefix>_<CellLine>_<replicate>.xlsx``, one row per well, nine metadata
columns (incl. ``treatment``, ``cell-line`` and ``cell-count``) followed by
one directional-EMD column per feature.
"""
from pathlib import Path

import numpy as np
import pandas as pd

CELL_LINES = ["HepG2", "MDAMB231", "PC3", "5637", "U2OS", "Panc1", "A549"]
COMPOUNDS = [
    "mefloquine", "vps34-in1", "ym-201636", "fccp", "rottlerin", "toyocamycin", "17-dmag",
    "6-aminonicotinamide", "brequinar", "erastin", "tms", "gsk1059615", "torin-1", "tmp-153",
    "6-thioguanine", "gemcitabine", "adaptaquin", "dorsomorphin", "hqno", "nk-252",
    "artemisinin", "methotrexate", "olaparib", "metformin", "oligomycin", "2-deoxyglucose",
    "auranofin", "plumbagin", "rotenone", "digitoxin", "berb", "actd", "mk2206",
]
TOXIC = {"auranofin": 0.15, "plumbagin": 0.3, "digitoxin": 0.35, "actd": 0.1, "berb": 0.4}
META = ["plate", "well-id", "treatment", "dosage", "density", "time", "cell-line",
        "other", "cell-count"]
N_FEATURES = 25


def make_frames(seed=0):
    """Return {file name: DataFrame} for 7 cell lines x 3 replicates."""
    rng = np.random.default_rng(seed)
    frames = {}
    effects = {c: rng.normal(0, rng.uniform(0.05, 0.6), N_FEATURES) for c in COMPOUNDS}
    for ci, cl in enumerate(CELL_LINES):
        base_count = rng.integers(900, 2000)
        line_tox = rng.uniform(0.6, 1.2)
        for rep in (1, 2, 3):
            rows = []
            wells = ["dmso"] * 28 + [c for c in COMPOUNDS for _ in range(2)] + ["empty"] * 4
            for wi, t in enumerate(wells):
                if t == "empty":
                    count, feats = 0, np.zeros(N_FEATURES)
                elif t == "dmso":
                    count = int(base_count * rng.uniform(0.9, 1.1))
                    feats = rng.normal(0, 0.05, N_FEATURES)
                    if wi in (3, 17):              # two failed DMSO wells (<100 cells)
                        count = int(rng.integers(20, 90))
                    if True:  # one constant feature per file (as in the real data)
                        feats[0] = 0.0
                else:
                    frac = TOXIC.get(t, 1.0)
                    frac = min(1.0, frac * line_tox) if t in TOXIC else rng.uniform(0.7, 1.05)
                    if t == "rotenone" and cl in ("U2OS", "A549"):
                        frac = 0.45
                    count = int(base_count * frac * rng.uniform(0.9, 1.1))
                    feats = effects[t] * rng.uniform(0.7, 1.3) + rng.normal(0, 0.05, N_FEATURES)
                rows.append([f"P{rep}", f"W{wi:03d}", t, 10, 3000, 24, cl, "", count, *feats])
            df = pd.DataFrame(rows, columns=META + [f"feature_{i:03d}" for i in range(N_FEATURES)])
            frames[f"distmap_{cl}_{rep}.xlsx"] = df
    return frames


def make(out_dir, seed=0):
    """Write the synthetic files to ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in make_frames(seed).items():
        df.to_excel(out / name, index=False)


if __name__ == "__main__":
    import sys
    make(sys.argv[1] if len(sys.argv) > 1 else "synthetic_data")


# --------------------------------------------------------------------------
# Single-cell features (SPACe Step4_Features), for the Figure 6 analysis
SC_TREATMENTS = {"dmso": 4, "toyocamycin": 2, "17-dmag": 2, "tms": 2, "artemisinin": 2,
                 "fccp": 2, "olaparib": 2}
SC_FEATURES = {
    "shape": ["Shape_Area_cell", "Shape_Area_cyto", "Shape_Cell_perimeter"],
    "intensity": ["mean_dapi", "mean_mito", "mean_rna"],
    "texture": ["contrast", "homogeneity", "flag"],
}


def make_single_cell(root, seed=0, cells_per_well=(30, 60)):
    """Write ``root/<plate>_<CellLine>-<rep>/Step4_Features/*.csv``."""
    rng = np.random.default_rng(seed)
    root = Path(root)
    effects = {(cl, t): 0.0 if t == "dmso" else rng.normal(0, 0.4)
               for cl in CELL_LINES for t in SC_TREATMENTS}   # same in every replicate
    for cl in CELL_LINES:
        for rep in (1, 2, 3):
            shift = {1: 0.0, 2: 0.08, 3: -0.05}[rep]       # replicate batch effect
            scale = {1: 1.0, 2: 1.0, 3: 1.4}[rep]
            meta, tabs = [], {k: [] for k in SC_FEATURES}
            w = 0
            for t, n_wells in SC_TREATMENTS.items():
                effect = effects[(cl, t)]
                for _ in range(n_wells):
                    w += 1
                    n = int(rng.integers(*cells_per_well))
                    for i in range(n):
                        meta.append([f"EXP{rep}", f"W{w:03d}", int(rng.integers(1, 4)), t, i])
                    for kind, names in SC_FEATURES.items():
                        for _ in range(n):
                            row = []
                            for name in names:
                                if name == "flag":
                                    row.append(0)
                                elif name.startswith("Shape_Area"):
                                    row.append(int(rng.normal(500 * (1 + effect + shift), 80 * scale)))
                                else:
                                    row.append(rng.normal(1 + effect + shift, 0.2 * scale))
                            tabs[kind].append(row)
            d = root / f"20250220-Plates_{cl}-{rep}" / "Step4_Features"
            d.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(meta, columns=["exp-id", "well-id", "fov", "treatment", "cell-id"]
                         ).to_csv(d / "metadata_features.csv")       # with index column
            for kind, names in SC_FEATURES.items():
                pd.DataFrame(tabs[kind], columns=names).to_csv(d / f"{kind}_features.csv",
                                                               index=False)
