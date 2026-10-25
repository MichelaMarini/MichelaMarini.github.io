"""Command line: regenerate all tables and figures.

    cpmetab-figures DISTANCE_MAPS_DIR [--single-cell-dir DIR] -o results/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from . import plots
from .config import DEFAULT_POSITIVE_CONTROLS, Settings, ordered
from .io import load_dataset
from .pipeline import Results, run_analysis


def _saver(out: Path, formats, dpi, written):
    if isinstance(formats, str):          # "png" -> ("png",)
        formats = (formats,)
    formats = [f.strip(". ").lower() for f in formats]
    fdir = out / "figures"
    fdir.mkdir(parents=True, exist_ok=True)

    def save(fig, stem):
        for fmt in formats:
            p = fdir / f"{stem}.{fmt}"
            fig.savefig(p, dpi=dpi, bbox_inches="tight")
            written.append(p)
        plt.close(fig)
    return save


def export(res: Results, out: Path, formats=("png", "svg"), dpi: int = 300) -> list[Path]:
    """Figures 3B, 4-5 and Table 2: tables and figures."""
    s = res.settings
    lines = res.cell_lines
    written: list[Path] = []
    tdir = out / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    save = _saver(out, formats, dpi, written)

    def table(df, name, index=True):
        df.to_csv(tdir / name, index=index)
        written.append(tdir / name)

    table(res.cell_counts, "cell_count_pct_dmso.csv", index=False)
    table(res.ed, "euclidean_distance_all.csv", index=False)
    table(res.correlations, "replicate_correlations.csv", index=False)
    table(res.correlation_summary, "replicate_correlation_summary.csv")
    table(res.control_ed, "controls_replicate_validation_ed.csv", index=False)
    for rep, m in res.toxicity.items():
        table(m, f"cell_count_heatmap_rep{rep}.csv")
        save(plots.toxicity_heatmap(m, lines), f"cell_count_heatmap_rep{rep}")
    for variant, mats in res.heatmaps.items():
        for rep, m in mats.items():
            table(m, f"ed_{variant}_rep{rep}.csv")
            save(plots.ed_heatmap(m, variant, s.ed_vmax), f"ed_{variant}_rep{rep}")
    if len(res.correlations):
        save(plots.correlation_boxplot(res.correlations, lines), "replicate_correlation_boxplot")
    for label in s.controls:
        if label in set(res.control_ed["control"]):
            save(plots.control_ed_plot(res.control_ed, label, lines),
                 f"controls_replicate_validation_{label}")
    p = out / "run_info.json"
    p.write_text(json.dumps(res.diagnostics, indent=2, default=str))
    written.append(p)
    return written


def export_distributions(dres, out: Path, cell_line_order=None, formats=("png", "svg"),
                         dpi: int = 300) -> list[Path]:
    """EMD feature deviation counts between biological replicates: tables
    (csv and xlsx), the per-feature summary, and bar plots."""
    from .distributions import feature_reproducibility, to_original_layout
    written: list[Path] = []
    tdir = out / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    save = _saver(out, formats, dpi, written)

    def both(df, stem):
        df.to_csv(tdir / f"{stem}.csv", index=False)
        df.to_excel(tdir / f"{stem}.xlsx", index=False)
        written.extend([tdir / f"{stem}.csv", tdir / f"{stem}.xlsx"])

    both(dres.distances, "wasserstein_distances")
    both(dres.counts, "emd_feature_deviation_counts")
    feat = feature_reproducibility(dres.distances, dres.settings.threshold)
    both(feat, "emd_feature_reproducibility")
    p = tdir / "wasserstein_distances_with_counts.xlsx"     # original layout
    with pd.ExcelWriter(p) as xw:
        for sheet, df in to_original_layout(dres.distances, dres.counts).items():
            df.to_excel(xw, sheet_name=sheet, index=False)
    written.append(p)
    treatments = [t for t in dres.settings.treatments if t in set(dres.counts["treatment"])]
    lines = ordered(dres.counts["cell-line"].unique())
    for rep in sorted(dres.counts["replicate"].astype(str).unique()):
        for i, cl in enumerate(lines):
            save(plots.outlier_bars(dres.counts, rep, cl, treatments, show_ylabel=i == 0,
                                    show_legend=i == len(lines) - 1),
                 f"emd_feature_deviation_counts_rep{rep}_{cl}")
    return written


def export_kde(cells, feats, out: Path, cell_line_order=None, formats=("png", "svg"),
               dpi: int = 300) -> list[Path]:
    """Single-cell density plots: per replicate (7a) and ridge plots (7b),
    for every cell line in the data."""
    from . import kde
    written: list[Path] = []
    save = _saver(out, formats, dpi, written)
    a, b = kde.FIG7A, kde.FIG7B
    lines = ordered(cells["cell-line"].unique())
    if a["feature"] in feats:
        for cl in lines:
            try:
                for k, fig in kde.replicate_density_figures(cells, cl).items():
                    save(fig, f"single_cell_density_{cl}_{a['treatment']}_{k}")
            except ValueError as e:
                print(f"  note: {e}")
    else:
        print(f"  note: density plots skipped ({a['feature']} not in data)")
    if b["feature"] in feats:
        for cl in lines:
            for rep in sorted(cells.loc[cells["cell-line"] == cl, "replicate"].unique()):
                try:
                    save(kde.ridge_plot(cells, cl, replicate=rep), f"single_cell_ridge_{cl}_rep{rep}")
                except ValueError as e:
                    print(f"  note: {e}")
    else:
        print(f"  note: ridge plots skipped ({b['feature']} not in data)")
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?",
                    help="folder with the wellwise-summary-distance-metrics_<CellLine>_<rep>.xlsx files")
    ap.add_argument("--single-cell-dir", default=None,
                    help="Step4_Features experiment folders, or the 'Distribution Features' "
                         "folder (Figures 6-7)")
    ap.add_argument("--positive-controls", nargs="+", default=list(DEFAULT_POSITIVE_CONTROLS),
                    help="positive controls for the replicate-validation plot "
                         "(default: ActD BerbCl MK2206; DMSO is always the negative control)")
    ap.add_argument("-o", "--out", default="results", help="output folder")
    ap.add_argument("--formats", nargs="+", default=["png", "svg"])
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args(argv)
    if not a.data_dir and not a.single_cell_dir:
        ap.error("give DATA_DIR and/or --single-cell-dir")
    import matplotlib
    matplotlib.use("Agg")
    out = Path(a.out)
    if a.data_dir:
        res = run_analysis(load_dataset(a.data_dir),
                           Settings(positive_controls=tuple(a.positive_controls)))
        if res.diagnostics["controls_not_in_data"]:
            print(f"  note: controls not in data: {res.diagnostics['controls_not_in_data']}")
        files = export(res, out, a.formats, a.dpi)
        print(f"Distance-map analysis: wrote {len(files)} files to {out}")
        if res.diagnostics["compounds_not_in_data"]:
            print(f"  note: compounds not in data: {res.diagnostics['compounds_not_in_data']}")
    if a.single_cell_dir:
        from . import kde
        from .distributions import (DistributionSettings, load_single_cell,
                                    run_distribution_analysis)
        ds = DistributionSettings()
        wanted = tuple(dict.fromkeys(ds.treatments + kde.FIG7B["treatments"]
                                     + (kde.FIG7A["treatment"],)))
        cells, feats = load_single_cell(a.single_cell_dir, wanted)
        dres = run_distribution_analysis(None, ds, cells[cells["treatment"].isin(ds.treatments)],
                                         feats)
        files = export_distributions(dres, out, None, a.formats, a.dpi)
        files += export_kde(cells, feats, out, None, a.formats, a.dpi)
        print(f"Single-cell analysis: wrote {len(files)} files to {out}")


if __name__ == "__main__":
    main()
