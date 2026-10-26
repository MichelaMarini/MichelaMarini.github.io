"""Web interface. Run locally with:  streamlit run app/streamlit_app.py"""
from __future__ import annotations

import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import io  # noqa: E402
import threading  # noqa: E402
import tempfile  # noqa: E402
import zipfile  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import streamlit as st  # noqa: E402

from cpmetab import Settings, load_dataset, run_analysis  # noqa: E402
from cpmetab import plots  # noqa: E402
from cpmetab.cli import export  # noqa: E402
from cpmetab import kde  # noqa: E402
from cpmetab.config import (DEFAULT_POSITIVE_CONTROLS, DMSO, NEGATIVE_CONTROL,  # noqa: E402
                            canon, ordered)
from cpmetab.distributions import (FIG6_TREATMENTS, DistributionSettings,  # noqa: E402
                                   count_outliers, feature_reproducibility,
                                   load_single_cell_files, read_wasserstein_table,
                                   run_distribution_analysis)
from cpmetab.synthetic import make_frames  # noqa: E402

st.set_page_config(page_title="Cell Painting metabolic screen", layout="wide")
#: Optional: a Wasserstein table shipped with the app (put it here at publication).
APP_DATA = Path(__file__).parent / "data"


@st.cache_data(show_spinner=False)
def example_files():
    out = []
    for name, df in make_frames().items():
        buf = io.BytesIO()
        df.to_excel(buf, index=False)
        out.append((name, buf.getvalue()))
    return out


@st.cache_data(show_spinner="Reading files…")
def treatments_in(files):
    data = load_dataset([(n, io.BytesIO(b)) for n, b in files])
    return sorted(t for t in data.wells["treatment"].unique() if t != DMSO)


@st.cache_data(show_spinner="Running analysis…")
def analyse(files, positive_controls=DEFAULT_POSITIVE_CONTROLS):
    return run_analysis(load_dataset([(n, io.BytesIO(b)) for n, b in files]),
                        Settings(positive_controls=tuple(positive_controls)))


@st.cache_data(show_spinner="Reading table…")
def read_table(name, data: bytes):
    return read_wasserstein_table(io.BytesIO(data), name)


@st.cache_data(show_spinner="Reading single-cell files (this can take a few minutes)…",
               max_entries=1)
def single_cells(files):
    return load_single_cell_files([(n, io.BytesIO(b)) for n, b in files])   # all compounds


@st.cache_data(show_spinner="Computing Wasserstein distances…", max_entries=2)
def wasserstein_from_cells(files, compounds):
    cells, feats = single_cells(files)
    ds = DistributionSettings(treatments=tuple(compounds))
    return run_distribution_analysis(None, ds, cells[cells["treatment"].isin(compounds)],
                                     feats).distances


def compound_order(names):
    """DMSO first, then the usual Figure 6 compounds, then the rest alphabetically."""
    known = [t for t in FIG6_TREATMENTS if t in names]
    return known + sorted(t for t in names if t not in known)


@st.cache_data(show_spinner="Preparing download…")
def results_zip(files, positive_controls) -> bytes:
    res = analyse(files, positive_controls)
    with tempfile.TemporaryDirectory() as tmp:
        export(res, Path(tmp), formats=("png", "svg"))
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for p in Path(tmp).rglob("*"):
                if p.is_file():
                    z.write(p, p.relative_to(tmp))
    return buf.getvalue()


def excel_bytes(df) -> bytes:
    buf = io.BytesIO()
    df.to_excel(buf, index=False)
    return buf.getvalue()


_PLOT_LOCK = threading.RLock()   # matplotlib is not thread-safe; draw one figure at a time


def show(make):
    """Create a figure with ``make()`` and display it, holding the plot lock
    so that reruns or other sessions cannot draw on it at the same time."""
    with _PLOT_LOCK:
        fig = make()
        st.pyplot(fig, width="content")
        plt.close(fig)


def show_many(make):
    with _PLOT_LOCK:
        figs = make()
        cols = st.columns(min(3, max(len(figs), 1)))
        for i, fig in enumerate(figs.values()):
            with cols[i % len(cols)]:
                st.pyplot(fig, width="content")
            plt.close(fig)


# ------------------------------------------------------------------ sidebar
st.sidebar.header("Distance maps")
uploads = st.sidebar.file_uploader(
    "The 21 wellwise-summary-distance-metrics files (.xlsx)", type="xlsx",
    accept_multiple_files=True)
use_example = st.sidebar.checkbox("Use synthetic example data", value=not uploads)

st.sidebar.header("Single-cell features (optional)")
sc_uploads = st.sidebar.file_uploader(
    "distribution_features_<CellLine>_<rep> files (.xlsx or .csv)", type=["xlsx", "csv"],
    accept_multiple_files=True, key="sc",
    help="Large files: best used on your own computer. Needed for the single-cell "
         "density plots; also used to compute the Wasserstein table.")
sc_files = tuple((f.name, f.getvalue()) for f in sc_uploads) if sc_uploads else None

st.title("Cell Painting screen of metabolic modulators")
st.caption("Uploaded files are processed in memory for this session only.")

res, files = None, None
if uploads:
    files = tuple((f.name, f.getvalue()) for f in uploads)
elif use_example:
    files = tuple(example_files())
    st.info("Distance-map tabs show **synthetic example data** (random values, not real results).")
positive = tuple(DEFAULT_POSITIVE_CONTROLS)
if files:
    try:
        options = treatments_in(files)
        st.sidebar.header("Controls")
        st.sidebar.caption(f"Negative control: **{NEGATIVE_CONTROL}**")
        default = [o for o in options
                   if canon(o) in {canon(c) for c in DEFAULT_POSITIVE_CONTROLS}]
        chosen = st.sidebar.multiselect(
            "Positive controls", options, default=default,
            help="Treatments shown in the replicate-validation plot. "
                 "Default: ActD, BerbCl, MK2206 (if present in the files).")
        names = {canon(c): c for c in DEFAULT_POSITIVE_CONTROLS}   # keep usual spelling
        positive = tuple(names.get(canon(c), c) for c in chosen)
        res = analyse(files, positive)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        st.error(str(e))

order = res.cell_lines if res is not None else []
tabs = st.tabs(["Controls: replicate validation", "Cell count", "Euclidean distance",
                "Pearson correlation", "EMD feature deviation: biological replicates",
                "Single-cell distributions"])


def need_maps():
    st.info("Add the distance-map files in the sidebar (or tick the example data).")


with tabs[0]:
    pos_txt = ", ".join(positive) if positive else "none selected"
    st.markdown(
        f"Negative control ({NEGATIVE_CONTROL}) and positive controls ({pos_txt}): Euclidean "
        "distance of each well's EMD (Wasserstein) profile, per cell line and biological "
        "replicate, to validate technical and biological replicates. Choose the positive "
        "controls in the sidebar.")
    if res is None:
        need_maps()
    else:
        for c in res.settings.controls:
            if c in set(res.control_ed["control"]):
                show(lambda c=c: plots.control_ed_plot(res.control_ed, c, order))
        if res.diagnostics["controls_not_in_data"]:
            st.caption("Not found in the files: " + ", ".join(res.diagnostics["controls_not_in_data"]))

with tabs[1]:
    if res is None:
        need_maps()
    else:
        rep = st.radio("Biological replicate", list(res.toxicity), horizontal=True, key="r1")
        st.markdown(f"**Cell count, biological replicate {rep}** (% of DMSO)")
        show(lambda: plots.toxicity_heatmap(res.toxicity[rep], order))
        st.dataframe(res.toxicity[rep].round(1))

with tabs[2]:
    if res is None:
        need_maps()
    else:
        variant = st.radio("Heatmap", list(res.heatmaps), horizontal=True)
        rep = st.radio("Biological replicate", list(res.toxicity), horizontal=True, key="r2")
        st.markdown(f"**Euclidean distance ({variant}), biological replicate {rep}**")
        show(lambda: plots.ed_heatmap(res.heatmaps[variant][rep], variant, res.settings.ed_vmax))
        st.dataframe(res.heatmaps[variant][rep].round(2))

with tabs[3]:
    if res is None:
        need_maps()
    elif res.correlations.empty:
        st.info("At least two biological replicates are needed for the correlations.")
    else:
        show(lambda: plots.correlation_boxplot(res.correlations, order))
        st.dataframe(res.correlation_summary.round(3))

with tabs[4]:
    st.markdown(
        "For each cell line, treatment and feature, the single-cell distribution of "
        "biological replicate 2 and 3 is compared with replicate 1 (signed Wasserstein "
        "distance, threshold ±0.15).")
    bundled = sorted(list(APP_DATA.glob("wasserstein_distances*.xlsx"))
                     + list(APP_DATA.glob("wasserstein_distances*.csv")))
    wfile = st.file_uploader("Wasserstein table (.xlsx or .csv)", type=["xlsx", "csv"],
                             key="w")
    with st.expander("Where do I get this table?"):
        st.markdown(
            "The table is computed from the single-cell feature files, which are too large "
            "to upload here.\n\n"
            "* **Download it** from the data repository, where it is deposited with the "
            "dataset (`wasserstein_distances.xlsx`), or\n"
            "* **compute it** from the single-cell files with the code in this repository "
            "(`scripts/run_in_spyder.py` or `cpmetab-figures --single-cell-dir ...`); it is "
            "written to `tables/wasserstein_distances.xlsx`.\n\n"
            "The original `wasserstein_distances_with_counts.xlsx` also works.")
    dist = None
    try:
        if wfile is not None:
            dist = read_table(wfile.name, wfile.getvalue())
        elif sc_files:
            avail = compound_order(sorted(single_cells(sc_files)[0]["treatment"].unique()))
            comp6 = st.multiselect("Compounds (compared between biological replicates)", avail,
                                   default=[t for t in FIG6_TREATMENTS if t in avail] or avail,
                                   key="comp6")
            dist = wasserstein_from_cells(sc_files, tuple(comp6)) if comp6 else None
            st.caption("Computed from the single-cell files in the sidebar.")
            if dist is not None:
                st.download_button("Download this Wasserstein table (.xlsx)", excel_bytes(dist),
                                   file_name="wasserstein_distances.xlsx")
        elif bundled:
            dist = read_table(bundled[0].name, bundled[0].read_bytes())
            st.caption(f"Using the table provided with the app: {bundled[0].name}")
    except ValueError as e:
        st.error(str(e))
    if dist is not None:
        in_table = compound_order(sorted(dist["treatment"].unique()))
        shown = st.multiselect(f"Compounds to show ({len(in_table)} in the table)", in_table,
                               default=in_table, key="shown6")
        dist = dist[dist["treatment"].isin(shown)]
    if dist is not None and dist.empty:
        st.info("Select at least one compound.")
        dist = None
    if dist is not None:
        counts = count_outliers(dist)
        treatments = [t for t in in_table if t in shown]
        rep = st.radio("Biological replicate (vs replicate 1)",
                       sorted(counts["replicate"].unique()), horizontal=True, key="r3")
        lines6 = ordered(counts["cell-line"].unique())
        st.markdown(f"**EMD feature deviation counts, biological replicate {rep} vs replicate 1**")
        cols = st.columns(min(4, len(lines6)))
        for i, cl in enumerate(lines6):
            with cols[i % len(cols)]:
                show(lambda: plots.outlier_bars(counts, rep, cl, treatments, show_ylabel=i == 0))
        c_rep = counts[counts["replicate"] == rep]
        st.caption("Number of deviating features (low = W ≤ −0.15, up = 0.15 ≤ W ≤ 2) "
                   f"per cell line and compound, replicate {rep} vs replicate 1:")
        st.dataframe(c_rep.pivot_table(index="cell-line", columns="treatment",
                                       values=["low_count", "up_count"], aggfunc="sum")
                     .reindex(lines6).reindex(columns=treatments, level=1))

        st.subheader("Features: consistent or deviating between biological replicates")
        n_comp = int(dist.groupby("feature")["signed_w"].count().max())
        pct = st.slider("A feature counts as consistent if it deviates in at most this "
                        f"% of the {n_comp} replicate comparisons", 0.0, 100.0, 0.0, 0.5)
        feat = feature_reproducibility(dist, 0.15, pct)
        ok = int((feat["status"] == "consistent").sum())
        st.markdown(f"**{ok}** of **{len(feat)}** features are consistent.")
        st.dataframe(feat, hide_index=True)
        c1, c2 = st.columns(2)
        c1.download_button("Download feature table (.xlsx)", excel_bytes(feat),
                           file_name="emd_feature_reproducibility.xlsx")
        c2.download_button("Download feature table (.csv)", feat.to_csv(index=False).encode(),
                           file_name="emd_feature_reproducibility.csv")

with tabs[5]:
    if not sc_files:
        st.info("Add distribution_features files in the sidebar (Single-cell features). "
                "They are large, so this works best when the app runs on your own computer.")
    else:
        try:
            cells, feats = single_cells(sc_files)
        except ValueError as e:
            st.error(str(e))
            cells = None
        if cells is not None:
            lines7 = ordered(cells["cell-line"].unique())
            feat_list = sorted(feats)
            avail = compound_order(sorted(cells["treatment"].unique()))
            cl = st.selectbox("Cell line", lines7)

            st.markdown("#### Density per well and field of view, per biological replicate")
            c1, c2 = st.columns(2)
            fa = c1.selectbox("Feature", feat_list, key="fa",
                              index=feat_list.index(kde.FIG7A["feature"])
                              if kde.FIG7A["feature"] in feat_list else 0,
                              help="Type to search the feature list.")
            in_cl = [t for t in avail if t in set(cells.loc[cells["cell-line"] == cl, "treatment"])]
            tr = c2.selectbox("Compound", in_cl, key="tr",
                              index=in_cl.index("DMSO") if "DMSO" in in_cl else 0)
            try:
                show_many(lambda: kde.replicate_density_figures(cells, cl, fa, tr, label=tr,
                                                                xlabel=fa))
            except ValueError as e:
                st.warning(str(e))

            st.markdown("#### Ridge plot of compounds")
            c3, c4 = st.columns(2)
            fb = c3.selectbox("Feature", feat_list, key="fb",
                              index=feat_list.index(kde.FIG7B["feature"])
                              if kde.FIG7B["feature"] in feat_list else 0,
                              help="Type to search the feature list.")
            reps7 = sorted(cells.loc[cells["cell-line"] == cl, "replicate"].unique())
            r7 = c4.selectbox("Biological replicate", reps7, key="r7")
            comps = st.multiselect("Compounds", in_cl, key="comps",
                                   default=[t for t in kde.FIG7B["treatments"] if t in in_cl]
                                   or in_cl[:6])
            if comps:
                try:
                    show(lambda: kde.ridge_plot(cells, cl, fb, tuple(comps), replicate=r7,
                                                xlabel=fb))
                except ValueError as e:
                    st.warning(str(e))

if res is not None:
    st.download_button("Download all distance-map tables and figures (.zip)",
                       results_zip(files, positive), file_name="cpmetab_results.zip",
                       mime="application/zip")
