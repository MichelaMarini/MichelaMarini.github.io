"""Constants taken from the original analysis scripts (see original_scripts/)."""
from __future__ import annotations

import re
from dataclasses import dataclass

DMSO = "DMSO"

#: Distance maps: features are ``df.columns[9:]`` (Figures 4-5) ...
FEATURE_START_COL = 9
#: ... and ``df.columns[16:-1]`` for the control wells (Figure 3B).
CONTROL_FEATURE_START_COL = 16

EXCLUDED_TREATMENTS = ("empty",)

#: Cell-line order used in the figures.
CELL_LINE_ORDER = ("HepG2", "MDAMB231", "PC3", "5637", "U2OS", "Panc1", "A549")

#: Figure 4-5, ED heatmap "with order" (``compounds_keep`` in the original).
FIGURE_COMPOUNDS = (
    'Mefloquine', 'Vps34-in1', 'Ym-201636', 'Fccp', 'Rottlerin', 'Toyocamycin', '17-dmag',
    '6-aminonicatinamide', 'Brequinar', 'Erastin', 'Tms', 'Gsk1059615', 'Torin-1', 'Tmp-153',
    '6-thioguanine', 'Gemcitabine', 'Adaptaquin', 'Dorsomorphin', 'Hqno', 'DMSO', 'Nk-252',
    'Spautin-1', 'Gdc-0068', 'Olaparib', 'Copanlisib', 'Methotrexate', 'Swainsonine', 'Gossypol',
    'Artemisinin', 'Mycophenolic-acid', 'Dipyridamole', 'Ro-48-8071',
)
#: ED heatmap 1 (``create_heatmap``): treatments excluded ...
HEATMAP1_EXCLUDE = (
    'Auranofin', 'Sf1670', 'Plumbagin', 'Pr-619', 'Cb-5083', 'Pfk158', 'Eeyarestatin-1',
    'Digitoxin', 'Tg101348', 'Actinomycin d', 'Actinomycin-d', 'Rotenone', 'Pfk15', 'Berb',
    'Actd', 'Berbcl', 'MK2206',
)
#: ... and treatments placed at the bottom, in this order.
HEATMAP1_LAST = ('Vps34-in1', 'Rotenone', 'Paclitaxel', 'Torin-1', 'Adaptaquin', 'Ym-201636',
                 '17-dmag', '6-thioguanine', 'Fccp')
#: ED heatmap 2 (``create_heatmap_2``) cell-line order.
HEATMAP2_CELL_ORDER = ('PC3', '5637', 'MDAMB231', 'Panc1', 'A549', 'U2OS', 'HepG2')

#: Default controls for the replicate-validation plot: the negative control
#: (DMSO) and the positive controls. Names are matched to the treatment names
#: in the files ignoring case, spaces and hyphens ("ActD" = "actd").
NEGATIVE_CONTROL = "DMSO"
DEFAULT_POSITIVE_CONTROLS = ("ActD", "BerbCl", "MK2206")
CONTROLS = (NEGATIVE_CONTROL,) + DEFAULT_POSITIVE_CONTROLS


def ordered(present, preferred=CELL_LINE_ORDER) -> list[str]:
    """The cell lines in ``present``: those in ``preferred`` first, in that
    order, then any others alphabetically. Only cell lines in the data appear."""
    present = list(dict.fromkeys(present))
    known = [c for c in preferred if c in present]
    return known + sorted(c for c in present if c not in known)


def canon(name) -> str:
    """Case-, space- and hyphen-insensitive key used to match compound names."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


@dataclass
class Settings:
    """Parameters of the original scripts (defaults = original values)."""
    min_dmso_cells: int = 100          # DMSO wells below this are dropped
    cell_loss_threshold: float = 50.0  # % of DMSO
    ed_threshold: float = 2.0          # treatments kept if ED >= this anywhere
    ed_vmax: float = 6.9               # ED heatmap colour-scale maximum
    mad_tolerance: float = 1e-6
    #: Order of the cell lines in the figures. None = the cell lines found in
    #: the data, in CELL_LINE_ORDER order (others appended).
    cell_line_order: tuple[str, ...] | None = None
    compounds: tuple[str, ...] = FIGURE_COMPOUNDS
    #: The original used the DMSO medians of the first file in directory
    #: order (alphabetical on Windows). None = first file in sorted order;
    #: a dict maps replicate -> cell line.
    reference_cell_line: str | dict[str, str] | None = None
    #: Positive controls for the replicate-validation plot (the negative
    #: control, DMSO, is always included).
    positive_controls: tuple[str, ...] = DEFAULT_POSITIVE_CONTROLS

    @property
    def controls(self) -> tuple[str, ...]:
        rest = [c for c in self.positive_controls if canon(c) != canon(NEGATIVE_CONTROL)]
        return (NEGATIVE_CONTROL, *rest)
