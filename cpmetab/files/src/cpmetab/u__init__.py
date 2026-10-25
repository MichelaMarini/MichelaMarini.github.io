"""cpmetab: analysis of the Cell Painting metabolic-modulator dataset.

Cell-count (toxicity) annotation, Euclidean distance (ED) of directional-EMD
profiles from DMSO, and replicate reproducibility, as used for the figures of
the accompanying Scientific Data descriptor.
"""
from .config import Settings
from .io import load_dataset
from .pipeline import run_analysis, Results

__all__ = ["Settings", "load_dataset", "run_analysis", "Results"]
__version__ = "0.1.0"
