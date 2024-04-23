#!/usr/bin/python
# -----------------
# Load datasets:
# Updated: 06/28/21
# -----------------
from pathlib import Path
import pandas as pd

HERE = Path(__file__).parent


def snap_colors():
    """Load a list of 88 categorical colors for plotting."""
    filename = HERE / "snap_colors.tsv"
    return pd.read_csv(filename, sep="\t", header=None)
