# -*- coding: utf-8 -*-
"""Paths to data."""
from pathlib import Path

# Top-level directory containing code and data (one level up)
topdir = Path(__file__).absolute().parent.parent

# Modelling results
datadir = topdir / "data"
extdatadir = Path("/media") / Path.home().parts[-1] / "Elements" / "exoplanets" / "thai" / "data"

# Plotting output
plotdir = topdir / "plots"
plotdir.mkdir(parents=True, exist_ok=True)
