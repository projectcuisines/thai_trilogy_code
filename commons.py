"""Definitions and objects commonly used between scripts."""
from aeolus.region import Region


# Tidally-locked setups
DAYSIDE = Region(-88.74, 91.24, -90, 90, "dayside")
NIGHTSIDE = Region(91.25, -88.75, -90, 90, "nightside")
_NIGHTSIDE_PREC = Region(91.25, -88.75, -90, 87, "isr_flux_zero")
