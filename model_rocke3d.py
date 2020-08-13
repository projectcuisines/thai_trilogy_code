# -*- coding: utf-8 -*-
"""Utilities for the Unified Model output."""
import dask.array as da

import xarray as xr

from grid import reverse_along_dim, roll_da_to_pm180


__all__ = ("adjust_rocke3d_grid")


def adjust_rocke3d_grid(darr, lon_name="lon", lat_name="lat"):
    """
    Adjust the grid of a ROCKE3D data array.
    Reverse the latitude dimension and shift the substellar coordinate
    from -180 degrees to 0 degree in longitude.
    """
    out = darr
    if lat_name in out.dims:
        out = reverse_along_dim(out, lat_name)
    if lon_name in out.dims:
        # Shift data along the longitude to center the substellar at (0,0)
        out = roll_da_to_pm180(
            out.assign_coords(**{lon_name: out[lon_name] + 180}), lon_name=lon_name
        )
    return out

