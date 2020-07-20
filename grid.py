# -*- coding: utf-8 -*-
"""Functionality related to coordinates and grids of xarray Data Arrays."""
import dask.array as da

import numpy as np

import xarray as xr


__all__ = (
    "add_cyclic_point_to_da",
    "reverse_along_dim",
    "roll_da_to_0360",
    "roll_da_to_pm180",
    "wrap_lons",
)


def add_cyclic_point_to_da(darr, coord_name):
    """
    Add a cyclic point to an `xarray.DataArray` and a corresponding
    coordinate. Mimics `cartopy.util.add_cyclic_point()`.

    Parameters
    ----------
    darr: xarray.DataArray
        An n-dimensional `DataArray` of data to add a cyclic point to.
    coord: str
        Coordinate name correspoinding to a 1-dimensional coordinate of `darr`.
        The coordinate values must be regularly spaced.

    Returns
    -------
    cyclic_darr: xarray.DataArray
        The array with a cyclic point added.
    """
    coord = darr[coord_name]
    axis = darr.get_axis_num(coord_name)
    if coord is not None:
        if coord.ndim != 1:
            raise ValueError("The coordinate must be 1-dimensional.")
        delta_coord = np.diff(coord)
        if not da.allclose(delta_coord, delta_coord[0]):
            raise ValueError("The coordinate must be equally spaced.")
        new_coord = xr.concat((coord, coord[-1:] + delta_coord[0]), dim=coord_name)

    slicer = [slice(None)] * darr.ndim
    slicer[axis] = slice(0, 1)
    cyclic_darr = xr.concat((darr, darr[tuple(slicer)]), dim=coord_name)
    cyclic_darr[coord_name] = new_coord
    return cyclic_darr


def reverse_along_dim(darr, dim_name):
    """Reverse `xarray.DataArray` along one of its dimensions."""
    out = darr.reindex(**{dim_name: list(reversed(darr[dim_name]))})
    return out


def roll_da_to_0360(darr, lon_name="longitude"):
    """Roll DataArray's data with longitudes from (-180, 180) to (0, 360)."""
    # Roll longitudes and corresponding data
    out = darr.roll(**{lon_name: darr[lon_name].shape[0] // 2}, roll_coords=True)
    # Reset western (negative) longitudes to values within (180, 360) range
    out[lon_name] = out[lon_name] - (out[lon_name] // 360) * 360
    return out


def roll_da_to_pm180(darr, lon_name="longitude"):
    """Roll DataArray's data with longitudes from (0, 360) to (-180, 180)."""
    # Roll longitudes and corresponding data
    out = darr.roll(**{lon_name: darr[lon_name].shape[0] // 2}, roll_coords=True)
    # Reset western (negative) longitudes to values within (180, 360) range
    out[lon_name] = wrap_lons(out[lon_name], -180, 360)
    return out


def wrap_lons(lons, base, period):
    """Wrap longitude values into the range between base and base+period."""
    # It is important to use 64bit floating precision
    lons = lons.astype(np.float64)
    return ((lons - base + period * 2) % period) + base
