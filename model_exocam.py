# -*- coding: utf-8 -*-
"""Utilities for the ExoCAM output."""
from grid import add_cyclic_point_to_da, reverse_along_dim


__all__ = ("adjust_exocam_grid", "calc_pres_exocam", )


def adjust_exocam_grid(darr, lat_name="lat", lon_name="lon"):
    """
    Adjust the grid of an ExoCAM data array.

    Reverse the latitude and vertical level dimensions, insert a
    cyclic column of data at 360 degrees longitude.
    """
    out = darr
    try:
        lev_name = "lev"
        out = reverse_along_dim(out, lev_name)
    except KeyError:
        try:
            lev_name = "ilev"
            # For variables like "hyai"
            out = reverse_along_dim(out, lev_name)
        except KeyError:
            pass

    if lat_name in darr.dims:
        out = reverse_along_dim(out, lat_name)
    if lon_name in darr.dims:
        # Add a cyclic point at 360 degree longitude.
        out = add_cyclic_point_to_da(out, lon_name)
        # Shift longitudes
        # (no need to shift the data, since they have the substellar point in the middle)
        out[lon_name] = out[lon_name] - 180
    return out


def calc_pres_exocam(ds):
    r"""
    Derive a 3D array of pressure for ExoCAM.

    .. math::
        p = a p_{0} + b * p_{s}

    Parameters
    ----------
    ds: xarray.Dataset
        ExoCAM dataset with the required variables: "hyam", "hybm", "P0", "PS".

    Returns
    -------
    darr: xarray.DataArray
        Array of air pressure [Pa].
    """
    assert ds.P0.units == ds.PS.units, "Units of pressure variables should be the same!"
    pres3d = ds.hyam * ds.P0 + ds.hybm * ds.PS
    pres3d.rename("air_pressure")
    pres3d.attrs = {
        "units": ds.P0.units,
        "long_name": "air_pressure"
    }
    return pres3d
