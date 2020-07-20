# -*- coding: utf-8 -*-
"""Utilities for the Unified Model output."""
import dask.array as da

import xarray as xr

from grid import reverse_along_dim, roll_da_to_pm180


__all__ = ("adjust_um_grid", "calc_um_ocean_frac", "open_mf_um")


def adjust_um_grid(darr):
    """
    Adjust the grid of a UM data array.

    Reverse the latitude dimension and shift the data to +/-180 degrees in longitude.
    """
    out = darr
    try:
        lat_name = [i for i in darr.dims if "latitude" in i.lower()][0]
        out = reverse_along_dim(out, lat_name)
    except IndexError:
        pass
    try:
        lon_name = [i for i in darr.dims if "longitude" in i.lower()][0]
        # Shift data along the longitude
        out = roll_da_to_pm180(out, lon_name=lon_name)
    except IndexError:
        pass
    return out


def calc_um_ocean_frac(um_ds, t_freeze=273.15):
    """Calculate the open ocean fraction from the UM dataset."""
    t_sfc = um_ds["STASH_m01s00i024"]  # Surface temperature
    ocean_frac = da.exp(-(t_freeze - t_sfc) / 2)
    ocean_frac = ocean_frac.where(ocean_frac < 1, 1)
    return ocean_frac


def open_mf_um(files, main_time, rad_time, **kw_open):
    """
    Open multiple UM output files and do the time dim concatenation.

    The function uses `xarray.open_dataset()` to open files, chunking
    the data along the main (non-radiation) time dimension; and interpolates
    data with the "radiation" time dimension to the main one.

    Parameters
    ----------
    files: list of str
        List of files to open.
    main_time: str
        Name of the main (non-radiation) time coordinate.
        For example, "hourly".
    rad_time: str
        Name of the radiation time coordinate.
        For example, "hourly_rad"
    **kw_open: dict, optional
        Keyword arguments passed to `xarray.open_dataset()`

    Returns
    -------
    xarray.Dataset
        Concatenated dataset with unified time dimension.
    """
    dsets = []
    for fname in files:
        ds = xr.open_dataset(fname, **kw_open,)
        proc_ds = {}
        target_time = ds[main_time]
        for d in ds.data_vars:
            if rad_time in ds[d].dims:
                var = ds[d].interp(
                    **{rad_time: target_time}, kwargs={"fill_value": "extrapolate"}
                )
                var = var.drop_vars(rad_time)
            else:
                var = ds[d]
            proc_ds[d] = var
        proc_ds = xr.Dataset(proc_ds)
        dsets.append(proc_ds)
    return xr.concat(dsets, main_time)
