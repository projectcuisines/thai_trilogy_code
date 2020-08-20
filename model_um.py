# -*- coding: utf-8 -*-
"""Utilities for the Unified Model output."""
import dask.array as da

import numpy as np

import xarray as xr

from grid import reverse_along_dim, roll_da_to_pm180


__all__ = (
    "adjust_um_grid",
    "calc_um_ocean_frac",
    "calc_um_rei",
    "calc_um_rel",
    "open_mf_um",
    "prep_um_ds",
)


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


def calc_um_ocean_frac(um_ds, t_freeze=271.35):
    """Calculate the open ocean fraction from the UM dataset."""
    t_sfc = um_ds["STASH_m01s00i024"]  # Surface temperature
    ocean_frac = da.exp(-(t_freeze - t_sfc) / 2)
    ocean_frac = ocean_frac.where(ocean_frac < 1, 1)
    return ocean_frac


def calc_um_rei(um_ds, t_thresh=193.15):
    """Calculate REI from the UM dataset using A. Baran's fit extended to low temperatures."""
    is_cld_ice = xr.where(
        um_ds.STASH_m01s00i012 > 0, 1, 0
    )  # select where is any cloud ice
    air_temp = um_ds.STASH_m01s16i004
    # Cap temperature by t_thresh
    air_temp = xr.where(air_temp > t_thresh, air_temp, t_thresh)
    # Calculate the R_eff
    rei = is_cld_ice * (-353.613 + 1.868 * air_temp) / 2
    # Convert to microns
    rei *= 1e-6
    rei.rename("ice_cloud_condensate_effective_radius")
    rei.attrs = {"long_name": "ice_cloud_condensate_effective_radius", "units": "micron"}
    return rei

def calc_um_rel(um_ds,rho = 1.225):   
	"""     Calculation of the effective radius of liquid water
	following: 
	r_eff = ( 3 * rho_air * qcl / ( 4 * pi * rho_water * 0.8 * 1e8 ))**(1./3.)
	corresponding to Eq 14 
	in https://journals.ametsoc.org/jas/article/51/13/1823/23387/The-Measurement-and-Parameterization-of-Effective
	Parameters
	----------
	liq: xarray.DataArray
		mass mixing ratio of liquid [kg/kg].
	mw_dryair : float, optional
	    Mean molecular weight of dry air [kg mol-1].
	
	 Returns
	-------
	rei: xarray.DataArray
	    Ice effective radius [um], dimension (time,lev,lat,lon).
	"""  
	liq = um_ds.STASH_m01s00i012
	
	rel = (3 * rho * liq / (4*np.pi*1000*0.8*1e8))**(1/3)

	rel.rename("liquid_cloud_condensate_effective_radius")
	rel.attrs = {"long_name": "liquid_cloud_condensate_effective_radius", "units": "micron"}
	
	return rel


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


def prep_um_ds(raw_ds, vert_lev_miss_val="drop"):
    """Prepare UM dataset: interpolate data to common grid, drop redundant coordinates."""
    # Reference grid: use theta-grid for all variables
    lat = "latitude_t"
    lon = "longitude_t"
    eta = "thlev_eta_theta"
    alt = "thlev_zsea_theta"
    coords_to_drop = [
        "rholev_eta_rho",
        "rholev_C_rho",
        "rholev_zsea_rho",
        "rholev_model_level_number",
        "latitude_cu",
        "longitude_cu",
        "latitude_cv",
        "longitude_cv",
    ]
    # Loop over data arrays and adjust their grids
    new_ds = {}
    for d in raw_ds.data_vars:
        if (raw_ds[d].ndim == 1) or (raw_ds[d].name.lower().endswith("_bounds")):
            # Skip `grid_crs` and bound arrays
            continue
        if d == "STASH_m01s00i407":
            # Skip pressure on rho levels
            continue
        if d == "STASH_m01s00i002":
            var = raw_ds[d].interp(
                latitude_cu=raw_ds[lat],
                longitude_cu=raw_ds[lon],
                rholev_eta_rho=raw_ds[eta],
                kwargs={"fill_value": "extrapolate"},
            )
            for coord in coords_to_drop:
                try:
                    var = var.drop_vars(coord)
                except ValueError:
                    pass
        elif d == "STASH_m01s00i003":
            var = raw_ds[d].interp(
                latitude_cv=raw_ds[lat],
                longitude_cv=raw_ds[lon],
                rholev_eta_rho=raw_ds[eta],
                kwargs={"fill_value": "extrapolate"},
            )
            for coord in coords_to_drop:
                try:
                    var = var.drop_vars(coord)
                except ValueError:
                    pass
        else:
            var = raw_ds[d]
        new_ds[d] = adjust_um_grid(var)
    ds = xr.Dataset(new_ds).swap_dims({eta: alt})
    # Extrapolate or drop missing values at the top model level
    new_ds = {}
    if vert_lev_miss_val == "extrapolate":
        for d in ds.data_vars:
            if d in ["STASH_m01s16i004", "STASH_m01s30i113"]:
                # if alt in ds[d].dims:
                new_ds[d] = ds[d].interpolate_na(
                    dim=alt, method="slinear", fill_value="extrapolate"
                )
            else:
                new_ds[d] = ds[d]
    elif vert_lev_miss_val == "drop":
        for d in ds.data_vars:
            if alt in ds[d].dims:
                # Drop the topmost level
                new_ds[d] = ds[d].isel(**{alt: slice(None, -1)})
            else:
                new_ds[d] = ds[d]
    return xr.Dataset(new_ds)
