# -*- coding: utf-8 -*-
"""Utilities for the ExoCAM output."""
import dask.array as da

from grid import add_cyclic_point_to_da, reverse_along_dim
from names import exocam

__all__ = ("adjust_exocam_grid", "calc_pres_exocam", "calc_virtual_temp_exocam")


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


def calc_pres_exocam(ds, pos="mid"):
    r"""
    Derive a 3D array of pressure for ExoCAM.

    .. math::
        p = a p_{0} + b * p_{s}

    Parameters
    ----------
    ds: xarray.Dataset
        ExoCAM dataset with the required variables: "hya*", "hyb*", "P0", "PS".
    pos: str, optional
        Which points to use, middle ("mid") or interface ("inter").

    Returns
    -------
    darr: xarray.DataArray
        Array of air pressure [Pa].
    """
    assert ds.P0.units == ds.PS.units, "Units of pressure variables should be the same!"
    if pos == "mid":
        coef_a = ds["hyam"]
        coef_b = ds["hybm"]
    elif pos == "inter":
        coef_a = ds["hyai"]
        coef_b = ds["hybi"]
    else:
        raise ValueError(f"`pos` should be one of 'middle', 'inter'; {pos} given")
    pres3d = coef_a * ds.P0 + coef_b * ds.PS
    pres3d = pres3d.rename("air_pressure")
    pres3d.attrs.update({"units": ds.P0.units, "long_name": f"air_pressure_at_{pos}_points"})
    return pres3d


def calc_virtual_temp_exocam(ds, mw_ratio=1.55423618, epsilon=287.058 / 461.52):
    r"""
    Calculate virtual temperature from ExoCAM data.

    .. math::
        T_v = T \frac{(1 + R_v / \epsilon)}{1 + R_v}

    Parameters
    ----------
    ds: xarray.Dataset
        Model dataset with the required variables: "Q" and "T".
    mw_ratio: float, optional
        Ratio of dry air and condensible species molecular weights.
        By default, assumes a N2-dominated atmosphere with H2O as the condensible species.
    epsilon: float, optional
        Ratio of dry air and condensible species gas constants.
        By default, uses 287.058 [J kg-1 K-1] divided by 461.52 [J kg-1 K-1].

    Returns
    -------
    darr: xarray.DataArray
        Array of virtual temperature [K].
    """
    # TODO: make this abstract or switch to MetPy
    # Convert specific humidity to mass mixing ratio
    mmr_dry = ds[exocam.sh] / (1.0 - ds[exocam.sh])
    # Convert to volume mixing ratio
    vmr_dry = mmr_dry * mw_ratio
    # Calculate volume mixing ratio of water relative to all species
    rv = vmr_dry / (1.0 + vmr_dry)
    # Calculate virtual temperature
    temp_v = ds.T * (1 + rv / epsilon) / (1 + rv)
    return temp_v


def calc_alt_exocam(
    ds,
    case="Ben1",
    mw_ratio=1.55423618,
    dry_air_gas_constant=287.058,
    condens_gas_constant=461.52,
    gravity=9.12,
):
    """
    Derive a 3D array of altiude for ExoCAM using the hypsometric equation.

    .. math::
        h = \frac{R T_v}{g} ln{ \frac{p_1}{p_2} }

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with the required variables: "T", "Q", "hyam", "hybm","hyai", "hybi", "P0", "PS".
    mw_ratio: float, optional
        Ratio of dry air and condensible species molecular weights.
        By default, assumes a N2-dominated atmosphere with H2O as the condensible species.
    dry_air_gas_constant: float, optional
        Dry air gas constant [J kg-1 K-1].
    condens_gas_constant: float, optional
        Condensible species gas constant 461.52 [J kg-1 K-1].
    g: float, optional
        Gravity constant [m s-2]. By default, the value for Trappist-1e is used.

    Returns
    -------
    darr: xarray.DataArray
        Array of altitude level [m].
    """
    if case in ["Ben1", "Ben2"]:
        # Use real temperature
        temp_v = ds[exocam.temp]
    else:
        # Calculate virtual temperature
        temp_v = calc_virtual_temp_exocam(
            ds, mw_ratio=mw_ratio, epsilon=dry_air_gas_constant / condens_gas_constant
        )
    # Calculate pressure at mid-level points
    pres_m = calc_pres_exocam(ds, pos="mid")
    # Calculate pressure and interface points
    pres_i = calc_pres_exocam(ds, pos="inter")
    p1 = pres_i[dict(ilev=slice(None, -1))]  # lower interface (higher pressure)
    p2 = pres_i[dict(ilev=slice(1, None))]  # upper interface (lower pressure)
    # Reassign the coordinate to mid-level points to be compatible
    # with the vertical coordinate of`temp_v`
    p1 = p1.rename(ilev=exocam.lev).assign_coords(lev=temp_v.lev)
    p2 = p2.rename(ilev=exocam.lev).assign_coords(lev=temp_v.lev)
    # Calculate the geopotential height thickness of each layer
    dz_int = (dry_air_gas_constant * temp_v / gravity) * da.log(p1 / p2)
    # Stack up layer thicknesses to get total height of lower interface layers
    ilev_alt = dz_int.cumsum(dim=exocam.lev)
    # Discard the last level height and insert zeros in the first level
    ilev_alt = ilev_alt.shift(lev=1).fillna(0.0)
    # Calculate geopotential height thickness between midpoints and interfaces
    dz_mid_int = (dry_air_gas_constant * temp_v / gravity) * da.log(p1 / pres_m)
    # Find midpoint height by adding half-layer thickness to the height of lower interface levels
    lev_alt = ilev_alt + dz_mid_int
    # Reverse level heights back to the original ordering
    lev_alt = lev_alt
    # Update metadata
    lev_alt = lev_alt.rename("altitude")
    lev_alt.attrs.update({"units": "m", "long_name": "altitude_at_mid_points"})
    return lev_alt
