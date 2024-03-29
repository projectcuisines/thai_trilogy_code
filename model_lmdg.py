# -*- coding: utf-8 -*-
"""Utilities for the ExoCAM output."""
import dask.array as da

from names import lmdg

__all__ = ("calc_alt_lmdg", "calc_virtual_temp_lmdg")


def calc_virtual_temp_lmdg(ds, mw_ratio=1.55423618, epsilon=287.058 / 461.52):
    r"""
    Calculate virtual temperature from LMDG data.

    .. math::
        T_v = T \frac{(1 + R_v / \epsilon)}{1 + R_v}

    Parameters
    ----------
    ds: xarray.Dataset
        Model dataset with the required variables: "h2o_vap" and "temp".
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
    # "specific humidity" is actually mass mixing ratio in LMD-G
    mmr_dry = ds[lmdg.sh]
    # Convert to volume mixing ratio
    vmr_dry = mmr_dry * mw_ratio
    # Calculate volume mixing ratio of water relative to all species
    rv = vmr_dry / (1.0 + vmr_dry)
    # Calculate virtual temperature
    temp_v = ds[lmdg.temp] * (1 + rv / epsilon) / (1 + rv)
    return temp_v


def calc_alt_lmdg(
    ds,
    case="Ben1",
    mw_ratio=1.55423618,
    dry_air_gas_constant=287.058,
    condens_gas_constant=461.52,
    gravity=9.12,
):
    r"""
    Derive a 3D array of altiude for LMD-G using the hypsometric equation.
    .. math::
        h = \frac{R T_v}{g} ln{ \frac{p_1}{p_2} }
    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with the required variables:
        "temp", "h2o_vap", "aps", "bps", "ap", "bp", "ps".
    case: str, optional
        THAI case.
    mw_ratio: float, optional
        Ratio of dry air and condensible species molecular weights.
        By default, assumes a N2-dominated atmosphere with H2O as the condensible species.
    dry_air_gas_constant: float, optional
        Dry air gas constant [J kg-1 K-1].
    condens_gas_constant: float, optional
        Condensible species gas constant 461.52 [J kg-1 K-1].
    gravity: float, optional
        Gravity constant [m s-2]. By default, the value for Trappist-1e is used.
    Returns
    -------
    darr: xarray.DataArray
        Array of altitude level [m].
    """
    if case in ["Ben1", "Ben2"]:
        # Use real temperature
        temp_v = ds[lmdg.temp]
    else:
        # Calculate virtual temperature
        temp_v = calc_virtual_temp_lmdg(
            ds, mw_ratio=mw_ratio, epsilon=dry_air_gas_constant / condens_gas_constant
        )
    # Calculate pressure at midlayers
    pres_m = ds.aps + ds.bps * ds.ps
    # Calculate pressure at interlayers
    pres_i = ds.ap + ds.bp * ds.ps
    p1 = pres_i[dict(interlayer=slice(None, -1))]  # lower interface (higher pressure)
    p2 = pres_i[dict(interlayer=slice(1, None))]  # upper interface (lower pressure)
    # Reassign the coordinate to midlayers to be compatible
    # with the vertical coordinate of `temp_v`
    p1 = p1.rename(interlayer=lmdg.z).assign_coords(altitude=temp_v.altitude)
    p2 = p2.rename(interlayer=lmdg.z).assign_coords(altitude=temp_v.altitude)
    # Calculate the geopotential height thickness of each layer
    dz_int = (dry_air_gas_constant * temp_v / gravity) * da.log(p1 / p2)
    # Stack up layer thicknesses to get total height of lower interface layers
    ilev_alt = dz_int.cumsum(dim=lmdg.z)
    # Discard the last level height and insert zeros in the first level
    ilev_alt = ilev_alt.shift(altitude=1).fillna(0.0)
    # Calculate geopotential height thickness between midpoints and interfaces
    dz_mid_int = (dry_air_gas_constant * temp_v / gravity) * da.log(p1 / pres_m)
    # Find midpoint height by adding half-layer thickness to the height of lower interface levels
    lev_alt = ilev_alt + dz_mid_int
    # Update metadata
    lev_alt = lev_alt.rename(lmdg.lev)
    lev_alt.attrs.update({"units": "m", "long_name": "altitude_at_mid_points"})
    return lev_alt
