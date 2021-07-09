# -*- coding: utf-8 -*-
"""Science diagnostics."""
from functools import partial

import dask.array as da
import iris
from iris.experimental import stratify
import numpy as np
import xarray as xr

from grid import EARTH_RADIUS, grid_cell_areas, reverse_along_dim
from names import names

__all__ = (
    "bond_albedo",
    "brunt_vaisala_frequency",
    "cloud_area_fraction",
    "cloud_mmr_ice",
    "cloud_mmr_liquid",
    "cloud_volume_fraction_ice",
    "cloud_volume_fraction_liquid",
    "cloud_volume_fraction_total",
    "cloud_path_ice",
    "cloud_path_liquid",
    "cloud_path_total",
    "cre_toa",
    "dayside_mean",
    "dry_lapse_rate",
    "get_time_rel_days",
    "global_mean",
    "greenhouse_effect",
    "integral",
    "hdiv",
    "mass_weighted_vertical_integral",
    "meridional_mean",
    "moist_static_energy",
    "nightside_mean",
    "nondim_rossby_deformation_radius",
    "potential_temperature",
    "rossby_deformation_radius_isothermal",
    "rossby_deformation_radius_stratified",
    "scale_height",
    "sfc_dn_lw_flux",
    "sfc_net_up_lw_flux",
    "sfc_temp",
    "specific_humidity",
    "spatial_mean",
    "spatial_sum",
    "terminator_mean",
    "time_mean",
    "toa_olr",
    "upper_atm_vap_mean",
    "vert_mer_mean_of_mse_flux",
    "wind_rot_div",
    "zonal_mass_streamfunction",
    "zonal_mean",
)

INTERPOLATOR = partial(
    stratify.stratify.interpolate,
    interpolation=stratify.stratify.INTERPOLATE_LINEAR,
    extrapolation=stratify.stratify.EXTRAPOLATE_LINEAR,
)


def _trapz(y, x, axis):
    if axis < 0:
        axis = y.ndim + axis
    x_sl1 = [slice(None)] * len(x.shape)
    x_sl1[axis] = slice(1, None)
    x_sl2 = [slice(None)] * len(x.shape)
    x_sl2[axis] = slice(None, -1)
    slice1 = (slice(None),) * axis + (slice(1, None),)
    slice2 = (slice(None),) * axis + (slice(None, -1),)
    dx = x[tuple(x_sl1)] - x[tuple(x_sl2)]
    integrand = dx * 0.5 * (y[tuple(slice1)] + y[tuple(slice2)])
    return xr.core.duck_array_ops.sum(integrand, axis=axis, skipna=False)


def _integrate_generic(ds, coord, dim):
    if dim not in ds.variables and dim not in ds.dims:
        raise ValueError(f"Dimension {dim} does not exist.")

    coord_var = coord.variable

    variables = {}
    coord_names = set()
    for k, v in ds.variables.items():
        if k in ds.coords:
            if dim not in v.dims:
                variables[k] = v
                coord_names.add(k)
        else:
            if k in ds.data_vars and dim in v.dims:
                integ = _trapz(v.data, coord_var.data, axis=v.get_axis_num(dim))
                v_dims = list(v.dims)
                v_dims.remove(dim)
                variables[k] = xr.Variable(v_dims, integ)
            else:
                variables[k] = v
    indexes = {k: v for k, v in ds.indexes.items() if k in variables}
    return ds._replace_with_new_dims(variables, coord_names=coord_names, indexes=indexes)


def bond_albedo(ds, model_key):
    r"""
    Calculate Bond albedo.

    .. math::
        \alpha_b = \frac{OSR_{TOA}}{ISR_{TOA}}

    Parameters
    ----------
    ds: xarray.Dataset
        Input dataset containing relevant variables.
    model_key: str,
        Model name.

    Returns
    -------
    xarray.DataArray
    """
    model_names = names[model_key]
    if model_key == "ExoCAM":
        toa_osr = ds[model_names.toa_isr] - ds[model_names.toa_net_sw]
    elif model_key == "LMDG":
        toa_osr = ds[model_names.toa_isr] - ds[model_names.toa_net_sw]
    elif model_key == "ROCKE3D":
        toa_osr = ds[model_names.toa_isr] - ds[model_names.toa_net_sw]
    elif model_key == "UM":
        toa_osr = ds[model_names.toa_osr]
    alb = toa_osr / ds[model_names.toa_isr]
    return alb


def brunt_vaisala_frequency(
    temp,
    press,
    gas_constant=287.058,
    c_p=1039,
    p_ref=100_000,
    gravity=9.80665,
    lev_name="altitude",
):
    """
    Calculate Brunt–Väisälä frequency.

    Parameters
    ----------
    temp : xarray.DataArray
        Atmospheric temperature [K].
    press : xarray.DataArray
        Atmospheric pressure [Pa].
    gas_constant : float, optional
        Specific gas constant [J kg-1 K-1].
    c_p: float, optional
        Dry air specific heat capacity [m2 s-2 K-1].
    p_ref : float, optional
        Standard reference pressure [Pa].
    gravity: float, optional
        Gravity constant [m s-2].
    lev_name: str, optional
        Name of the vertical coordinate.

    Returns
    -------
    bv_freq: xarray.DataArray
        Brunt-Väisälä frequency [s-1]
    """
    # Compute potential temperature from real temperature and pressure
    theta = potential_temperature(
        temp,
        press,
        gas_constant=gas_constant,
        c_p=c_p,
        p_ref=p_ref,
    )

    bv_freq = ((gravity / theta) * theta.differentiate(lev_name)) ** 0.5
    bv_freq = bv_freq.rename("brunt_vaisala_frequency")
    bv_freq.attrs.update({"units": "s-1"})

    return bv_freq


def cloud_area_fraction(ds, model_key):
    """Extract cloud fraction from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "ROCKE3D":
        out = ds[model_names.caf]
    else:
        # input in [0-1]
        out = ds[model_names.caf] * 100
    return out


def cloud_mmr_ice(ds, model_key):
    """Extract ice cloud MMR on levels from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "LMDG":
        t_min = 258
        t_max = 273
        scaling = (np.clip(ds[model_names.temp], t_min, t_max) - t_min) / (t_max - t_min)
        out = ds[model_names.cld_ice_mf] * (1 - scaling)
    else:
        out = ds[model_names.cld_ice_mf]
    return out


def cloud_mmr_liquid(ds, model_key):
    """Extract liquid cloud MMR on levels from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "LMDG":
        t_min = 258
        t_max = 273
        scaling = (np.clip(ds[model_names.temp], t_min, t_max) - t_min) / (t_max - t_min)
        out = ds[model_names.cld_ice_mf] * scaling
    else:
        out = ds[model_names.cld_liq_mf]
    return out


def cloud_volume_fraction_total(ds, model_key):
    """Extract cloud fraction on levels from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "ROCKE3D":
        out = ds[model_names.cld_liq_v] + ds[model_names.cld_ice_v]
    else:
        out = ds[model_names.cld_v]
    return out * 100


def cloud_volume_fraction_ice(ds, model_key):
    """Extract ice cloud fraction on levels from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "ExoCAM":
        out = ds[model_names.cld_v]
    elif model_key == "LMDG":
        out = ds[model_names.cld_v]
    else:
        out = ds[model_names.cld_ice_v]
    return out * 100


def cloud_volume_fraction_liquid(ds, model_key):
    """Extract liquid cloud fraction on levels from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "ExoCAM":
        out = ds[model_names.cld_v]
    elif model_key == "LMDG":
        out = ds[model_names.cld_v]
    else:
        out = ds[model_names.cld_liq_v]
    return out * 100


def cloud_path_ice(ds, model_key):
    """Extract ice water path from a THAI dataset."""
    model_names = names[model_key]
    if model_key in ["ExoCAM", "ROCKE3D"]:
        # input in [g m-2]
        out = ds[model_names.iwp] / 1000
    elif model_key == "LMDG":
        out = ds[model_names.cwp]
    elif model_key == "UM":
        out = ds[model_names.iwp]
    return out


def cloud_path_liquid(ds, model_key):
    """Extract ice water path from a THAI dataset."""
    model_names = names[model_key]
    if model_key in ["ExoCAM", "ROCKE3D"]:
        # input in [g m-2]
        out = ds[model_names.lwp] / 1000
    elif model_key == "LMDG":
        out = ds[model_names.cwp]
    elif model_key == "UM":
        out = ds[model_names.lwp]
    return out


def cloud_path_total(ds, model_key):
    """Extract total cloud water path from a THAI dataset."""
    model_names = names[model_key]
    if model_key in ["ExoCAM", "ROCKE3D"]:
        # input in [g m-2]
        out = (ds[model_names.lwp] + ds[model_names.iwp]) / 1000
    elif model_key == "LMDG":
        out = ds[model_names.cwp]
    elif model_key == "UM":
        out = ds[model_names.lwp] + ds[model_names.iwp]
    return out


def cre_toa(ds, model_key, kind="total"):
    r"""
    Calculate domain-average TOA cloud radiative effect (CRE).

    .. math::
        CRE_{TOA} = F_{up,clear-sky} - F_{up,all-sky}

    Parameters
    ----------
    ds: xarray.Dataset
        Input dataset containing relevant variables.
    model_key: str,
        Model name.
    kind: str, optional
        Shortwave ('sw'), longwave ('lw'), or 'total' CRE.

    Returns
    -------
    xarray.DataArray
    """
    name = f"toa_cloud_radiative_effect_{kind}"
    model_names = names[model_key]
    if kind == "sw":
        if model_key == "ExoCAM":
            out = -(ds[model_names.toa_net_sw_cs] - ds[model_names.toa_net_sw])
        elif model_key == "LMDG":
            out = -(ds[model_names.toa_net_sw_cs] - ds[model_names.toa_net_sw])
        elif model_key == "ROCKE3D":
            # ds.swup_toa_clrsky - (ds.incsw_toa - ds.srnf_toa)
            out = ds[model_names.toa_osr_cs] - (
                ds[model_names.toa_isr] - ds[model_names.toa_net_sw]
            )
        elif model_key == "UM":
            out = ds[model_names.toa_osr_cs] - ds[model_names.toa_osr]
    elif kind == "lw":
        if model_key == "ExoCAM":
            out = ds[model_names.toa_net_lw_cs] - ds[model_names.toa_net_lw]
        elif model_key == "LMDG":
            out = ds[model_names.toa_olr_cs] - ds[model_names.toa_olr]
        elif model_key == "ROCKE3D":
            out = ds[model_names.toa_crf_lw]
        elif model_key == "UM":
            out = ds[model_names.toa_olr_cs] - ds[model_names.toa_olr]
    elif kind == "total":
        sw = cre_toa(ds, model_key, "sw")
        lw = cre_toa(ds, model_key, "lw")
        out = sw + lw
        out = out.rename(name)
        return out

    out = out.rename(name)
    return out


def dry_lapse_rate(ds, model_key):
    """Compute a lapse rate from an n-dimensional THAI dataset."""
    model_names = names[model_key]
    if model_key == "ExoCAM":
        alt = ds[model_names.z]
        coord = model_names.lev
    elif model_key == "LMDG":
        alt = ds[model_names.lev]
        coord = model_names.z
    elif model_key == "ROCKE3D":
        alt = ds[model_names.z]
        coord = model_names.lev
    elif model_key == "UM":
        alt = ds[model_names.z]
        coord = model_names.z

    lr = ds[model_names.temp].differentiate(coord) / alt.differentiate(coord)
    return lr


def greenhouse_effect(ds, model_key, const, kind="all_sky"):
    r"""
    Calculate the greenhouse effect [K].

    Parameters
    ----------
    ds: xarray.Dataset
        Input dataset containing relevant variables.
    model_key: str,
        Model name.
    kind: str, optional
        Type of GHE:  "all_sky" or "clear_sky"

    Returns
    -------
    xarray.DataArray
    """
    if kind == "all_sky":
        if model_key == "ExoCAM":
            olr = ds[names[model_key].toa_net_lw]
        elif model_key == "LMDG":
            olr = ds[names[model_key].toa_olr]
        elif model_key == "ROCKE3D":
            olr = ds[names[model_key].toa_olr_cs] - ds[names[model_key].toa_crf_lw]
        elif model_key == "UM":
            olr = ds[names[model_key].toa_olr]
    elif kind == "clear_sky":
        if model_key == "ExoCAM":
            olr = ds[names[model_key].toa_net_lw_cs]
        elif model_key == "LMDG":
            olr = ds[names[model_key].toa_olr_cs]
        elif model_key == "ROCKE3D":
            olr = ds[names[model_key].toa_olr_cs]
        elif model_key == "UM":
            olr = ds[names[model_key].toa_olr_cs]

    t_sfc = ds[names[model_key].t_sfc]
    if model_key == "ROCKE3D":
        t_sfc = t_sfc.copy() + const.t_melt  # convert from degC to K

    out = t_sfc - (olr / const.stefan_boltzmann) ** 0.25
    return out


def hdiv(i_arr, j_arr, lon_name="longitude", lat_name="latitude", r_planet=EARTH_RADIUS):
    r"""
    Calculate horizontal divergence of two components of a vector as `xarray.DataArray`s.

    Parameters
    ----------
    i_arr: xarray.DataArray
        i-th component.
    j_arr: xarray.DataArray
        j-th component.
    lon_name: str, optional
        Name of x-coordinate
    lat_name: str, optional
        Name of y-coordinate
    r_planet: float, optional
        Radius of the planet (m). Default is Earth's radius.

    Returns
    -------
    h_div: xarray.DataArray
        Array of horizontal divergence.

    Notes
    -----
    Divergence in spherical coordinates is defined as

    .. math::

        \nabla\cdot \vec A = \frac{1}{r cos \phi} (
        \frac{\partial \vec A_\lambda}{\partial \lambda}
        + \frac{\partial}{\partial \phi}
        (\vec A_\phi cos \phi))

    where \lambda is longitude, \phi is latitude.
    """
    lon_rad = da.deg2rad(i_arr[lon_name])
    lat_rad = da.deg2rad(i_arr[lat_name])
    cos_lat = da.cos(lat_rad)

    # i-component: \frac{\partial \vec A_\lambda}{\partial \lambda}
    di_dlambda = i_arr.diff(lon_name) / lon_rad.diff(lon_name)

    # j-component: \frac{\partial}{\partial \phi} (\vec A_\phi cos \phi))
    djcos_dphi = (j_arr * cos_lat).diff(lat_name) / lat_rad.diff(lat_name)

    # Sum the components and divide by {r cos \phi}
    h_div = (di_dlambda + djcos_dphi) / (r_planet * cos_lat)
    # Interpolate to the original grid for consistency
    h_div = h_div.interp(
        **{lat_name: i_arr[lat_name], lon_name: i_arr[lon_name]},
        kwargs={"fill_value": "extrapolate"},
    )
    h_div = h_div.rename("horizontal_divergence")
    return h_div


def integral(xr_da, dim, coord=None, datetime_unit=None):
    """
    Integrate an `xarray.DataArray` over its dimension(s) or an external N-dim coordinate.

    A hack to extend `xarray.DataArray.integrate()` to a more general case.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Array to integrate.
    dim: hashable, or a sequence of hashable
        Dimension(s) used for the integration.
    coord: xarray.DataArray, optional
        External N-dimensional coordinate for integration.
    datetime_unit: str, optional
        Can be used to specify the unit if datetime coordinate is used.
        One of {'Y', 'M', 'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', 'ps',
        'fs', 'as'}

    Returns
    -------
    result: xarray.DataArray

    See also
    --------
    xarray.DataArray.integrate: xarray function used when `coord` is None
    """
    if coord is None:
        return xr_da.integrate(coord=dim, datetime_unit=datetime_unit)
    else:
        name = xr_da.name
        coord_name = coord.name
        units = xr_da.attrs.get("units", None)
        coord_units = coord.attrs.get("units", None)
        tmp_ds = xr_da._to_temp_dataset()
        if isinstance(dim, (list, tuple)):
            raise ValueError(
                f"Only 1 dim is allowed when using an external array for integration, {dim} given"
            )
        if dim not in coord.dims:
            raise ValueError(f"{coord} does not have {dim} dimension.")
        if datetime_unit is not None:
            raise ValueError(f"Using {coord} with {datetime_unit} is not allowed.")
        result = _integrate_generic(tmp_ds, coord, dim)
        result = result.to_array().squeeze().drop_vars("variable")
        if name is not None and coord_name is not None:
            result = result.rename(f"integral_of_{name}_wrt_{coord_name}")
        if units is not None and coord_units is not None:
            result.attrs["units"] = f"{units} {coord_units}"
        return result


def mass_weighted_vertical_integral(
    xr_da, dim, coord=None, coord_type=None, rho=None, gravity=None
):
    """
    Calculate a vertical integral with mass-weighting.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Array to integrate.
    dim: hashable
        Dimension to use for the integration.
    coord: xarray.DataArray, optional
        Array of a coordinate to use for vertical integration.
    coord_type: str, optional
        Type of vertical coordinate ("height" or "pressure").
    rho: xarray.DataArray, optional
        Array of air density [kg m-3]. Required if `zcoord_type="height"`.
    gravity: float
        Gravity constant [m s-2]. Required if `zcoord_type="pressure"`.

    Returns
    -------
    integ: xarray.DataArray
        Vertical integral.
    """
    # Do the vertical integration
    if coord_type == "height":
        # Integrate along the height coordinate
        if rho is None:
            # weight by air density
            raise ValueError("`rho` array is required to do weighting for 'height'-coordinate")
        # if isinstance(coord, collections.abc.Hashable):
        integ = integral(rho * xr_da, dim=dim, coord=coord)
        # integ /= integral(rho, dim=dim, coord=coord)
    elif coord_type == "pressure":
        # Integrate along the pressure coordinate
        if gravity is None:
            raise ValueError("`gravity` is required to do weighting for 'pressure'-coordinate")
        integ = -integral(xr_da, dim=dim, coord=coord) / gravity
    return integ


def meridional_mean(xr_da, lat_name="latitude"):
    """
    Calculate a meridional average of an `xarray.DataArray`.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data array with a latitude coordinate.
    lat_name: str, optional
        Name of y-coordinate

    Returns
    -------
    xarray.DataArray
        Array averaged over the latitudes.
    """
    coslat = da.cos(da.deg2rad(xr_da[lat_name]))
    xr_da_mean = xr_da.weighted(coslat).mean(dim=lat_name)
    # xr_da_mean = (xr_da * coslat).sum(dim=lat_name) / (coslat.sum(lat_name))
    return xr_da_mean


def moist_static_energy(
    cmpnt="all",
    temp=None,
    alt=None,
    spec_hum=None,
    c_p=None,
    gravity=None,
    latent_heat=None,
):
    """
    Calculate moist static energy or its components.

    .. math::
        MSE = DSE + LSE = (c_p T + g z) + L_v q

    Parameters
    ----------
    cmpnt: str, optional
        Component of MSE to output: "dry" | "latent" | "moist"
        By default, outputs all three of them: DSE, LSE, and their sum, MSE.
    temp: xarray.DataArray, optional
        Array of temperature [K].
    alt: xarray.DataArray, optional
        Array of level heights [m].
    spec_hum: xarray.DataArray, optional
        Array of specific humidity [kg kg-1].
    c_p: float, optional
        Dry air specific heat capacity [m2 s-2 K-1].
    gravity: float, optional
        Gravity constant [m s-2].
    latent_heat: float, optional
        Latent heat of vaporization [J kg-1].

    Returns
    -------
    xarray.Dataset
        Data arrays of moist static energy or its components.
    """
    if cmpnt in ["dry", "moist", "all"]:
        # Geopotential height
        ghgt = gravity * alt
        # Dry component: c_p T + g z
        dse = c_p * temp + ghgt
        if cmpnt == "dry":
            return xr.Dataset({"dry": dse})
    if cmpnt in ["latent", "moist", "all"]:
        # latent component :
        lse = latent_heat * spec_hum
        if cmpnt == "latent":
            return xr.Dataset({"latent": lse})
    if cmpnt in ["moist", "all"]:
        # dry and latent components
        mse = dse + lse
        if cmpnt == "moist":
            return xr.Dataset({"moist": mse})
        elif cmpnt == "all":
            return xr.Dataset({"dry": dse, "latent": lse, "moist": mse})


def vert_mer_mean_of_mse_flux(
    u,
    v,
    temp=None,
    alt=None,
    spec_hum=None,
    zcoord=None,
    rho=None,
    zcoord_type="height",
    lon_name="longitude",
    lat_name="latitude",
    z_name="level_height",
    cmpnt="all",
    opt="finite_diff",
    truncation=None,
    skiprows=None,
    c_p=1005,
    gravity=9.80665,
    latent_heat=2_501_000,
    r_planet=6_371_200,
):
    """
    Vertical and meridional integral of DSE, LSE and MSE fluxes.

    Wrapper-function to calculate the horizontal divergence of the dry static energy,
    latent static energy and moist static energy fluxes
    integrated over latitudes and in the vertical.

    Parameters
    ----------
    u: xarray.DataArray
        Array of zonal wind component [m s-1].
    v: xarray.DataArray, optional
        Array of meridional wind component [m s-1].
    temp: xarray.DataArray, optional
        Array of temperature [K].
    alt: xarray.DataArray, optional
        Array of model level heights (to calculate geopotential) [m].
    spec_hum: xarray.DataArray, optional
        Array of specific humidity [kg kg-1].
    zcoord: xarray.DataArray, optional
        Array of a coordinate to use for vertical integration.
    rho: xarray.DataArray, optional
        Array of air density [kg m-3]. Required if `zcoord_type="height"`.
    zcoord_type: str, optional
        Type of vertical coordinate ("height" or "pressure").
    lon_name: str, optional
        Name of x-coordinate.
    lat_name: str, optional
        Name of y-coordinate.
    z_name: str, optional
        Name of z-coordinate.
    cmpnt: str, optional
        Component of MSE to output: "dry" | "latent" | "moist"
        By default, outputs all three of them: DSE, LSE, and their sum, MSE.
    opt: str, optional
        Choose how to calculate the horizontal divergence: "spectral" | "finite_diff"
        spectral - use `windspharm` (with specified truncation)
        finite_diff - use finite differences
    truncation: int, optional
        Spectral truncation parameter passed to `windspharm`.
    skiprows: int, optional
        Omit this number of latitude points close to each pole to avoid spurious values.
    c_p: float
        Dry air specific heat capacity [m2 s-2 K-1].
    gravity: float
        Gravity constant [m s-2].
    latent_heat: float
        Latent heat of vaporization [J kg-1].
    r_planet: float, optional
        Radius of the planet [m]. Default is Earth's radius.

    Returns
    -------
    xarray.Dataset
        Data arrays of the flux divergence of moist static energy or its components.

    See also
    --------
    hdiv, moist_static_energy
    """
    # Calculate DSE
    mse_cmpnts = moist_static_energy(
        cmpnt=cmpnt,
        temp=temp,
        alt=alt,
        spec_hum=spec_hum,
        c_p=c_p,
        gravity=gravity,
        latent_heat=latent_heat,
    )
    results = {}
    for key, mse_cmpnt in mse_cmpnts.items():
        # Calculate horizontal fluxes (zonal and meridional components)
        # and their horizontal divergence in spherical coordinates
        if zcoord_type == "height":
            flux_x = u * mse_cmpnt * rho
            flux_y = v * mse_cmpnt * rho
        elif zcoord_type == "pressure":
            flux_x = u * mse_cmpnt
            flux_y = v * mse_cmpnt
        if opt == "finite_diff":
            result = hdiv(
                flux_x,
                flux_y,
                lon_name=lon_name,
                lat_name=lat_name,
                r_planet=r_planet,
            )
        elif opt == "spectral":
            from windspharm.xarray import VectorWind  # noqa

            vec = VectorWind(flux_x, flux_y, rsphere=r_planet)
            result = vec.divergence(truncation=truncation)
        # Do the vertical integration
        # result = mass_weighted_vertical_integral(
        #     result,
        #     z_name,
        #     coord=zcoord,
        #     coord_type=zcoord_type,
        #     rho=rho,  # XXX
        #     gravity=gravity,
        # )
        if zcoord_type == "height":
            result = integral(result, dim=z_name, coord=zcoord)
        elif zcoord_type == "pressure":
            result = -integral(result, dim=z_name, coord=zcoord) / gravity
        # Do the meridional averaging
        if isinstance(skiprows, int):
            result = result.isel(**{lat_name: slice(skiprows, -skiprows)})
        result = meridional_mean(result, lat_name=lat_name)
        results[key] = result
    return xr.Dataset(results)


def nondim_rossby_deformation_radius(
    method,
    temp,
    press=None,
    r_planet=EARTH_RADIUS,
    period=86_400,
    gravity=9.81,
    mw_dryair=28.97 * 1e-3,
    mgas_constant=8.314462,
    lon_name="longitude",
    lat_name="latitude",
    lev_name="level_height",
    time_name="time",
):

    """
    Estimate the circulation regime via the non-dimensional Rossby radius of deformation.

    For details, see eq. (1) in https://iopscience.iop.org/article/10.3847/1538-4357/ab9a4b

    Parameters
    ----------
    method: str
        Method to calculate the rossby deformation radius.
        "isothermal": eq. (1) in https://iopscience.iop.org/article/10.3847/1538-4357/ab9a4b
        "stratified": eq. (2) in https://iopscience.iop.org/article/10.3847/1538-4357/ab9a4b
    temp : xarray.DataArray
        Temperature proxy, e.g. surface temperature [K].
    press : xarray.DataArray, optional
        Atmospheric pressure [Pa]. Required for method="stratified".
    r_planet: float, optional
        Radius of the planet [m]. Default is Earth's radius.
    period: float, optional
        Period of the rotation [s]. Default is Earth's rotation period.
    gravity: float
        Gravity constant [m s-2].
    mw_dryair : float, optional
        Mean molecular weight of dry air [kg mol-1].
    mgas_constant : float, optional
        Molecular gas constant [J kg-1 mol-1].
    lon_name: str, optional
        Name of x-coordinate.
    lat_name: str, optional
        Name of y-coordinate.
    lev_name: str, optional
        Name of z-coordinate.
    time_name: str, optional
        Name of t-coordinate.

    Returns
    -------
    ratio: float
        Rossby radius of deformation divided by the radius of the planet.
    """
    if method == "isothermal":
        temp_mean = spatial_mean(temp, lon_name=lon_name, lat_name=lat_name).mean(dim=time_name)
        rossby_def_rad = rossby_deformation_radius_isothermal(
            temp=temp_mean,
            period=period,
            r_planet=r_planet,
            gravity=gravity,
            mw_dryair=mw_dryair,
            mgas_constant=mgas_constant,
        )
    elif method == "stratified":
        rossby_def_rad = rossby_deformation_radius_stratified(
            temp=temp,
            press=press,
            period=period,
            r_planet=r_planet,
            gravity=gravity,
            mw_dryair=mw_dryair,
            mgas_constant=mgas_constant,
            lev_name=lev_name,
        )
        rossby_def_rad = spatial_mean(rossby_def_rad, lon_name=lon_name, lat_name=lat_name).mean(
            dim=[time_name, lev_name]
        )

    return rossby_def_rad / r_planet


def potential_temperature(
    temp,
    press,
    gas_constant=287.058,
    c_p=1039,
    p_ref=100_000,
):
    """
    Calculate potential temperature
    ----------
    temp : xarray.DataArray
        Atmospheric temperature [K].
    press : xarray.DataArray
        Atmospheric pressure [Pa].
    gas_constant : float, optional
        Specific gas constant [J kg-1 K-1].
    c_p: float, optional
        Dry air specific heat capacity [m2 s-2 K-1].
    p_ref: float, optional
        Standard reference pressure [Pa].

    Returns
    -------
    theta: xarray.DataArray
        Atmospheric potential temperature [K]
    """
    theta = temp * (p_ref / press) ** (gas_constant / c_p)
    theta = theta.rename("air_potential_temperature")
    theta.attrs.update({"units": "K"})
    return theta


def rossby_deformation_radius_isothermal(
    temp,
    period=86_400,
    r_planet=EARTH_RADIUS,
    gravity=9.81,
    mw_dryair=28.97 * 1e-3,
    mgas_constant=8.314462,
):
    """
    Calculate the Rossby radius of deformation.

    For details, see eq. (1) in https://iopscience.iop.org/article/10.3847/1538-4357/ab9a4b

    Parameters
    ----------
    temp : xarray.DataArray
        Estimate of temperature, e.g. mean surface temperature [K].
    r_planet: float, optional
        Radius of the planet [m]. Default is Earth's radius.
    period: float, optional
        Period of the rotation [s]. Default is Earth's radius.
    gravity: float
        Gravity constant [m s-2].
    mw_dryair : float, optional
        Mean molecular weight of dry air [kg mol-1].
    mgas_constant : float, optional
        Molecular gas constant [J kg-1 mol-1].

    Returns
    -------
    radius: xarray.DataArray
        Rossby radius of deformation of shape (time).
    """
    omega = 2 * np.pi / period
    beta = 2 * omega / r_planet

    sc_hgt = scale_height(
        temp,
        r_planet=r_planet,
        gravity=gravity,
        mw_dryair=mw_dryair,
        mgas_constant=mgas_constant,
    )
    radius = ((gravity * sc_hgt) ** 0.5 / (2 * beta)) ** 0.5
    return radius


def rossby_deformation_radius_stratified(
    temp,
    press,
    period=86_400,
    r_planet=EARTH_RADIUS,
    gravity=9.81,
    mw_dryair=28.97 * 1e-3,
    mgas_constant=8.314462,
    c_p=1039,
    p_ref=100_000,
    lev_name="level_height",
):
    """
    Calculate the Rossby radius of deformation.
    For details, see eq. (2) in https://iopscience.iop.org/article/10.3847/1538-4357/ab9a4b

    Parameters
    ----------
    temp : xarray.DataArray
        Estimate of temperature, e.g. mean surface temperature [K].
    press : xarray.DataArray
        Atmospheric pressure [Pa].
    r_planet: float, optional
        Radius of the planet [m]. Default is Earth's radius.
    period: float, optional
        Period of the rotation [s]. Default is Earth's radius.
    gravity: float
        Gravity constant [m s-2].
    mw_dryair : float, optional
        Mean molecular weight of dry air [kg mol-1].
    mgas_constant : float, optional
        Molecular gas constant [J kg-1 mol-1].
    c_p: float, optional
        Dry air specific heat capacity [m2 s-2 K-1].
    p_ref : float, optional
        Standard reference pressure [Pa].
    lev_name: str, optional
        Name of z-coordinate.

    Returns
    -------
    radius: xarray.DataArray
        Rossby radius of deformation.
    """
    omega = 2 * np.pi / period
    beta = 2 * omega / r_planet

    sc_hgt = scale_height(
        temp,
        r_planet=r_planet,
        gravity=gravity,
        mw_dryair=mw_dryair,
        mgas_constant=mgas_constant,
    )
    bv_freq = brunt_vaisala_frequency(
        temp,
        press,
        gas_constant=mgas_constant / mw_dryair,
        c_p=c_p,
        p_ref=p_ref,
        gravity=gravity,
        lev_name=lev_name,
    )
    radius = ((bv_freq * sc_hgt) / (2 * beta)) ** 0.5
    return radius


def spatial_mean(xr_da, lon_name="longitude", lat_name="latitude"):
    """
    Perform averaging on an `xarray.DataArray` with latitude weighting.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data to average
    lon_name: str, optional
        Name of x-coordinate
    lat_name: str, optional
        Name of y-coordinate

    Returns
    -------
    xarray.DataArray
        Spatially averaged xarray.DataArray.
    """
    weights = da.cos(da.deg2rad(xr_da[lat_name]))
    res = xr_da.weighted(weights).mean(dim=[lon_name, lat_name])
    return res


def scale_height(
    temp,
    r_planet=EARTH_RADIUS,
    gravity=9.81,
    mw_dryair=28.97 * 1e-3,
    mgas_constant=8.314462,
):
    r"""
    Calculate the atmopheric scale height.

    .. math::
        H = R * T / (M * g)

    Parameters
    ----------
    temp : xarray.DataArray
        Estimate of temperature, e.g. mean surface temperature [K].
    r_planet: float, optional
        Radius of the planet [m]. Default is Earth's radius.
    gravity: float
        Gravity constant [m s-2].
    mw_dryair : float, optional
        Mean molecular weight of dry air [kg mol-1].
    mgas_constant : float, optional
        Molecular gas constant [J kg-1 mol-1].

    Returns
    -------
    xarray.DataArray
        Atmospheric scale height [m] of the same shape as `temp`.
    """
    return temp * mgas_constant / (mw_dryair * gravity)


def sfc_temp(ds, model_key, const):
    """Extract surface temperature from a THAI dataset."""
    model_names = names[model_key]
    out = ds[model_names.t_sfc].copy()
    if model_key == "ROCKE3D":
        out += const.t_melt  # convert from degC to K
    return out


def sfc_dn_lw_flux(ds, model_key, const):
    """Calculate the downward longwave radiative flux at the surface."""
    net_flux = sfc_net_up_lw_flux(ds, model_key, const)
    up_flux = const.stefan_boltzmann * sfc_temp(ds, model_key, const) ** 4
    return up_flux - net_flux


def sfc_net_up_lw_flux(ds, model_key, const):
    r"""
    Calculate the net downward longwave radiative flux at the surface.

    Parameters
    ----------
    ds: xarray.Dataset
        Input dataset containing relevant variables.
    model_key: str,
        Model name.

    Returns
    -------
    xarray.DataArray
    """
    if model_key == "ExoCAM":
        out = ds[names[model_key].sfc_net_down_lw]
    elif model_key == "LMDG":
        out = -ds[names[model_key].sfc_net_down_lw]
    elif model_key == "ROCKE3D":
        out = -(ds[names[model_key].sfc_dn_lw] - ds[names[model_key].sfc_up_lw])
    elif model_key == "UM":
        out = -ds[names[model_key].sfc_net_down_lw]
    return out


def spatial_sum(xr_da, lon_name="longitude", lat_name="latitude", r_planet=EARTH_RADIUS):
    """
    Calculate spatial integral of xarray.DataArray with grid cell weighting.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data to average
    lon_name: str, optional
        Name of x-coordinate
    lat_name: str, optional
        Name of y-coordinate
    r_planet: float
        Radius of the planet [metres], currently assumed spherical (not important anyway)

    Returns
    -------
    xarray.DataArray
        Spatially averaged xarray.DataArray.
    """
    lon = xr_da[lon_name].values
    lat = xr_da[lat_name].values

    area_weights = grid_cell_areas(lon, lat, r_planet=r_planet)

    return (xr_da * area_weights).sum(dim=[lon_name, lat_name])


def specific_humidity(ds, model_key):
    """Extract specific humidity from a THAI dataset."""
    names[model_key]
    if model_key == "LMDG":
        # LMD-G outputs mixing ratio instead of specific humidity,
        # so here MR is converted to SH.
        out = 1 / (ds.sh + 1)
    else:
        out = ds.sh
    return out


def time_mean(xr_da, time_name="time"):
    """
    Calculate a time average of an `xarray.DataArray`.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data array with a time coordinate.
    time_name: str, optional
        Name of t-coordinate

    Returns
    -------
    xarray.DataArray
        Array averaged over the time dimension.
    """
    xr_da_mean = xr_da.mean(dim=time_name)
    xr_da_mean.attrs.update(xr_da.attrs)
    return xr_da_mean


def toa_olr(ds, model_key, case, const):
    """Extract top-of-the-atmosphere outgoing longwave radiation from a THAI dataset."""
    model_names = names[model_key]
    if model_key == "ExoCAM":
        out = ds[model_names.toa_net_lw]
    elif model_key == "LMDG":
        out = ds[model_names.toa_olr]
    elif model_key == "ROCKE3D":
        out = ds[model_names.toa_olr_cs] - ds[model_names.toa_crf_lw]
    elif model_key == "UM":
        out = ds[model_names.toa_olr]
    return out


def upper_atm_vap_mean(ds, model_key, const, pres_levels=[100]):
    """Estimate the mean mixing ratio of water vapor in the upper atmosphere."""
    model_names = names[model_key]
    spec_hum = ds[model_names.sh].to_iris()
    if model_key == "LMDG":
        mix_ratio = ds[model_names.sh].to_iris()
    else:
        spec_hum = ds[model_names.sh].to_iris()
        mix_ratio = spec_hum / (1 - spec_hum)
    mix_ratio *= const.mw_ratio
    pres = ds[model_names.pres].to_iris()
    pres.convert_units("Pa")
    mix_ratio_on_pres_lev = stratify.relevel(
        mix_ratio,
        pres,
        pres_levels,
        axis=mix_ratio.coords(axis="z")[0].name(),
        interpolator=INTERPOLATOR,
    )
    mix_ratio_on_pres_lev.coord(pres.name()).attributes = {}
    mix_ratio_on_pres_lev = iris.util.squeeze(mix_ratio_on_pres_lev)
    return xr.DataArray.from_iris(mix_ratio_on_pres_lev)


def wind_rot_div(u, v, truncation=None, const=None):
    """Split the wind field into divergent and zonal mean and eddy rotational components."""
    from windspharm.xarray import VectorWind

    vec = VectorWind(u, v, rsphere=const.rplanet_m)
    div_cmpnt_u, div_cmpnt_v, rot_cmpnt_u, rot_cmpnt_v = vec.helmholtz(truncation=truncation)
    out = {}
    out["u_total"] = u
    out["v_total"] = v
    out["u_div"] = div_cmpnt_u.rename("irrotational_eastward_wind")
    out["v_div"] = div_cmpnt_v.rename("irrotational_northward_wind")
    out["u_rot"] = rot_cmpnt_u.rename("non_divergent_eastward_wind")
    out["v_rot"] = rot_cmpnt_v.rename("non_divergent_northward_wind")

    out["u_rot_zm"] = xr.broadcast(zonal_mean(rot_cmpnt_u), u)[0].rename(
        "zonal_mean_of_non_divergent_eastward_wind"
    )
    out["v_rot_zm"] = xr.broadcast(zonal_mean(rot_cmpnt_v), v)[0].rename(
        "zonal_mean_of_non_divergent_northward_wind"
    )
    out["u_rot_eddy"] = (rot_cmpnt_u - out["u_rot_zm"]).rename(
        "zonal_deviation_of_non_divergent_eastward_wind"
    )
    out["v_rot_eddy"] = (rot_cmpnt_v - out["v_rot_zm"]).rename(
        "zonal_deviation_of_non_divergent_northward_wind"
    )
    return xr.Dataset(out)


def zonal_mass_streamfunction(
    u,
    press,
    time_name="time",
    lev_name="level_height",
    lat_name="latitude",
    lon_name="longitude",
    g_planet=9.12,
    r_planet=5988740.0,
):
    r"""
    Calculate mean zonal mass streamfunction.

    .. math::
        \Psi_Z = \frac{2\pi a}{g} \int_{p_{sfc}}^{p_{top}} \overline{u}^* dp

    Parameters
    ----------
    u: xarray.DataArray
        Zonal wind [m s-1].
    press : xarray.DataArray
        Atmospheric pressure [Pa].
    lon_name: str, optional
        Name of x-coordinate.
    lat_name: str, optional
        Name of y-coordinate.
    lev_name: str, optional
        Name of z-coordinate.
    time_name: str, optional
        Name of t-coordinate.
    g_planet: float, optional
        Gravity constant [m s-2].
    r_planet: float, optional
        Radius of the planet [m].

    References
    ----------
    * Haqq-Misra & Kopparapu (2015), eq. 5;
    * Hartmann (1994), Global Physical Climatology, eq. 6.21

    Examples
    --------
    >>> mzsf = zonal_mass_streamfunction(
        u_lmdg,
        press_lmdg,
        time_name="time",
        lev_name="lev",
        lat_name="lat",
        lon_name="lon",
    )
    """
    const = 2 * np.pi * r_planet / g_planet
    # calulate the zonal average value minus the time-average component of zonal wind
    u_tmean = u.mean(dim=time_name)
    u_lattmean = zonal_mean(u_tmean, lon_name=lon_name)
    u_mean_star = u_lattmean - u_tmean

    deltap = press.mean(dim=time_name).diff(dim=lev_name)
    walker = u_mean_star * deltap.interp(**{lev_name: u_mean_star[lev_name]})
    walker = reverse_along_dim(walker, lev_name)
    walker = meridional_mean(walker.cumsum(dim=lev_name), lat_name=lat_name)
    return const * reverse_along_dim(walker, lev_name)


def zonal_mean(xr_da, lon_name="longitude"):
    """
    Calculate a zonal average of an `xarray.DataArray`.

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data array with a longitude coordinate.
    lon_name: str, optional
        Name of x-coordinate

    Returns
    -------
    xarray.DataArray
        Array averaged over the longitudes.
    """
    xr_da_mean = xr_da.mean(dim=lon_name)
    return xr_da_mean


def global_mean(arr, model_names):
    """Find a day-side average of THAI data array."""
    return spatial_mean(arr, model_names.x, model_names.y)


def dayside_mean(arr, model_names):
    """Find a day-side average of THAI data array."""
    # Assume the longitude span is -180 to +180 and SS point is at 0deg
    out = arr.where(abs(arr[model_names.x]) < 90)
    return spatial_mean(out, model_names.x, model_names.y)


def nightside_mean(arr, model_names):
    """Find a night-side average of THAI data array."""
    # Assume the longitude span is -180 to +180 and SS point is at 0deg
    out = arr.where(abs(arr[model_names.x]) > 90)
    return spatial_mean(out, model_names.x, model_names.y)


def terminator_mean(arr, model_names):
    """Find terminators average of THAI data array."""
    # Assume the longitude span is -180 to +180
    out = arr.interp(**{model_names.x: [-90, 90]})
    return spatial_mean(out, model_names.x, model_names.y)


def get_time_rel_days(da_time):
    """Convert points of a time coordinate to relative number of days."""
    vals = da_time.values.astype(float)
    rel_vals = vals - vals[0]
    units = da_time.units
    if isinstance(units, xr.DataArray):
        units = str(da_time.units.values)
    if units.lower().startswith("day"):
        return rel_vals
    elif units.lower().startswith("hour"):
        indicator = "h"
    elif units.lower().startswith("second"):
        indicator = "s"
    rel_vals = rel_vals.astype(f"timedelta64[{indicator}]")
    # convert to days
    return rel_vals / np.timedelta64(1, "D")
