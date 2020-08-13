"""Science diagnostics."""
import dask.array as da

import numpy as np

import xarray as xr

from grid import EARTH_RADIUS, meridional_mean, spatial_mean


__all__ = (
    "brunt_vaisala_frequency",
    "integral",
    "hdiv",
    "mass_weighted_vertical_integral",
    "moist_static_energy",
    "nondim_rossby_deformation_radius",
    "potential_temperature",
    "rossby_deformation_radius_isothermal",
    "rossby_deformation_radius_stratified",
    "scale_height",
    "vert_mer_mean_of_mse_flux",
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
    return ds._replace_with_new_dims(
        variables, coord_names=coord_names, indexes=indexes
    )


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
        Name of y-coordinate.

    Returns
    -------
    bv_freq: xarray.DataArray
        Brunt-Väisälä frequency [s-1]
    """
    # Compute potential temperature from real temperature and pressure
    theta = potential_temperature(
        temp, press, gas_constant=gas_constant, c_p=c_p, p_ref=p_ref,
    )

    bv_freq = ((gravity / theta) * theta.differentiate(dim=lev_name)) ** 0.5
    bv_freq = bv_freq.rename("air_potential_temperature")
    bv_freq.attrs.update({"units": "K"})
    return bv_freq


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
        return xr_da.integrate(dim=dim, datetime_unit=datetime_unit)
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


def hdiv(
    i_arr, j_arr, lon_name="longitude", lat_name="latitude", r_planet=EARTH_RADIUS
):
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
            raise ValueError(
                "`rho` array is required to do weighting for 'height'-coordinate"
            )
        # if isinstance(coord, collections.abc.Hashable):
        integ = integral(rho * xr_da, dim=dim, coord=coord)
        # integ /= integral(rho, dim=dim, coord=coord)
    elif coord_type == "pressure":
        # Integrate along the pressure coordinate
        if gravity is None:
            raise ValueError(
                "`gravity` is required to do weighting for 'pressure'-coordinate"
            )
        integ = -integral(xr_da, dim=dim, coord=coord) / gravity
    return integ


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
    list
        List of data arrays of moist static energy or its components.
    """
    if cmpnt in ["dry", "moist", "all"]:
        # Geopotential height
        ghgt = gravity * alt
        # Dry component: c_p T + g z
        dse = c_p * temp + ghgt
        if cmpnt == "dry":
            return [dse]
    if cmpnt in ["latent", "moist", "all"]:
        # latent component :
        lse = latent_heat * spec_hum
        if cmpnt == "latent":
            return [lse]
    if cmpnt in ["moist", "all"]:
        # dry and latent components
        mse = dse + lse
        if cmpnt == "moist":
            return [mse]
        elif cmpnt == "all":
            return dse, lse, mse


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
    list
        List of data arrays of the flux divergence of moist static energy or its components.

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
    results = []
    for mse_cmpnt in mse_cmpnts:
        # Calculate horizontal fluxes (zonal and meridional components)
        # and their horizontal divergence in spherical coordinates
        result = hdiv(
            u * mse_cmpnt,
            v * mse_cmpnt,
            lon_name=lon_name,
            lat_name=lat_name,
            r_planet=r_planet,
        )
        # Do the vertical integration
        result = mass_weighted_vertical_integral(
            result,
            z_name,
            coord=zcoord,
            coord_type=zcoord_type,
            rho=rho,
            gravity=gravity,
        )
        # Do the meridional averaging
        result = meridional_mean(result, lat_name=lat_name)
        results.append(result)
    return results


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
        temp_mean = spatial_mean(temp, lon_name=lon_name, lat_name=lat_name).mean(
            dim=time_name
        )
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
        )
        rossby_def_rad = spatial_mean(
            rossby_def_rad, lon_name=lon_name, lat_name=lat_name
        ).mean(dim=[time_name, lev_name])

    return rossby_def_rad / r_planet


def potential_temperature(
    temp, press, gas_constant=287.058, c_p=1039, p_ref=100_000,
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
    theta = theta.rename("brunt_vaisala_frequency")
    theta.attrs.update({"units": "s-1"})
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
