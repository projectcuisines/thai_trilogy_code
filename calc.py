"""Science diagnostics."""
import dask.array as da

import xarray as xr

from grid import EARTH_RADIUS, meridional_mean


__all__ = (
    "integral",
    "hdiv",
    "mass_weighted_vertical_integral",
    "moist_static_energy",
    "vert_mer_mean_of_dse_flux",
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
    h_div = h_div.rename("horizontal_divergence")
    h_div.attrs = {"units": "s-1", "long_name": "horizontal_divergence"}
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


def moist_static_energy(temp, alt, spec_hum, c_p, gravity, latent_heat):
    """
    Calculate moist static energy and its components.

    .. math::
        MSE = DSE + LSE = (c_p T + g z) + L_v q

    Parameters
    ----------
    temp: xarray.DataArray
        Array of temperature [K].
    alt: xarray.DataArray
        Array of level heights [m].
    spec_hum: xarray.DataArray
        Array of specific humidity [kg kg-1].
    c_p: float
        Dry air specific heat capacity [m2 s-2 K-1].
    gravity: float
        Gravity constant [m s-2].
    latent_heat: float
        Latent heat of vaporization [J kg-1].

    Returns
    -------
    dse: xarray.DataArray
        Array of dry static energy.
    lse: xarray.DataArray
        Array of latent static energy.
    mse: xarray.DataArray
        Array of moist static energy (Sum of DSE and LSE).
    """
    # Geopotential height
    ghgt = gravity * alt
    # Dry component: c_p T + g z
    dse = c_p * temp + ghgt
    # latent component :
    lse = latent_heat * spec_hum
    # dry and latent components
    mse = dse + lse
    return dse, lse, mse


def vert_mer_mean_of_dse_flux(
    temp,
    alt,
    spec_hum,
    u,
    v,
    zcoord=None,
    rho=None,
    zcoord_type="height",
    lon_name="longitude",
    lat_name="latitude",
    z_name="level_height",
    c_p=1005,
    gravity=9.80665,
    latent_heat=2_501_000,
    r_planet=6_371_200,
):
    """
    Vertical and meridional integral of DSE flux.

    Wrapper-function to calculate the horizontal divergence of the dry static energy flux,
    integrated over latitudes and in the vertical.

    Parameters
    ----------
    temp: xarray.DataArray
        Array of temperature [K].
    alt: xarray.DataArray
        Array of model level heights (to calculate geopotential) [m].
    spec_hum: xarray.DataArray
        Array of specific humidity [kg kg-1].
    u: xarray.DataArray
        Array of zonal wind component [m s-1].
    v: xarray.DataArray
        Array of meridional wind component [m s-1].
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
    flux_div_zm_int: xarray.DataArray
        Array of vertical and meridional integral of DSE flux divergence.
    """
    # Calculate DSE
    dse, _, _ = moist_static_energy(
        temp=temp,
        alt=alt,
        spec_hum=spec_hum,
        c_p=c_p,
        gravity=gravity,
        latent_heat=latent_heat,
    )
    # Calculate horizontal fluxes of DSE (zonal and meridional components)
    # and their horizontal divergence in spherical coordinates
    flux_div = hdiv(
        u * dse, v * dse, lon_name=lon_name, lat_name=lat_name, r_planet=r_planet,
    )

    # Do the vertical integration
    flux_div_z_int = mass_weighted_vertical_integral(
        flux_div,
        z_name,
        coord=zcoord,
        coord_type=zcoord_type,
        rho=rho,
        gravity=gravity,
    )

    # Do the meridional averaging
    flux_div_zm_int = meridional_mean(flux_div_z_int, lat_name=lat_name)
    return flux_div_zm_int
