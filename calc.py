"""Science diagnostics."""
import dask.array as da

from grid import EARTH_RADIUS


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
    h_div = h_div.rename("horizontal_divergence")
    h_div.attrs = {
        "units": "s-1",
        "long_name": "horizontal_divergence"
    }
    return h_div
