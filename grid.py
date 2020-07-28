# -*- coding: utf-8 -*-
"""Functionality related to cartesian geographical grids of xarray Data Arrays."""
import dask.array as da

import numpy as np

import xarray as xr


EARTH_RADIUS = 6371000.0  # m


__all__ = (
    "add_cyclic_point_to_da",
    "grid_cell_areas",
    "meridional_mean",
    "reverse_along_dim",
    "roll_da_to_0360",
    "roll_da_to_pm180",
    "spatial_integral",
    "spatial_mean",
    "wrap_lons",
    "zonal_mean",
)


def _guess_bounds(points, bound_position=0.5):
    """
    Guess bounds of grid cells.

    Simplified function from iris.coord.Coord.

    Parameters
    ----------
    points: numpy.array
        Array of grid points of shape (N,).
    bound_position: float, optional
        Bounds offset relative to the grid cell centre.

    Returns
    -------
    numpy.array
        Array of shape (N, 2).
    """
    diffs = np.diff(points)
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])

    min_bounds = points - diffs[:-1] * bound_position
    max_bounds = points + diffs[1:] * (1 - bound_position)

    return np.array([min_bounds, max_bounds]).transpose()


def _quadrant_area(radian_lat_bounds, radian_lon_bounds, r_planet):
    r"""
    Calculate spherical segment areas.

    Taken from the scitools-iris library.

    Area weights are calculated for each lat/lon cell as:
        .. math::
            r^2 (lon_1 - lon_0) ( sin(lat_1) - sin(lat_0))

    The resulting array will have a shape of
    *(radian_lat_bounds.shape[0], radian_lon_bounds.shape[0])*
    The calculations are done at 64 bit precision and the returned array
    will be of type numpy.float64.

    Parameters
    ----------
    radian_lat_bounds: numpy.array
        Array of latitude bounds (radians) of shape (M, 2).
    radian_lon_bounds: numpy.array
        Array of longitude bounds (radians) of shape (N, 2).
    r_planet: float
        Radius of the planet (currently assumed spherical).

    Returns
    -------
    numpy.array
       Array of grid cell areas of shape (M, N).
    """
    # ensure pairs of bounds
    if (
        radian_lat_bounds.shape[-1] != 2
        or radian_lon_bounds.shape[-1] != 2
        or radian_lat_bounds.ndim != 2
        or radian_lon_bounds.ndim != 2
    ):
        raise ValueError("Bounds must be [n,2] array")

    # fill in a new array of areas
    radius_sqr = r_planet ** 2
    radian_lat_64 = radian_lat_bounds.astype(np.float64)
    radian_lon_64 = radian_lon_bounds.astype(np.float64)

    ylen = np.sin(radian_lat_64[:, 1]) - np.sin(radian_lat_64[:, 0])
    xlen = radian_lon_64[:, 1] - radian_lon_64[:, 0]
    areas = radius_sqr * np.outer(ylen, xlen)

    # we use abs because backwards bounds (min > max) give negative areas.
    return np.abs(areas)


def add_cyclic_point_to_da(xr_da, coord_name):
    """
    Add a cyclic point to an `xarray.DataArray` and a corresponding
    coordinate. Mimics `cartopy.util.add_cyclic_point()`.

    Parameters
    ----------
    xr_da: xarray.DataArray
        An n-dimensional `DataArray` of data to add a cyclic point to.
    coord: str
        Coordinate name correspoinding to a 1-dimensional coordinate of `xr_da`.
        The coordinate values must be regularly spaced.

    Returns
    -------
    cyclic_xr_da: xarray.DataArray
        The array with a cyclic point added.
    """
    coord = xr_da[coord_name]
    axis = xr_da.get_axis_num(coord_name)
    if coord is not None:
        if coord.ndim != 1:
            raise ValueError("The coordinate must be 1-dimensional.")
        delta_coord = np.diff(coord)
        if not da.allclose(delta_coord, delta_coord[0]):
            raise ValueError("The coordinate must be equally spaced.")
        new_coord = xr.concat((coord, coord[-1:] + delta_coord[0]), dim=coord_name)

    slicer = [slice(None)] * xr_da.ndim
    slicer[axis] = slice(0, 1)
    cyclic_xr_da = xr.concat((xr_da, xr_da[tuple(slicer)]), dim=coord_name)
    cyclic_xr_da[coord_name] = new_coord
    return cyclic_xr_da


def grid_cell_areas(lon1d, lat1d, r_planet=EARTH_RADIUS):
    """
    Calculate grid cell areas given 1D arrays of longitudes and latitudes
    for a planet with the given radius.

    Parameters
    ----------
    lon1d: numpy.array
        Array of longitude points [degrees] of shape (M,)
    lat1d: numpy.array
        Array of latitude points [degrees] of shape (M,)
    r_planet: float, optional
        Radius of the planet [metres] (currently assumed spherical)

    Returns
    -------
    numpy.array
        Array of grid cell areas [metres**2] of shape (M, N).
    """
    lon_bounds_radian = np.deg2rad(_guess_bounds(lon1d))
    lat_bounds_radian = np.deg2rad(_guess_bounds(lat1d))
    area = _quadrant_area(lat_bounds_radian, lon_bounds_radian, r_planet)
    return area


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


def reverse_along_dim(xr_da, dim_name):
    """Reverse `xarray.DataArray` along one of its dimensions."""
    out = xr_da.reindex(**{dim_name: list(reversed(xr_da[dim_name]))})
    return out


def roll_da_to_0360(xr_da, lon_name="longitude"):
    """Roll DataArray's data with longitudes from (-180, 180) to (0, 360)."""
    # Roll longitudes and corresponding data
    out = xr_da.roll(**{lon_name: xr_da[lon_name].shape[0] // 2}, roll_coords=True)
    # Reset western (negative) longitudes to values within (180, 360) range
    out[lon_name] = out[lon_name] - (out[lon_name] // 360) * 360
    return out


def roll_da_to_pm180(xr_da, lon_name="longitude"):
    """Roll DataArray's data with longitudes from (0, 360) to (-180, 180)."""
    # Roll longitudes and corresponding data
    out = xr_da.roll(**{lon_name: xr_da[lon_name].shape[0] // 2}, roll_coords=True)
    # Reset western (negative) longitudes to values within (180, 360) range
    out[lon_name] = wrap_lons(out[lon_name], -180, 360)
    return out


def spatial_integral(
    xr_da, lon_name="longitude", lat_name="latitude", r_planet=EARTH_RADIUS
):
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


def wrap_lons(lons, base, period):
    """Wrap longitude values into the range between base and base+period."""
    # It is important to use 64bit floating precision
    lons = lons.astype(np.float64)
    return ((lons - base + period * 2) % period) + base


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
