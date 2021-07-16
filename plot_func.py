# -*- coding: utf-8 -*-
"""Various plotting-related functions."""
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import numpy as np

from aeolus.model import um
from aeolus.plot import (
    GeoAxesGrid,
    label_global_map_gridlines,
    subplot_label_generator,
    unit_format,
)
from calc import spatial_mean
from grid import add_cyclic_point_to_da

# from cartopy.util import add_cyclic_point


KW_CART = dict(transform=ccrs.PlateCarree())
KW_SBPLT_LABEL = dict(fontsize="xx-large", fontweight="bold", pad=5, loc="left")
KW_MAIN_TTL = dict(fontsize="xx-large", pad=5, loc="center")
KW_AUX_TTL = dict(fontsize="x-large", pad=5, loc="right")
# Axes grid specs
KW_AXGR = dict(
    axes_pad=(0.65, 0.5),
    cbar_location="right",
    cbar_mode="single",
    cbar_pad=0.2,
    cbar_size="1.5%",
    label_mode="",
)
KW_CBAR_TTL = dict(fontsize="small", pad=5)

# Locations of grid lines on maps
XLOCS = np.arange(-180, 181, 90)
YLOCS = np.arange(-90, 91, 30)


def make_map_figure(ncols, nrows, rect=111, **axgr_kw):
    """
    Make a figure with a grid of cartopy axes with the Robinson projection.

    Parameters
    ----------
    ncols: int
        Number of columns
    nrows: int
        Number of rows
    axgr_kw: dict, optional
        Parameters passed to `aeolus.plot.cart.GeoAxesGrid`.

    Returns
    -------
    matplotlib.figure.Figure, aeolus.plot.cart.GeoAxesGrid
        The figure and axes grid.
    """

    iletters = subplot_label_generator()

    fig = plt.figure(figsize=(8 * ncols, 4 * nrows))

    axgr = GeoAxesGrid(fig, rect, projection=ccrs.Robinson(), nrows_ncols=(nrows, ncols), **axgr_kw)
    for ax in axgr.axes_all:
        label_global_map_gridlines(
            fig, ax, XLOCS[1:-1], YLOCS[1:-1], degree=True, size="medium", xoff=-20
        )
        ax.gridlines(xlocs=XLOCS, ylocs=YLOCS, crs=ccrs.PlateCarree())
        ax.set_title(f"{next(iletters)}", **KW_SBPLT_LABEL)

    return fig, axgr


def draw_scalar(
    xr_arr,
    ax,
    method="contourf",
    cax=None,
    tex_units=None,
    cbar_ticks=None,
    use_cyclic=True,
    model_names=um,
    **plt_kw,
):
    """
    Plot a cube on a map.

    Parameters
    ----------
    xr_arr: xarray.DataArray
        Data array.
    ax: matplotlib.axes._subplots.AxesSubplot
        Cartopy axes.
    method: str, optional
        Method of plotting, e.g. "contour", "pcolormesh", etc.
    cax: matplotlib.axes._subplots.AxesSubplot or similar
        Axes for the colorbar.
    tex_units: str, optional
        TeX string of units to be attached to the colorbar.
    cbar_ticks: sequence, optional
        Colorbar ticks.
    use_cyclic: bool, optional
        Use `cartopy.utils.add_cyclic_point` for the data.
    model_names: base.Model, optional
        Container with model-specific names and coordinates.
    plt_kw: dict, optional
        Keywords for the plotting method.

    Returns
    -------
    Output of the plotting method.
    """
    if use_cyclic:
        xr_arr = add_cyclic_point_to_da(xr_arr, model_names.x)
    lon2d, lat2d = np.meshgrid(xr_arr[model_names.x], xr_arr[model_names.y])

    h = getattr(ax, method)(lon2d, lat2d, xr_arr, **plt_kw, **KW_CART)
    if cax is not None:
        cb = ax.figure.colorbar(h, cax=cax, aspect=80)
        if tex_units is not None:
            cb.ax.set_title(f"[{tex_units}]", **KW_CBAR_TTL)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
    return h


def draw_vector(
    u,
    v,
    ax,
    cax=None,
    tex_units=None,
    cbar_ticks=None,
    mag=(),
    xstride=1,
    ystride=1,
    model_names=um,
    qk_ref_wspd=None,
    kw_quiver={},
    kw_quiverkey={},
    quiverkey_xy=(0.17, 0.87),
):
    """
    Plot vectors of two cubes on a map.

    Parameters
    ----------
    u: xarray.DataArray
        X-component of the vector.
    v: xarray.DataArray
        Y-component of the vector.
    ax: matplotlib.axes._subplots.AxesSubplot
        Cartopy axes.
    cax: matplotlib.axes._subplots.AxesSubplot or similar
        Axes for the colorbar.
    tex_units: str, optional
        TeX string of units to be attached to the colorbar.
    cbar_ticks: sequence, optional
        Colorbar ticks.
    mag: tuple, optional
        Tuple of numpy arrays to be used for colour-coding the vectors.
    xstride: int, optional
        Stride x-component data.
    ystride: int, optional
        Stride y-component data.
    model_names: base.Model, optional
        Container with model-specific names and coordinates.
    qk_ref_wspd: float, optional
        Reference vector magnitude (wind speed).
        If given, a reference arrow (quiver key) is added to the figure.
    kw_quiver: dict, optional
        Keywords passed to quiver().
    kw_quiverkey: dict, optional
        Keywords passed to quiverkey().
    quiverkey_xy: tuple, optional
        Quiver key position.
    """
    lon2d, lat2d = np.meshgrid(u[model_names.x], u[model_names.y])
    skip = (slice(xstride, -xstride, xstride), slice(ystride, -ystride, ystride))

    h = ax.quiver(
        lon2d[skip],
        lat2d[skip],
        u[skip].values,
        v[skip].values,
        *[i[skip] for i in mag],
        **kw_quiver,
        **KW_CART,
    )
    if cax is not None and mag:
        cb = ax.figure.colorbar(h, cax=cax, aspect=80)
        if tex_units is not None:
            cb.ax.set_title(f"[{tex_units}]", **KW_CBAR_TTL)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)

    if qk_ref_wspd is not None:
        ax.quiverkey(
            h,
            *quiverkey_xy,
            qk_ref_wspd,
            fr"${qk_ref_wspd}$" + r" $m$ $s^{-1}$",
            **kw_quiverkey,
        )


def figsave(fig, imgname, **kw_savefig):
    """Save figure and print relative path to it."""
    save_dir = imgname.absolute().parent
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(imgname, **kw_savefig)
    pth = Path.cwd()
    rel_path = None
    pref = ""
    for par in pth.parents:
        pref += ".." + pth.anchor
        try:
            rel_path = f"{pref}{imgname.relative_to(par)}"
            break
        except ValueError:
            pass
    if rel_path is not None:
        print(f"Saved to {rel_path}.{plt.rcParams['savefig.format']}")


def linspace_pm1(n):
    """Return 2n evenly spaced numbers from -1 to 1, always skipping 0."""
    seq = np.linspace(0, 1, n + 1)
    return np.concatenate([-seq[1:][::-1], seq[1:]])


def darr_stats_string(
    darr, lon_name, lat_name, sep=" | ", eq_sign="=", fmt="auto", **kw_unit_format
):
    """Return min, mean and max of an `xarray.DataArray` as a string."""
    # Compute the stats
    _min = darr.min()
    # _mean = darr.mean()
    _mean = spatial_mean(darr, lon_name=lon_name, lat_name=lat_name)
    _max = darr.max()
    # Assemble a string
    txts = []
    for label, arr in zip(["min", "mean", "max"], [_min, _mean, _max]):
        if fmt == "auto":
            if (np.log10(abs(_mean)) < 0) or (np.log10(abs(_mean)) > 5):
                _str = f"{label}{eq_sign}{arr.values:.0e}"
            else:
                _str = f"{label}{eq_sign}{np.round(arr.values):.0f}"
        elif fmt == "pretty":
            _str = f"{label}{eq_sign}{unit_format(float(arr.values), **kw_unit_format)}"
        else:
            _str = f"{label}{eq_sign}{arr.values:{fmt}}"
        txts.append(_str)

    return sep.join(txts)


def set_alpha_in_cmap(cmap, alpha_min=0, alpha_max=1):
    """Set linearly spaced opacity channel in a matplotlib colormap."""
    cmap = plt.cm.get_cmap(cmap)

    # Get the colormap colors
    my_cmap = cmap(np.arange(cmap.N))

    # Set alpha
    my_cmap[:, -1] = np.linspace(alpha_min, alpha_max, cmap.N)

    # Create new colormap
    return mcol.ListedColormap(my_cmap)
