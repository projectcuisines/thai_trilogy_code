# -*- coding: utf-8 -*-
"""Functions to load and clean up THAI data."""
import dask
import numpy as np
import xarray as xr

from model_exocam import adjust_exocam_grid
from model_rocke3d import adjust_rocke3d_grid
from model_um import open_mf_um, prep_um_ds
import mypaths
from names import exocam, lmdg, rocke3d

THAI_DIR = mypaths.extdatadir


def load_exocam(case):
    """Load and cleanup ExoCAM data for THAI."""
    files = str(THAI_DIR / "ExoCAM" / f"{case}_ExoCAM*.nc")
    ds = xr.open_mfdataset(
        files,
        decode_times=False,
        combine="nested",
        concat_dim=exocam.t,
        chunks=LOAD_CONF["ExoCAM"]["chunks"],
    )
    _tmp_ds = {}
    for d in ds.data_vars:
        _tmp_ds[d] = adjust_exocam_grid(ds[d])
    ds = xr.Dataset(_tmp_ds)
    ds = ds.rename({"lat": exocam.y, "lon": exocam.x})
    return ds


def load_lmdg(case):
    """Load and cleanup LMDG data for THAI."""
    files = str(THAI_DIR / "LMDG" / f"{case}_LMDG*.nc")
    ds = xr.open_mfdataset(
        files,
        decode_times=False,
        combine="nested",
        concat_dim=lmdg.t,
        chunks=LOAD_CONF["LMDG"]["chunks"],
    )
    # Assign CF-compliant time units
    # Otherwise the time has no units and resets to 0 in each of the files.
    delta = np.diff(ds[lmdg.t].values)
    delta[delta < 0] = delta[0]
    new_vals = np.cumsum(np.hstack([0, delta])) + delta[0]
    ds = ds.assign_coords(**{lmdg.t: new_vals})
    ds[lmdg.t].attrs["units"] = "days since 2000-01-01 00:00:00"
    return ds


def load_rocke3d(case):
    """Load and cleanup ROCKE3D data for THAI."""
    files = str(THAI_DIR / "ROCKE3D" / f"{case}_ROCKE3D*.nc")
    ds = xr.open_mfdataset(
        files,
        decode_times=False,
        chunks=LOAD_CONF["ROCKE3D"]["chunks"],
    )
    _tmp_ds = {}
    for d in ds.data_vars:
        _tmp_ds[d] = adjust_rocke3d_grid(ds[d])
    ds = xr.Dataset(_tmp_ds)
    ds = ds.rename({"lat": rocke3d.y, "lon": rocke3d.x})
    return ds


def load_um(case):
    """Load and cleanup UM data for THAI."""
    files = sorted((THAI_DIR / "UM").glob(f"{case}_UM*.nc"))
    ds = open_mf_um(
        files,
        main_time="hourly",
        rad_time="hourly_rad",
        decode_times=False,
        chunks=LOAD_CONF["UM"]["chunks"],
    )
    with dask.config.set(**{"array.slicing.split_large_chunks": False}):
        ds = prep_um_ds(ds, vert_lev_miss_val="drop")
    return ds


LOAD_CONF = {
    "ExoCAM": {
        "chunks": {"time": 61},
        "loader": load_exocam,
    },
    "LMDG": {
        "chunks": {"Time": 61},
        "loader": load_lmdg,
    },
    "ROCKE3D": {
        "chunks": {"time": 50},  # ROCKE3D has 2500 time steps
        "loader": load_rocke3d,
    },
    "UM": {
        "chunks": {"hourly": 61},
        "loader": load_um,
    },
}
