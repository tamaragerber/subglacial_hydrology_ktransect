#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue May 19 11:42:38 2026
@author: tamara
"""

import xarray as xr
import numpy as np
import pandas as pd


# ============================================================
# BASAL VELOCITY
# ============================================================

def compute_basal_velocity(doy, stations, v_surf, basal_sliding_fraction):
    """
    calculate basal sliding velocity from basal sliding fraction of gnss velocity, 
    then interpolate onto the grid for each timestep
    """
    vb_base = basal_sliding_fraction * v_surf.copy()
    vb = vb_base.copy()

    for st in stations.values():

        if np.isnan(st["gnss_series"][doy - 1]):
            continue

        target = st["gnss_series"][doy - 1]

        if "gnss_reference" not in st:
            valid = np.where(np.isfinite(st["gnss_series"]))[0]
            st["gnss_reference"] = st["gnss_series"][valid[0]]

        dv_needed = target - st["gnss_reference"]
        vb += dv_needed * st["W"]

    return vb


# ============================================================
# INTERPOLATION
# ============================================================

def interpolate_vsurf_to_dem_grid(v_surf, x, y, dem):
    """
    Interpolate surface velocity onto DEM grid of model domain.
    """

    da = xr.DataArray(
        v_surf,
        coords={"y": y, "x": x}
    )

    da_interp = da.interp(
        x=dem.x,
        y=dem.y,
        method="linear"
    )

    x_new = dem.x.values
    y_new = dem.y.values
    xx, yy = np.meshgrid(x_new, y_new)

    return da_interp.values, x_new, y_new, xx, yy


# ============================================================
# GNSS PROCESSING
# ============================================================

def implement_gnss_velocities(config, xx, yy, x, y):

    stations = config["gnss_stations"]
    year = config["meta"]["year"]

    full_doy = np.arange(1, 366)

    for st in stations.values():

        df = pd.read_csv(st["gnss_file"])

        df["date"] = pd.to_datetime(df["date"])
        df = df[df["date"].dt.year == year]

        # remove leap day if present
        if len(df) > 365:
            df = df[~((df["date"].dt.month == 2) & (df["date"].dt.day == 29))]

        df["doy"] = df["date"].dt.dayofyear
        gnss_doy = df["doy"].values

        df["v_ms"] = df["v_myr"] / (365 * 24 * 60 * 60)

        # basal sliding from GNSS
        v_basal = df["v_ms"].values * config["basal_sliding"]["basal_sliding_fraction"]

        background = np.nan
        v_interp = np.full_like(full_doy, background, dtype=float)

        start = gnss_doy[1] - 1
        v_interp[start:] = v_basal[1:]

        v_smooth = pd.Series(v_interp).rolling(
            window=10,
            center=True,
            min_periods=5
        ).mean().values

        st["gnss_series"] = v_smooth

        decay_radius = 5 * st["ice_thickness"]

        dist = np.sqrt((xx - st["x"])**2 + (yy - st["y"])**2)
        st["W"] = np.exp(-dist / decay_radius)

        st["ix"] = np.argmin(np.abs(x - st["x"]))
        st["iy"] = np.argmin(np.abs(y - st["y"]))

    return stations


# ============================================================
# UTILITIES
# ============================================================

def get_column_fallback(df, candidates):
    """
    Return first available column among candidate names.
    """
    for c in candidates:
        if c in df.columns:
            return df[c].to_numpy(dtype=float)
    return None