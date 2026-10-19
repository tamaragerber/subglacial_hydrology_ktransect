#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data loading utilities for DEM, melt, and surface velocity fields.
"""

import pandas as pd
import numpy as np
import xarray as xr
import rioxarray as rxr


# ============================================================
# DEM
# ============================================================

def load_dem(fn):
    dem = rxr.open_rasterio(fn)
    return dem.squeeze()


# ============================================================
# MELT
# ============================================================

def load_melt(fn):

    def index_to_dt(df, add_2000=False):
        offset = 2000 if add_2000 else 0

        df.index = pd.to_datetime(
            ((offset + df["Year"]) * 1000 + df["Day"]) * 100 + df["Hour"],
            format="%Y%j%H",
        )
        return df

    melt = pd.read_csv(fn, sep=r"\s+")
    melt = index_to_dt(melt)
    melt = melt.drop(columns=["Year", "Day", "Hour"])

    # Convert to metres water equivalent
    melt *= 0.001

    # Convert elevation columns to float
    cols_as_float = [float(c[:-1]) for c in melt.columns]
    melt = melt.rename(columns=dict(zip(melt.columns, cols_as_float)))

    return melt


# ============================================================
# SURFACE VELOCITY
# ============================================================

def load_surface_velocity(fn):
    ds = xr.open_dataset(fn)

    v_surf = ds["Band1"].values
    x = ds["x"].values
    y = ds["y"].values

    xx, yy = np.meshgrid(x, y)

    return v_surf, x, y, xx, yy