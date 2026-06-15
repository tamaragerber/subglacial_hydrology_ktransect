#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate daily subglacial hydrology forcing files for GlaDS / Elmer/Ice.

# This script converts elevation-dependent surface runoff / melt (from SEBM or
# similar models) into spatially distributed daily water input fields on an ice
# sheet DEM. It produces NetCDF forcing files suitable for subglacial hydrology
# and ice dynamics models (e.g. GlaDS / Elmer/Ice).

# The workflow includes:

# 1. INPUT DATA PROCESSING
#    - Reads hourly SEBM melt/runoff time series (mm w.e.)
#    - Converts to metres water equivalent
#    - Aggregates to daily totals
#    - Removes leap days to enforce 365-day year consistency
#    - Interpolates elevation-binned melt onto a 2D DEM using linear regression

# 2. ELEVATION-DEPENDENT HYDROLOGY PARAMETERISATION
#    - Builds spatial fields of:
#      a) cumulative melt threshold for hydrofracture initiation (Tm)
#      b) fraction of accumulated water that can drain vertically (Th)
#    - Both vary linearly between user-defined elevation bounds (defined in the config file).
#

# 3. SURFACE WATER STORAGE & HYDROFRACTURE LOGIC
#    - Surface melt is accumulated on each grid cell
#    - Water is retained until a cumulative threshold is reached
#    - Once threshold is exceeded:
#        → all stored water is released instantaneously to the sheet
#    - After first drainage event:
#        → melt is routed directly to the bed (no further storage)

# 4. ELEVATION-DEPENDENT MELT SCALING
#    - Optional correction applied to SEBM melt inputs
#    - Linear scaling between two elevation bounds (scale_lower–scale_upper)
#    - Used to compensate for systematic overestimation in high-melt years

# 5. BASAL SLIDING COUPLING (GNSS CONSTRAINED)
#    - GNSS stations are used to define time-varying surface velocity signals
#    - Signals are spatially spread using Gaussian weighting
#    - Converted into basal velocity perturbations added to background flow

# 6. OPTIONAL LAKE DRAINAGE ("BLISTER" FORCING)
#    - Optional synthetic lake drainage events can be injected
#    - Controlled via config:
#          lake_drainage:
#              drainage: true / false
#    - Each lake event injects a spatially distributed pressure/melt pulse
#    - Implemented using a parabolic "blister" function
#    - Allows representation of rapid supraglacial lake drainage to bed

# 7. OUTPUTS
#    - Daily NetCDF forcing fields:
#          water_input_m_sec
#          v_basal
#    - Optional diagnostic time series:
#          accumulation state
#          drainage state
#          total supply
#    - Optional diagnostic plots of elevation dependence and time series

# ===============================================================================
# NOTES
# ===============================================================================

# - Designed for a single melt season (no multi-year handling)
# - Requires consistent 365-day indexing
# - Elevation bin melt input must come from SEBM-style datasets
# - Hydrofracture thresholds should be calibrated per year using observations
#   (e.g. lake drainage timing from satellite imagery)

# ===============================================================================
# AUTHORS
# ===============================================================================

# Tamara Gerber, Andrew Tedstone
"""

import os
import shutil
from copy import deepcopy
from pathlib import Path

import argparse
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import statsmodels.api as sm
import matplotlib.pyplot as plt

from load_data import (
    load_dem,
    load_melt,
    load_surface_velocity,
)

from velocities import (
    compute_basal_velocity,
    interpolate_vsurf_to_dem_grid,
    implement_gnss_velocities,
)

from lake_drainage import (
    generate_blister_on_grid,
    parse_lake_drainage,
)

from surface_hydrology import (
    create_melt_threshold_raster,
    create_hydrofracture_mult_raster,
)

plt.close("all")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    args = argparse.Namespace(yaml="config_2023.yaml")

    with open(args.yaml) as f:
        config = yaml.safe_load(f.read())

    print("Running:", config["meta"]["name"])

    output_to = Path(config["outputs"]["path"]) / config["meta"]["name"]
    output_to.mkdir(parents=True, exist_ok=True)

    shutil.copy(args.yaml, output_to / args.yaml)
    print("Outputs:", output_to)

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------
    dem = load_dem(config["input_files"]["dem"])
    melt = load_melt(config["input_files"]["seb"])

    assert config["date_range"]["start"].year == config["date_range"]["finish"].year
    assert config["melt_thresholds"]["lower_elev"] < config["melt_thresholds"]["upper_elev"]

    v_surf, x, y, xx, yy = load_surface_velocity(
        config["input_files"]["surface_velocity"]
    )

    lake_cfg = config.get("lake_drainage", {})
    use_lakes = lake_cfg.get("drainage", False)

    # --------------------------------------------------------
    # Velocity preprocessing
    # --------------------------------------------------------
    v_surf, x_new, y_new, xx, yy = interpolate_vsurf_to_dem_grid(
        v_surf, x, y, dem
    )

    gnss_stations = implement_gnss_velocities(
        config, xx, yy, x, y
    )

    # --------------------------------------------------------
    # Surface hydrology setup
    # --------------------------------------------------------
    Tm = create_melt_threshold_raster(
        dem,
        config["melt_thresholds"]["lower_elev"],
        config["melt_thresholds"]["lower_melt"],
        config["melt_thresholds"]["upper_elev"],
        config["melt_thresholds"]["upper_melt"],
    )

    Th = create_hydrofracture_mult_raster(
        dem,
        config["melt_thresholds"]["lower_elev"],
        config["melt_thresholds"]["upper_elev"],
    )

    if config["outputs"]["diagnostic_plots"]:
        plt.figure()
        plt.imshow(Tm)
        plt.colorbar()
        plt.title("Melt threshold (m w.e.)")
        plt.savefig(output_to / "threshold_mwe.pdf")

        plt.figure()
        plt.imshow(Th)
        plt.colorbar()
        plt.title("Hydrofracture fraction")
        plt.savefig(output_to / "threshold_frac.pdf")

    # --------------------------------------------------------
    # Prepare melt input
    # --------------------------------------------------------
    melt = melt.loc[
        config["date_range"]["start"] : config["date_range"]["finish"]
    ]

    melt_resampled_sum = melt.resample("1D").sum()

    melt_resampled_sum = melt_resampled_sum[
        ~(
            (melt_resampled_sum.index.month == 2)
            & (melt_resampled_sum.index.day == 29)
        )
    ]

    full_index = pd.date_range(
        start=config["date_range"]["start"],
        end=config["date_range"]["finish"],
        freq="D",
    )

    full_index = full_index[
        ~((full_index.month == 2) & (full_index.day == 29))
    ]

    melt_resampled_sum = melt_resampled_sum.reindex(
        full_index, fill_value=0
    )

    X = sm.add_constant(melt_resampled_sum.columns)

    accum_sum = np.zeros(dem.shape)
    drained = np.zeros(dem.shape)

    lakes = parse_lake_drainage(config)

    # --------------------------------------------------------
    # Storage for diagnostics
    # --------------------------------------------------------
    store_accum_state = []
    store_drained = []
    store_supply_total = []

    # ========================================================
    # DAILY LOOP
    # ========================================================
    for i, (ix, data) in enumerate(melt_resampled_sum.iterrows()):

        doy = i + 1
        year = ix.year

        f = sm.OLS(data, X).fit()

        melt_sum_grid = (
            float(f.params["x1"]) * dem
            + float(f.params["const"])
        )

        melt_sum_grid = np.maximum(melt_sum_grid, 0)

        # ----------------------------------------------------
        # Melt scaling
        # ----------------------------------------------------
        scale_lower = 1300
        scale_upper = 1800
        scale_factor = config["melt_thresholds"]["melt_scaling"]

        slope = (scale_factor - 1) / (scale_upper - scale_lower)

        scale = np.where(
            (dem > scale_lower) & (dem < scale_upper),
            1 + slope * (dem - scale_lower),
            1.0,
        )

        scale = np.where(dem >= scale_upper, scale_factor, scale)

        melt_sum_grid *= scale

        # ----------------------------------------------------
        # Accumulation + hydrofracture
        # ----------------------------------------------------
        accum_sum += np.where(drained == False, melt_sum_grid, 0)

        hydro = np.where(accum_sum >= Tm, accum_sum, 0)
        accum_sum -= hydro

        hydro_ms = hydro * Th / (60 * 60 * 24)

        melt_ms = melt_sum_grid / (60 * 60 * 24)
        melt_ms = np.where(drained == True, melt_ms, 0)

        today_supply_total = hydro_ms + melt_ms + 3.17e-11

        # ----------------------------------------------------
        # Lake forcing
        # ----------------------------------------------------
        if use_lakes:
            for lake in lakes:
                if lake["doy"] == doy:

                    ix_lake = np.argmin(np.abs(x - lake["x"]))
                    iy_lake = np.argmin(np.abs(y - lake["y"]))
                    H_local = dem.values[iy_lake, ix_lake]

                    blister = generate_blister_on_grid(
                        lake["V"],
                        H_local,
                        lake["x"],
                        lake["y"],
                        xx,
                        yy,
                    )

                    today_supply_total += blister

        drained = np.where(
            (hydro_ms > 0) | (drained == True),
            True,
            False,
        )

        vb = compute_basal_velocity(
            doy,
            gnss_stations,
            v_surf,
            config["basal_sliding"]["basal_sliding_fraction"],
        )

        if config["outputs"]["diagnostic_netcdfs"] or config["outputs"]["diagnostic_plots"]:
            store_accum_state.append(deepcopy(accum_sum))
            store_supply_total.append(today_supply_total)
            store_drained.append(drained)

        if config["outputs"]["daily_netcdfs"]:
            ds = xr.Dataset(
                {
                    "water_input_m_sec": (("y", "x"), today_supply_total),
                    "v_basal": (("y", "x"), vb),
                },
                coords={"y": dem.y, "x": dem.x},
            )

            suffix = "_with_lakes" if use_lakes else ""
            fn = os.path.join(output_to, f"forcing{suffix}_{year}_{doy}.nc")
            ds.to_netcdf(fn)

    # ========================================================
    # DIAGNOSTICS
    # ========================================================
    if (
        config["outputs"]["diagnostic_netcdfs"]
        or config["outputs"]["diagnostic_plots"]
    ):

        accum_ts = xr.DataArray(
            store_accum_state,
            coords={"time": melt_resampled_sum.index, "y": dem.y, "x": dem.x},
        )

        supply_ts = xr.DataArray(
            store_supply_total,
            coords={"time": melt_resampled_sum.index, "y": dem.y, "x": dem.x},
        )

        drain_ts = xr.DataArray(
            store_drained,
            coords={"time": melt_resampled_sum.index, "y": dem.y, "x": dem.x},
        )

    if config["outputs"]["diagnostic_netcdfs"]:
        supply_ts.name = "water_input_m_sec"
        supply_ts.to_netcdf(output_to / "water_input_m_sec_timeseries.nc")

        drain_ts.name = "drainage_open_bool"
        drain_ts.to_netcdf(output_to / "drainage_timeseries.nc")

    print("Finished.")


