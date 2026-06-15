#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

YEAR = 2012

STATIONS = {
    "lev4": (-175032.199, -2503893.928),
    "lev5": (-149294.207, -2502594.744),
    "lev6": (-111456.854, -2503397.163),
}

SECONDS_PER_YEAR = 365 * 24 * 3600

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------


def load_gnss_csv(path):
    """Load GNSS velocity file and ensure first column is datetime."""
    df = pd.read_csv(path, parse_dates=[0])
    df.iloc[:, 0] = pd.to_datetime(df.iloc[:, 0], format="%Y-%m-%d")
    return df


def extract_year_doy(df, date_col, year):
    """Return DOY array and dataframe for a given year."""
    df = df.copy()
    df["DOY"] = df[date_col].dt.dayofyear

    df_year = df[df[date_col].dt.year == year]

    return df_year["DOY"], df_year


def get_gnss_data(year):
    """Load GNSS data appropriate for the selected year."""

    if year in {2009, 2010, 2011, 2012}:
        files = {
            "lev4": "../data/gnss_velocities/gps_2009_2013_v3/level2/lev4/"
                    "lev4_2009_127_2013_131_disp_velocity_24h_legacy.csv",
            "lev5": "../data/gnss_velocities/gps_2009_2013_v3/level2/lev5/"
                    "lev5_2009_127_2013_131_disp_velocity_24h_legacy.csv",
            "lev6": "../data/gnss_velocities/gps_2009_2013_v3/level2/lev6/"
                    "lev6_2009_127_2013_131_disp_velocity_5d_legacy_tg.csv",
        }

    elif year in {2021, 2022, 2023}:
        files = {
            "lev4": "../data/its_live_velocities/lev4_2021-2025_24h.csv",
            "lev5": "../data/gnss_velocities/hislide-level2/lev5/"
                    "lev5_2021_129_2024_122_disp_velocity_24h_epochs.csv",
            "lev6": "../data/gnss_velocities/hislide-level2/lev6/"
                    "lev6_2021_127_2024_121_disp_velocity_24h_epochs.csv",
        }

    else:
        raise ValueError(f"No GNSS data configuration available for {year}")

    return {station: load_gnss_csv(path) for station, path in files.items()}


# -----------------------------------------------------------------------------
# MODEL OUTPUT FILES
# -----------------------------------------------------------------------------

combined_dir = f"{YEAR}_s/nolake_run_{YEAR}"
ntimesteps = 365

files = [
    os.path.join(combined_dir, f"forcing_{YEAR}_{day}.nc")
    for day in range(1, ntimesteps + 1)
]

# -----------------------------------------------------------------------------
# LOAD GNSS DATA
# -----------------------------------------------------------------------------

gnss_raw = get_gnss_data(YEAR)

# -----------------------------------------------------------------------------
# GRID AND STATION INDICES
# -----------------------------------------------------------------------------

with Dataset(files[0]) as ds:
    x = ds.variables["x"][:]
    y = ds.variables["y"][:]

station_idx = {}

for name, (xs, ys) in STATIONS.items():
    ix = np.argmin(np.abs(x - xs))
    iy = np.argmin(np.abs(y - ys))
    station_idx[name] = (iy, ix)

# -----------------------------------------------------------------------------
# EXTRACT MODEL TIMESERIES
# -----------------------------------------------------------------------------

time = np.arange(ntimesteps)

v_basal_ts = {station: [] for station in STATIONS}
water_ts = {station: [] for station in STATIONS}

for i, file in enumerate(files):

    with Dataset(file) as ds:
        v_basal = ds.variables["v_basal"][:]
        water = ds.variables["water_input_m_sec"][:]

    for station, (iy, ix) in station_idx.items():

        # convert m/s → m/yr 
        v_basal_ts[station].append(
            v_basal[iy, ix] * SECONDS_PER_YEAR
        )

        water_ts[station].append(
            water[iy, ix]
        )

# -----------------------------------------------------------------------------
# PREPARE GNSS TIMESERIES
# -----------------------------------------------------------------------------

gnss = {}

for station, df in gnss_raw.items():
    gnss[station] = extract_year_doy(df, "date", YEAR)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------

fig, axes = plt.subplots(
    len(STATIONS),
    1,
    figsize=(10, 8),
    sharex=True,
)

for ax, station in zip(axes, STATIONS):

    doy, df = gnss[station]

    ax.plot(
        time,
        v_basal_ts[station],
        label="model basal velocity",
    )

    ax.step(
        doy,
        df["v_myr"],
        label="GNSS velocity",
    )

    ax.set_title(station)
    ax.set_ylabel("Velocity (m yr$^{-1}$)")

    ax2 = ax.twinx()
    ax2.plot(
        time,
        water_ts[station],
        "--",
        label="Water input",
    )
    ax2.set_ylabel("Water input")

ax.set_xlabel("Day of year")
axes[0].legend(loc="upper left")

plt.tight_layout()
plt.show()