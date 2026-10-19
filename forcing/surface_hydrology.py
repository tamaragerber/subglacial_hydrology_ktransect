
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:46:59 2026

@author: tamara
"""

import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd

def create_melt_threshold_raster(dem, th_elev_lower, th_melt_lower, th_elev_higher, th_melt_higher):
    Tm = np.zeros(dem.shape)
    # Set threshold value up to the lower elevation threshold
    Tm = np.where(dem < th_elev_lower, th_melt_lower, Tm)
    # Set threshold value above the upper elevation threshold
    Tm = np.where(dem >= th_elev_higher, th_melt_higher, Tm)
    # Set threshold between these two elevations
    fit = sm.OLS([th_melt_lower, th_melt_higher], sm.add_constant([th_elev_lower, th_elev_higher])).fit()
    Tm = np.where((dem >= th_elev_lower) & (dem < th_elev_higher), (fit.params[1] * dem + fit.params[0]), Tm) 
    return Tm

def create_hydrofracture_mult_raster(dem, th_elev_lower, th_elev_higher):
    # The 'allow hydrofacture drainage' switch/multiplier
    Th = np.zeros(dem.shape)
    # Set threshold value up to the lower elevation threshold
    Th = np.where(dem < th_elev_lower, 1, Th)
    # Set threshold value above the upper elevation threshold
    Th = np.where(dem >= th_elev_higher, 0, Th)
    # Set threshold between these two elevations
    fit = sm.OLS([1, 0], sm.add_constant([th_elev_lower, th_elev_higher])).fit()
    Th = np.where((dem >= th_elev_lower) & (dem < th_elev_higher), (fit.params[1] * dem + fit.params[0]), Th) 
    return Th


def prepare_melt_processing(dem, melt, config, output_to):
    """
    Prepare melt thresholds and daily melt data for interpolation.

    Parameters
    ----------
    dem : ndarray or xarray.DataArray
        DEM raster.
    melt : pandas.DataFrame
        Melt time series with datetime index.
    config : dict
        Configuration dictionary.
    output_to : pathlib.Path
        Output directory for plots.

    Returns
    -------
    Tm : ndarray
        Melt threshold raster.
    Th : ndarray
        Hydrofracture multiplier raster.
    melt_resampled_sum : pandas.DataFrame
        Daily melt sums.
    X : ndarray or DataFrame
        Predictor matrix with constant term added.
    accum_sum : ndarray
        Accumulated melt raster initialized to zeros.
    drained : ndarray
        Drainage-state raster initialized to zeros.
    """

    # ------------------------------------------------------------------
    # Calculate elevation-dependent threshold maps
    Tm = create_melt_threshold_raster(
        dem,
        config['melt_thresholds']['lower_elev'],
        config['melt_thresholds']['lower_melt'],
        config['melt_thresholds']['upper_elev'],
        config['melt_thresholds']['upper_melt']
    )

    Th = create_hydrofracture_mult_raster(
        dem,
        config['melt_thresholds']['lower_elev'],
        config['melt_thresholds']['upper_elev'],
    )

    # ------------------------------------------------------------------
    # Diagnostic plots
    if config['outputs']['diagnostic_plots']:

        plt.figure()
        plt.imshow(Tm)
        plt.colorbar()
        plt.title('Melt (m w.e.) threshold for "hydrofracture"')
        plt.savefig(output_to / 'threshold_mwe.pdf')

        plt.figure()
        plt.imshow(Th)
        plt.colorbar()
        plt.title('Frac. of cum(melt) which can drain by hydrof.')
        plt.savefig(output_to / 'threshold_frac_cummelt.pdf')

    # ------------------------------------------------------------------
    # Prepare SEB melt data

    # Subset requested date range
    melt = melt.loc[
        config['date_range']['start']:
        config['date_range']['finish']
    ]

    # Hourly → daily sums
    melt_resampled_sum = melt.resample('1D').sum()

    # Remove leap day
    melt_resampled_sum = melt_resampled_sum[
        ~(
            (melt_resampled_sum.index.month == 2) &
            (melt_resampled_sum.index.day == 29)
        )
    ]

    # ------------------------------------------------------------------
    # Pad to full daily index

    full_index = pd.date_range(
        start=config['date_range']['start'],
        end=config['date_range']['finish'],
        freq='D'
    )

    # Remove leap day from target index
    full_index = full_index[
        ~(
            (full_index.month == 2) &
            (full_index.day == 29)
        )
    ]

    # Fill missing days with zero melt
    melt_resampled_sum = melt_resampled_sum.reindex(
        full_index,
        fill_value=0
    )

    # ------------------------------------------------------------------
    # Set up interpolation matrix

    X = melt_resampled_sum.columns
    X = sm.add_constant(X)

    # ------------------------------------------------------------------
    # Initialize storage arrays

    accum_sum = np.zeros(dem.shape)

    # Whether pixel has drained yet
    drained = np.zeros(dem.shape)

    return (
        Tm,
        Th,
        melt_resampled_sum,
        X,
        accum_sum,
        drained
    )


def daily_sheet_input(
    data,
    X,
    dem,
    accum_sum,
    drained,
    Tm,
    Th,
    config
):
    """
    Compute daily sheet input for subglacial hydrology forcing.

    Parameters
    ----------
    data : array-like
        Daily melt observations used for regression.
    X : pandas.DataFrame or ndarray
        Predictor matrix with constant term.
    dem : ndarray
        DEM grid.
    accum_sum : ndarray
        Accumulated melt prior to hydrofracture.
    drained : ndarray (bool)
        Boolean grid indicating whether each pixel has drained before.
    Tm : ndarray
        Melt threshold raster for hydrofracture triggering.
    Th : ndarray
        Hydrofracture multiplier raster.
    config : dict
        Configuration dictionary.

    Returns
    -------
    today_supply_total : ndarray
        Total daily water supply (m/s).
    accum_sum : ndarray
        Updated accumulated melt grid.
    melt_sum_grid : ndarray
        Daily melt sum raster.
    today_hydrofracture_supply_ms : ndarray
        Hydrofracture supply rate (m/s).
    today_melt_ms : ndarray
        Distributed melt supply rate (m/s).
    """

    # ------------------------------------------------------------------
    # Spread today's melt over the domain using linear regression

    f = sm.OLS(data, X).fit()

    melt_sum_grid = (
        float(f.params['x1']) * dem +
        float(f.params['const'])
    )

    # Prevent negative melt
    melt_sum_grid = np.maximum(melt_sum_grid, 0)

    # ------------------------------------------------------------------
    # Elevation-dependent melt scaling

    scale_lower = 1300
    scale_upper = 1800

    scale_factor = config['melt_thresholds']['melt_scaling']

    slope = (
        (scale_factor - 1) /
        (scale_upper - scale_lower)
    )

    scale = np.ones_like(dem, dtype=float)

    # Linear scaling between bounds
    scale = np.where(
        (dem > scale_lower) & (dem < scale_upper),
        1 + slope * (dem - scale_lower),
        scale
    )

    # Constant scaling above upper elevation
    scale = np.where(
        dem >= scale_upper,
        scale_factor,
        scale
    )

    # Apply scaling
    melt_sum_grid *= scale

    # ------------------------------------------------------------------
    # Update accumulated melt

    accum_sum += np.where(
        drained == False,
        melt_sum_grid,
        0
    )

    # ------------------------------------------------------------------
    # Hydrofracture drainage

    today_hydrofracture_supply = np.where(
        accum_sum >= Tm,
        accum_sum,
        0
    )

    # Remove drained water from accumulation
    accum_sum -= today_hydrofracture_supply

    # Convert daily amount to m/s
    today_hydrofracture_supply_ms = (
        (today_hydrofracture_supply * Th) /
        (60 * 60 * 24)
    )

    # ------------------------------------------------------------------
    # Standard melt delivery

    today_melt_ms = melt_sum_grid / (60 * 60 * 24)

    # Only supply where drainage already active
    today_melt_ms = np.where(
        drained == True,
        today_melt_ms,
        0
    )

    # ------------------------------------------------------------------
    # Total daily supply

    today_supply_total = (
        today_hydrofracture_supply_ms +
        today_melt_ms +
        3.17e-11
    )

    return (
        today_supply_total,
        accum_sum,
        melt_sum_grid,
        today_hydrofracture_supply_ms,
        today_melt_ms
    )
