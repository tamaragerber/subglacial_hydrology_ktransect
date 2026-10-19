#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:45:14 2026

@author: tamara
"""
import numpy as np

def parse_lake_drainage(config):
    """
    reads lake drainage config
    """
    
    lake_cfg = config.get("lake_drainage", {})

    use_lakes = lake_cfg.get("drainage", False)

    lakes = []
    if use_lakes:
        for key in lake_cfg:
            if key.endswith("_volume"):
                prefix = key.replace("_volume", "")
                V = lake_cfg[key]
                doy = lake_cfg.get(f"{prefix}_doy", None)
                x0 = lake_cfg.get(f"{prefix}x", None)
                y0 = lake_cfg.get(f"{prefix}y", None)

                if V > 0 and doy is not None:
                    lakes.append({
                        "V": V,
                        "doy": doy,
                        "x": x0,
                        "y": y0
                    })

    return lakes
    

def generate_blister_on_grid(V, H_ice, x0, y0, xx, yy):
    """
    Returns sheet input field (m/s) for a blister
    """

    # Reference scaling
    R0 = 2000.0
    V0 = 1e7

    # Radius scaling
    R = R0 * (V / V0)**(1/6) * (H_ice / 1000.0)**(1/3)

    # Peak thickness
    H0 = 2 * V / (np.pi * R**2)

    # Distance from lake center
    r = np.sqrt((xx - x0)**2 + (yy - y0)**2)

    # Parabolic shape
    h = np.zeros_like(r)
    inside = r < R
    h[inside] = H0 * (1 - (r[inside]**2 / R**2))

    # Convert to m/s (1 day timestep)
    q = h / 86400.0

    return q
