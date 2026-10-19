# subglacial_hydrology_ktransect

## Forcing/generate_forcing_filey.py

**Dependencies:** `load_data.py`, `velocities.py`, `lake_drainage.py`, `surface_hydrology.py` 

Generate forcing files for Elmer/Ice GlaDS simulations over the k-transect, SW Greenland. 
The script reads in a config file (e.g. config/config_2009.yaml) which contains all relevant parameters and data, including:
  - paths to input files, namely digital elevation model (DEM), surface runoff data, surface velocity files, GNSS data files for each station
  - basal sliding fraction
  - melt thresholds
  - lake drainage parameters

The generated forcing files contain basal sliding velocity and sheet input (surface water routed to the distributed subglacial system). 
Basal sliding velocity is assumed to be a defined fraction from surface flow velocities and are locally tuned using GNSS observations. The local transient basal sliding can be turned on and off. 
The sheet input consists of a constant background value accounting for geothermal heat, surface water that hydrofractures and (optional) lake drainage events. Melt water accumulates at the surface until a hydrofracture threshold is reached, after which it is directly routed into the sheet. Below the lower elevation threshold all water is reaching the bed, above the higher elevation threshold, none of the water reaches the bed. In between the two, the fraction of water reaching the bed is linearly scaled.

To reproduce forcing files you'll need to adjust the paths for input files in the config.yaml file to following datasets:
-  digital elevation model; Porter, Claire, et al., 2023, “ArcticDEM, Version 4.1”, https://doi.org/10.7910/DVN/3VDC4W, Harvard Dataverse, V1, [12.06.2026].
-  modeled surface runoff; van As, D., Bech Mikkelsen, A., Holtegaard Nielsen, M., Box, J. E., Claesson Liljedahl, L., Lindbäck, K., Pitcher, L., and Hasholt, B.: Hypsometric amplification and routing moderation of Greenland ice sheet meltwater release, The Cryosphere, 11, 1371–1386, https://doi.org/10.5194/tc-11-1371-2017, 2017. 
- surface velocities (converted to m/s); Joughin, I., Smith, B., Howat, I., & Scambos, T. (2016). Measures multi-year greenland ice sheet velocity mosaic, version 1. NASA National Snow and Ice DataCenter Distributed Active Archive Center. doi: 10.5067/QUA5Q9SVMSJG
- gnss surface velocities; Tedstone, A., Nienow, P., Sole, A., Bartholomew, I., Cowton, T., Mair, D., & King, M. (2026). Greenland Ice Sheet displacement on the K-Transect, south-west Greenland Ice Sheet from in-situ GNSS measurements (Version 3.0) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/0adac071-a2fa-41a4-952f-4a418fa32e2e and https://zenodo.org/records/18611668?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MTc3MzgzNTcwMSwiZXhwIjoxNzk4Njc1MTk5fQ.eyJpZCI6IjUxOTA4ZWE2LTY0MGEtNGY4OC1iY2I2LTc2ODAwNzEyNzQ4OCIsImRhdGEiOnt9LCJyYW5kb20iOiJmYTZmYjg5NDgxY2FlMTk4N2JlOTJlNTljMmYwNDYyMSJ9.zrTNreM4ftSFPw3aEJjPbkjWQKLcKl5Iwz_Q3xxiQKyP6YXqsaNSs8We3ThrxrrihHMO2oPM8D5uhNwwPrbqjg


Adjust line 128 in generate_forcing_files.py to generate forcing files for different years. 
  


## Example_sif_files
This folder contains sif files to reproduce the simulations with Elmer/Ice. The domain is split into 8 partitions, the "domain" folder contains the mesh files read by Elmer/Ice. 

- 00_spinup.sif: initializes the domain. This is not a real simulation but maps the basal and surface geometry as well as surace velocities to the mesh. The paths for topography (Bedmachine v3, 10.5067/2CIX82HUV88Y) and surface velocity (Measures multi-year mosaic, converted to m/s; doi: 10.5067/QUA5Q9SVMSJG) data need to be adjusted in read_parameters.in
- 01_spinup_08_6mm.sif: is the file for spinup using constant sheet input of 6 mm/yr.
- 08_run_2009_ii.sif: runs a simulation for the year 2009, restarting from the final state of the spinup. The paths for the forcing files in Solver 4 need to be adjusted.
- parameters_mks_hkc_llc.lua contains the model parameters.
  
