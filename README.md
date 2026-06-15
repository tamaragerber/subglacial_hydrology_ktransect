# subglacial_hydrology_ktransect

## generate_forcing_files.py
Generate forcing files for Elmer/Ice GlaDS simulations over the k-transect, SW Greenland. 
The script reads in a config file (e.g. config_2009.yaml) which contains all relevant parameters and data, including
  - paths to input files, namely digital elevation model (DEM), surface runoff data, surface velocity files, GNSS data files for each station
  - basal sliding fraction
  - melt thresholds
  - lake drainage parameters

The generated forcing files contain basal sliding velocity and sheet input (surface water routed to the distributed subglacial system). 
Basal sliding velocity is assumed to be a defined fraction from surface flow velocities and are locally tuned using GNSS observations. 
The sheet input consists of a constant background value accounting for geothermal heat, surface water that hydrofractures and (optional) lake drainage events. Melt water accumulates at the surface until a hydrofracture threshold is reached, after which it is directly routed into the sheet. Below the lower elevation threshold all water is reaching the bed, above the higher elevation threshold, none of the water reaches the bed. In between the two, the fraction of water reaching the bed is linearly scaled. 

