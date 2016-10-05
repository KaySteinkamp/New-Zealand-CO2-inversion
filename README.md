# Bayesian inversion of CO_2 observations in New Zealand

The provided MATLAB code performs a Bayesian inversion of daily CO_2 records from sites in New Zealand.

Computed are geographically distributed, regional sources and sinks for CO_2 with predetermined within-region
spatial pattern and weekly (optionally monthly) temporal resolution. The observations are assimilated for the years 2011-2013, or optionally a subset thereof.

The input data is the observational response for each source (from Lagrangian disperion backward modelling), a series of prior estimates 
for each source (along with variances) and a record of observations (daily afternoon CO_2 measurements at up to 3 sites through 2011-13).

The inversion uses Green's functions, which are are precalculated by Lagrangian dispersion model running backward mode, and using hourly
meteorology from a 10km x 10km weather forecast model.

The computational core of the code is the routine bayesinv.m which is a general-purpose Bayesian inverse solver.

# Usage
There is nothing to install, all necessary data are available as mat files already. Set you MATLAB path appropriately after cloning the repo. Then adjust the local environment variables in the config and doinv.m files, and you are good to go.

Run doinv.m in MATLAB. Adjust settings in config.

Output folder, along with results and figures, are automatically created in the local folder.

