# Metadata ===============================================================
# Description:Refactor of MaxEnt pipeline from PhD Chapter 1 for more generalised use
##
## R version: 4.3.2 for Windows
## Date: 2026-08-31 10:29:18.536453
## License: GPL3
## Author: Alan Vincent (alan.vincent@anu.edu.au)
#

## Features include:
# Creates background areas from a intersect between a set of polygons (like Biome polygon)
# SDM's require sensible background polygons, i.e. don't include coastal water areas for desert species. Overfits models otherwise.

# Creates and masks by a bias layer to control for data set level sampling bias.
# Considering including species by species level bias if it matters.

## Data Assumptions:
# No duplicated data points
# Data has been cleaned of spurious records
# Minimal co-linearity among env.vars variables if possible
# Environmental data is at appropriate resolution

## ENMEval controls for:
# Occurrence records sharing same grid cell. Removes extras to remove pseudo-replication
# Occurrence and background points with NA predictor (env) variables are removed.

## Parameters for final function
# Cleaned occurrence records of species of interest: "occ.data.sf"
# Terra Raster stack of environmental variables: "env.vars"
# An output directory.

## MUST HAVE JAVA INSTALLED ON THE MACHINE TO WORK (if you want a Maxent.jar model)
# See java_setup.R
