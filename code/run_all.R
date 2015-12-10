## Description
## -----------
## Run all analyses and make all figures for host-parasite IPM analysis
##
## Author: Mark Wilber
#
## Notes : The entire analysis takes ~30 minutes on a Macbook
## with 2.6 GHz Intel Core i5

# Set the working directory to ipm_models/code
# NOTE: THIS WILL NEED TO BE CHANGED FOR ANY PARTICULAR MACHINE

setwd("/Users/mqwilber/Repos/ipm_models/code")

# Check for necessary packages and install them if necessary
list.of.packages <- c("ggplot2", "MASS", "reshape2", "RColorBrewer", "grid",
                        "fields", "popbio", "Matrix", "AICcmodavg", "boot",
                        "nlme", "lme4", "plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Calculate parameters
print("Calculating vital rate parameters...")
source("ipm_vital_rate_parameter_analysis.R")
print("...done calculating parameters")

# Analyze IPM
print("Analyzing the host-parasite IPM...")
source("ipm_bd_analysis_for_manuscript.R")
print("...done with the IPM analysis")

# # Lower level sensitivity analysis
print("Propagating error and running a sensitivity analysis...")
source("lower_level_sensitivity_analysis.R")
print("...done with error propagation and sensitivity")

# # Density-dependent simulation
print("Running density-dependent transmission simulation...")
source("transmission_simulation.R")
print("...completed density-dependent simulation")

# # Make density dependent plots
print("Making density-dependent figures...")
source("transmission_figures.R")
print("...done with transmission figures")

# # Make trajectory plots
print("Making trajectory plots...")
source("trajectories_plot.R")
print("...finished trajectory plots")

print("Completed IPM analyses. Results are saved in ../results")
