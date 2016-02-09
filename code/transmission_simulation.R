## Description
## -----------
## This script simulates the density-dependent IPM with an environmental
## reservoir.  The simulation loads in the vital rate parameters specified
## for the IPM and then runs the density-dependent simulation for temperatures
## specified summer_temps. The results of the simulations are saved as csvs
## for easy plotting with the script `transmission_figures .R`.

## Author: Mark Wilber

#############################################################################

# Source the IPM functions
source("IPM_functions_for_R.R")
trans_folder = "../results/transmission_results/"

# Load in the parameters for the IPM
params_full = readRDS("../results/IPM_parameters_no_26_linear.rds")

# Set the basic matrix parameters for the simulation

# Length of a summer is four months...120 days...at 3 days per time step...40 time steps
TIME_STEPS = 40 # length of summer

# Lower and upper bounds
min_size = -5
max_size = 18

# Number of cells in the discretized matrix
bins = 100
matrix_params = set_discretized_values(min_size, max_size, bins)
bnd = matrix_params$bnd
y = matrix_params$y
h = matrix_params$h

# Run for multiple summer temps
summer_temps = 12:20
linear_temp = TRUE

# Transmission function parameters
num_vals = 40
beta_vals = seq(0, 0.00117, len=num_vals)

# For visualization purposes, just calculate up to 0.3 for the environmental
# reservoir
#inf_vals = seq(0.01, .6, len=num_vals)
inf_vals = seq(0.01, .3, len=num_vals)

#### BEGIN SIMULATION ####

for(SUMMER_TEMP in summer_temps){

    print(paste("Beginning temperature ", SUMMER_TEMP))

    params = set_temp_params(SUMMER_TEMP, params_full, linear_temp)

    # Get the time invariant portion of the kernel
    kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=F)
    P = kernels$P

    # Empty arrays to hold results
    pop_size_post_summer = array(NA, dim=c(num_vals, num_vals))
    time_to_prev = array(NA, dim=c(num_vals, num_vals))
    prop_uninfected = array(NA, dim=c(num_vals, num_vals))

    # Loop through different transmission coefficients
    for(k in 1:length(beta_vals)){

        params$beta = beta_vals[k]

        # Loop through different environmental infection values

        for(j in 1:length(inf_vals)){

            # Set environmental infection probability
            params$prob_inf = inf_vals[j]

            # Array for results
            results = array(NA, dim=c(bins + 1, TIME_STEPS + 1))
            all_inf_probs = array(NA, dim=TIME_STEPS)

            # Start with an uninfected population
            init = 100
            n_0 = c(init, rep(0, bins))
            results[, 1] = n_0

            # Iterate population through the summer

            stored_prev = FALSE

            for(i in 1:TIME_STEPS){

                n_t = as.matrix(results[, i], nrow=bins + 1, ncol=1)

                # Perform the IPM simulation
                one_step_out = within_season_simulation(n_t, params, bins, bnd, y, h, kernel)

                results[, i + 1] = n_t_new = one_step_out$nt
                all_inf_probs[i] = one_step_out$inf_prob

                # Set time to prevalence
                prev = 1 - (n_t_new[1] / sum(n_t_new))

                if((stored_prev == FALSE) & (prev >= 0.95)){

                    time_to_prev[k, j] = i
                    stored_prev = TRUE

                }

            }

            if(stored_prev == FALSE){
                time_to_prev[k, j] = TIME_STEPS
            }

            pop_size_post_summer[k, j] = colSums(results)[dim(results)[2]]
            prop_unifs = results[1, ] / colSums(results)
            prop_uninfected[k, j] = prop_unifs[dim(results)[2]]

        }

    }

    # Save the results of the parameter exploration
    write.csv(data.frame(beta_vals), paste(trans_folder, "beta_vals.csv", sep=""))
    write.csv(data.frame(inf_vals), paste(trans_folder, "inf_vals.csv", sep=""))

    write.csv(data.frame(pop_size_post_summer),
        paste(trans_folder, "pop_size_post_summer_", SUMMER_TEMP, ".csv", sep=""))

    write.csv(data.frame(time_to_prev),
        paste(trans_folder, "time_to_prev_", SUMMER_TEMP, ".csv", sep=""))

    write.csv(data.frame(prop_uninfected),
        paste(trans_folder, "prop_uninfected_", SUMMER_TEMP, ".csv", sep=""))

}





