library(fields)

# Source my IPM functions
source("/Users/mqwilber/Repos/ipm_models/code/IPM_functions_for_R.R")
trans_folder = "/Users/mqwilber/Repos/ipm_models/results/transmission_results/"

# Load in the parameters for the IPM
params_full = read.csv("/Users/mqwilber/Repos/ipm_models/results/IPM_parameters_no_26_linear.rds")

params_full$surv_int = 10.8241
# Get the basic matrix parameters for the simulation

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
SUMMER_TEMP = 12
linear_temp = TRUE

# Transmission function parameters
num_vals = 40
beta_vals = c(0) #seq(0, 1, len=num_vals)
inf_vals =  c(0.3) #seq(0.01, .99, len=num_vals)

#### BEGIN SIMULATION ####

params = set_temp_params(SUMMER_TEMP, params_full, linear_temp)

# Get the time invariant portion of the kernel
kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=F)
P = kernels$P

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

            results[, i + 1] =  one_step_out$nt
            all_inf_probs[i] = one_step_out$inf_prob


        }

    }

}



# Plots for an individual run over the course of the summer

# Going to plot population size and probability of infection


par(mfrow=c(2, 2))
time = seq(0:40) * 3
plot(time, colSums(results), type="l")

# ZE_total = matrix(c(0, exp(y)), nrow=1, ncol=bins + 1) %*% results
# plot(c(ZE_total), type="l")
sz = dim(results)[2]

# Plot distribution of zoospores on infected hosts
trun = sz - 1 #YEAR_LENGTH
Ns = colSums(results[2:(bins + 1), (sz - trun):sz])
matplot(y, t(t(results[2:(bins + 1), (sz - trun):sz]) / Ns), type='l',
                xlab="ln(Zoospore Load)", ylab="Proportion")

norm_dist = t(t(results[2:(bins + 1), (sz - trun):sz]) / Ns)
sz2 = dim(norm_dist)[2]
norm_means = array(NA, dim=sz2)
norm_vars = array(NA, dim=sz2)

for(i in 1:sz2){
    norm_means[i] = sum(norm_dist[, i] * y)
    norm_vars[i] = sum(norm_dist[, i] * y^2) - norm_means[i]^2
}

plot(time, norm_vars / norm_means)
plot(time, results[1, ] / colSums(results), type="l", ylab="Prop. uninfected",
                ylim=c(0, 1))


