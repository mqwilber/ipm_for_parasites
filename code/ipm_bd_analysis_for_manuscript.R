## Description
## ------------
## This script used the parameters calculated in the script
## `ipm_vital_rate_parameter_analysis.R` and builds the host-parasite IPM model.
## It then runs some basic IPM analysis such as mean time to extinction,
## population growth rate, stable distributions, etc.

## Notes: All results from this script are stored in ../results.  This script
## assumes that the working directory is ipm_models/code

## Author: Mark Wilber

#############################################################################

# Source the relevant functions (CHECK YOUR WORKING DIRECTORY!)
source("eviction_functions.R")
source("IPM_functions_for_R.R")

# Read in the fitted parameters that were obtained by running
# ipm_vital_rate_parameter_analysis.R

# Parameters for linear growth function
params_no_26_l = readRDS("../results/IPM_parameters_no_26_linear.rds")

# Parameters for nonlinear growth function
params_no_26_nl = readRDS("../results/IPM_parameters_no_26_nonlinear.rds")

# Decision parameter. If FALSE, considers a non-linear growth curve
linear = TRUE

# Specify parameters based on temperature relationship
if(linear){
  tparams = params_no_26_l
} else{
  tparams = params_no_26_nl
}

all_temps = 4:20 # Vector of temperatures
#all_temps = c(4, 12, 20)
days = 3 # length of time between swabs

# Lower and upper bounds (chosen to minimize eviction)
min_size = -5
max_size = 18

# Number of cells in the discretized matrix. Adding more just gives more resolution
bins = 100

# Arrays to save results
stable_dist_results = array(NA, dim=c(length(all_temps), bins + 1))
sd_means = array(NA, dim=length(all_temps))
sd_vars = array(NA, dim=length(all_temps))
absorption_results = array(NA, dim=c(length(all_temps), bins + 1))
absorption_var_results = array(NA, dim=c(length(all_temps), bins + 1))
lambdas = array(NA, dim=length(all_temps))
evict_dlambdas = array(NA, dim=length(all_temps))
evict_max = array(NA, dim=length(all_temps))
elasticity_results = array(NA, dim=c(length(all_temps), bins, bins))

# For loops to perform calculations for multiple temperatures
for(i in 1:length(all_temps)){

  desired_temp = all_temps[i]

  # Specific parameters for each temperature
  params = set_temp_params(desired_temp, tparams, linear)

  ##############################################################################

  # Exploring the IPM model

  # In this section we analyze the IPM model by discretizing it and treating it
  # like a matrix population model. Using functions defined in
  # "IPM_functions_for_R"

  # Get the relevant matrix parameters: midpoints, cell edges, and cell width
  matrix_params = set_discretized_values(min_size, max_size, bins)
  bnd = matrix_params$bnd
  y = matrix_params$y
  h = matrix_params$h

  # Calculate and plot eviction.  Eviction is most important when e is large
  # and. Using the script provided by Williams et al. for checked for eviction

  k_xpx = function(xp, x, params){

    # Define the full kernel (for eviction calculations)

    return(g_xpx(xp, x, params) * s_x(x, params) * (1 - r_x(x, params)))

  }

  evict_results = evictionMeasuresFC.Iter(g_xpx, k_xpx, s_x, min_size, max_size, params)
  evict_dlambdas[i] = evict_results$dlambda
  evict_max[i] = max(evict_results$evict)
  #plot(evict_results$y, evict_results$evict)

  # Get the kernel of the IPM without the zero class
  kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=T,
        save=paste("../results/kernel_plots_temp", all_temps[i], sep=""))
  P = kernels$P

  ##############################################################################

  # Add on the uninfected class

  # Compute the uninfected to infected column.  The zero to zero transition is
  # given by surviving with a load of zero ($s(0)$) and not becoming infected (1
  # - $p_{inf}$). The other transitions are given by surviving the zero class
  # ($s(0)$) and becoming infected ($p_{inf}$) and gaining a given load $x'$
  # ($\phi(x')$), which needs to be calculated with the midpoint rule.

  # The zero -> zero transition
  zero_trans = params$class_zero_surv * (1 - params$prob_inf)

  # The zero to nonzero transitions: Transitioning to zero to infected
  prob_clump = h * dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * params$prob_inf * prob_clump
  col1 = c(zero_trans, ztnz)


  # Compute the infected to uninfected row. This is just the probability of
  # surviving with load $x$ ($s(x)$) and then lossing an infection $r(x)$.  We
  # just need to multiply $s(x) r(x)$

  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)

  # Now we just need to stick these extra columns onto the **P** matrix. The
  # most likely thing that is going to happen is you are going to stay in the
  # uninfected class.

  full_P = get_full_P(P, row1, col1, min_size, y, plot_it=F)

  ##############################################################################

  # Now we need to get some of the predictions from the model

  # Get the mean absorption times and the variance
  absorp_times = absorption_times(full_P)
  absorption_results[i, ] = absorp_times$mean
  absorption_var_results[i, ] = absorp_times$var

  # Get the sensitivities and elasticities of lambda on each tranistion
  # probability.

  # The dominant eigenvalue
  lam = Re(eigen(full_P)$values[1])
  lambdas[i] = lam

  # Get the elasticity and sensitivity
  es_res = get_elasticity_and_sensitivity(P, h, y,
                plot_it=T,
                save=paste("../results/elas_temp",
                all_temps[i], sep=""))

  elasticity_results[i, , ] = es_res$elas

  # Stable distribution of individuals
  stable_dist = get_stable_dist(full_P)
  stable_dist_results[i, ] = stable_dist

  # Get mean and variance conditional on infection
  cond_stable_dist = stable_dist[2:length(stable_dist)] / (1 - stable_dist[1])
  cond_mean = sum(cond_stable_dist * y)
  cond_var = sum(cond_stable_dist * y^2) - (cond_mean)^2
  sd_means[i] = cond_mean
  sd_vars[i] = cond_var

} # End for loop for all temperatures

#######

# Now just making a bunch of plots using the results from the above analysis

#######

# Get the stable distributions conditional on infection
cond_dists = array(NA, dim=c(length(all_temps), bins))
for(i in 1:dim(cond_dists)[1]){
  cond_dists[i, ] = stable_dist_results[i, 2:(bins + 1)] /
          (1 - stable_dist_results[i, 1])
}

# Plot and save the stable distributions conditional on infection
df_stable = data.frame(t(cond_dists))
colnames(df_stable) = all_temps
stacked_stable = stack(df_stable)
colnames(stacked_stable) = c("prob", "Temperature")
x_vals = rep(y, length(all_temps))

library(ggplot2)
stacked_stable$Temperature = factor(stacked_stable$Temperature, levels=c(4, 8, 12, 16, 20))
ind = !is.na(stacked_stable$Temperature)
gdist_plot = ggplot(data=stacked_stable[ind, ], aes(x=x_vals[ind], y=prob)) +
            geom_line(aes(linetype=Temperature, color=Temperature)) + xlab("Log zoospore load") +
            ylab("Probabiity") + theme_bw()
gdist_plot$labels$colour = "Temperature, C"
gdist_plot$labels$linetype = "Temperature, C"
gdist_plot
ggsave("../results/conditional_stable_distributions.pdf", width=6, height=5)

# Plot and save the mean variance relationship
# pdf("../results/mean_var_relationship.pdf")
# plot(all_temps, sd_vars / sd_means, pch=19,
#                     xlab="Temperature", ylab="Variance / Mean")
# lines(all_temps, sd_vars / sd_means)
# dev.off()

# Plot and save the absoprtion times
df_absorp = data.frame(t(absorption_results))
trunc = floor(0.2 * bins):dim(df_absorp)[1]
df_absorp = df_absorp[trunc, ]
colnames(df_absorp) = all_temps
stacked_absorp = stack(df_absorp)
colnames(stacked_absorp) = c("mean_time", "Temperature")

# # Add in std
# df_absorp_var = data.frame(t(absorption_var_results))
# df_absorp_var = df_absorp_var[trunc, ]
# colnames(df_absorp_var) = all_temps
# stacked_absorp_var = stack(df_absorp_var)
# colnames(stacked_absorp_var) = c("var_time", "Temperature")

# stacked_absorp['std_time'] = sqrt(stacked_absorp_var$var_time)

x_vals = rep(c(min_size - 1, y)[trunc], length(all_temps))

library(ggplot2)
stacked_absorp$Temperature = factor(stacked_absorp$Temperature, levels=c(4, 8, 12, 16, 20))
ind = !is.na(stacked_absorp$Temperature)
stacked_absorp_trunc = stacked_absorp[ind, ]
x_vals = x_vals[ind]

ggplot(data=stacked_absorp_trunc, aes(x=x_vals, y=mean_time * days)) +
            geom_line(aes(linetype=Temperature)) + xlab("Log zoospore load") +
            scale_y_log10() +
            ylab("Mean time to death (days)") + theme_bw()
ggsave("../results/mean_time_to_death.pdf", width=6, height=5)


# Plot and save lambdas
lambda_data = data.frame(cbind(all_temps, lambdas))
ggplot(lambda_data, aes(x=all_temps, y=lambdas)) + geom_line() +
                  xlab("Temperature (C)") + ylab(expression("Population growth rate, " *lambda)) + theme_bw()
ggsave("../results/lambda_plot.pdf")

# Plot and save conditional means and variances
agg_data = data.frame(cbind(all_temps, sd_vars, sd_means))
ggplot(agg_data, aes(x=all_temps, y=sd_means)) + geom_line() +
              xlab("Temperature (C)") + ylab("Mean log zoospore load") +
              theme_bw()
ggsave("../results/mean_plot.pdf")

# Plot and save conditional means and variances
agg_data = data.frame(cbind(all_temps, sd_vars, sd_means))
ggplot(agg_data, aes(x=all_temps, y=sd_vars)) + geom_line() +
              xlab("Temperature (C)") + ylab("Variance in log zoospore load") +
              theme_bw()
ggsave("../results/var_plot.pdf")

# Plot and save conditional means and variances
agg_data = data.frame(cbind(all_temps, sd_vars, sd_means))
ggplot(agg_data, aes(x=all_temps, y=sd_vars / sd_means)) + geom_line() +
              xlab("Temperature (C)") + ylab("Variance : mean ratio") +
              theme_bw()
ggsave("../results/var_mean_plot.pdf", width=6, height=5)





