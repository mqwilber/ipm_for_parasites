## Description
## -----------
## Script runs two primary analyses. First, it propagates the uncertainty
## in the parameter estimates to calculate the uncertainty in the population
## growth rate.  Second, it calculates the sensitivity of the population level
## growth rate and variance to mean ratio to the various
## lower level parameters used in the host-parasite IPM and calculates the
## uncertainty in these sensitivities.

## Notes: Script assumes that the working direction is ipm_models/code

## Author: Mark Wilber


source("IPM_functions_for_R.R")
library(ggplot2)
library(plyr)

params = readRDS("../results/IPM_parameters_no_26_linear.rds")
linear=TRUE

### Growth rate uncertainty ###

all_temps = 4:20
desired_temp = 20

# Lower and upper bounds
min_size = -5
max_size = 18
delta = 0.001
SAMPS = 1000

# Number of cells in the discretized matrix Smaller matrix = Faster calculation,
# but worse approximation. I have tested it with both 100 and 50 and the results
# are the same. Using 50 for faster simulation.
bins = 50

lambdas_ests = array(NA, dim=c(SAMPS, length(all_temps)))
uq = function(x){quantile(x, 0.75)}
lq = function(x){quantile(x, 0.25)}

# Lambda uncertainty
for(i in 1:length(all_temps)){

    desired_temp = all_temps[i]
    for(j in 1:SAMPS){

        # Sample the params
        samp_params = monte_carlo_params(params)

        # Specific parameters for each temperature
        params_base = set_temp_params(desired_temp, samp_params, linear)

        # Build full_P
        full_P_base = build_full_P(min_size, max_size, bins, params_base)

        lam_base = Re(eigen(full_P_base)$values[1])
        lambdas_ests[j, i] = lam_base

    }

}

### Plot the lambdas and the first 1, 2, and 3 quartiles ###
median_lambda = apply(lambdas_ests, 2, median)
lambda_upper = apply(lambdas_ests, 2, uq)
lambda_lower = apply(lambdas_ests, 2, lq)
lambda_data = data.frame(temp=all_temps, median=median_lambda, upper=lambda_upper, lower=lambda_lower)
lplot = ggplot(data=lambda_data, aes(x=all_temps, y=median_lambda)) + geom_point()
lplot = lplot + geom_errorbar(aes(ymax=upper, ymin=lower), alpha=0.2) + theme_bw()
lplot + ylim(0.80, 1) + xlab("Temperature, C") + ylab(expression("Population growth rate, " *lambda))

ggsave("../results/population_growth_rate_error.pdf", width=6, height=5)

# Plot expected time to pop below 0.1 %
median_time = (log(0.1) / log(median_lambda)) * 3
upper_time = (log(0.1) / log(lambda_upper)) * 3
lower_time = (log(0.1) / log(lambda_lower)) * 3

line_data = data.frame(temp=3:21, line=120)

time_data = data.frame(temp=all_temps, median=median_time, upper=upper_time,
                            lower=lower_time)
tplot = ggplot(data=time_data, aes(x=temp, y=median)) + geom_point() + scale_y_log10()
tplot = tplot + geom_errorbar(aes(ymax=upper, ymin=lower), alpha=0.2) + theme_bw()
tplot = tplot + xlab("Temperature, C") + ylab("Expected # of days to 90% decline")
tplot = tplot + geom_line(data=line_data, aes(x=temp, y=line), color="red", linetype="dashed")
tplot

ggsave("../results/extinction_times.pdf", width=6, height=5)

##############################################################################

### Lower level parameter sensitivity and uncertainty ###

param_list = c('surv_int', 'surv_slope', 'growth_int', 'growth_temp',
                    'growth_size', 'loss_int', 'loss_size', 'loss_temp',
                    'prob_inf_int', 'prob_inf_temp',
                    'clump_int', 'clump_temp', "class_zero_surv", "growth_sigma2",
                    'growth_sigma2_exp', 'clump_sigma2', 'clump_sigma2_exp')

# Dataframe to hold results for growth rate
full_elas = data.frame(param=character(), temp=numeric(), lower=numeric(),
                                            median=numeric(), upper=numeric())

# Dataframe to hold results for variance to mean ratio
full_elas_vm = data.frame(param=character(), temp=numeric(), lower=numeric(),
                                            median=numeric(), upper=numeric())

SAMPS = 1000
delta = 0.001
all_temps = c(10, 16, 20)
y = set_discretized_values(min_size, max_size, bins)$y

for(k in 1:length(all_temps)){

    desired_temp = all_temps[k]

    # Arrays for lambda elas
    elas_est = array(NA, dim=c(SAMPS, length(param_list)))
    sens_est = array(NA, dim=c(SAMPS, length(param_list)))

    # Arrays for V:M elas
    elas_est_vm = array(NA, dim=c(SAMPS, length(param_list)))
    sens_est_vm = array(NA, dim=c(SAMPS, length(param_list)))
    print(paste("Working on temp", desired_temp))

    for(i in 1:length(param_list)){

        param = param_list[i]

        for(j in 1:SAMPS){

            # Sample the params
            samp_params = monte_carlo_params(params)

            # Specific parameters for each temperature
            params_base = set_temp_params(desired_temp, samp_params, linear)

            # Build full_P
            full_P_base = build_full_P(min_size, max_size, bins, params_base)

            lam_base = Re(eigen(full_P_base)$values[1])
            vm_base = get_vm_ratio(full_P_base, y)

            # Alter the given parameter and calc sensitivity and elastiicty
            # Using the elasticity formula given in Merow et al. 2014
            params_alter = samp_params
            params_alter[[param]] = params_alter[[param]] + (delta * params_alter[[param]])
            params_elas = set_temp_params(desired_temp, params_alter, linear)

            full_P_elas = build_full_P(min_size, max_size, bins, params_elas)
            lam_elas = Re(eigen(full_P_elas)$values[1])
            vm_elas = get_vm_ratio(full_P_elas, y)

            sens_lam = (lam_base - lam_elas) / delta
            elas_lam = (1 / lam_base) * sens_lam

            sens_vm = (vm_base - vm_elas) / delta
            elas_vm = (1 / vm_base) * sens_vm

            # Save lambda sensitivities
            sens_est[j, i] = abs(sens_lam)
            elas_est[j, i] = abs(elas_lam)

            # Save variance:mean sensitivities
            sens_est_vm[j, i] = abs(sens_vm)
            elas_est_vm[j, i] = abs(elas_vm)

        }

    }

    # Save lambda results
    q_lam = compute_quartiles(elas_est)
    sum_dat_lam = data.frame(param=param_list, temp=desired_temp, lower=q_lam$lower,
                                        median=q_lam$median, upper=q_lam$upper)
    full_elas = rbind(full_elas, sum_dat_lam)

    # Save vm results
    q_vm = compute_quartiles(elas_est_vm)
    sum_dat_vm = data.frame(param=param_list, temp=desired_temp, lower=q_vm$lower,
                                        median=q_vm$median, upper=q_vm$upper)
    full_elas_vm = rbind(full_elas_vm, sum_dat_vm)

}

param_list_clean = c(expression("s(x):" ~ b["0,0"]),
                     "s(x): load",
                     expression("G(x', x):" ~ b["0,1"]),
                     "G(x', x, T): temperature",
                    "G(x', x, T): load",
                    expression("l(x):" ~ b["0,2"]),
                    'l(x, T): load',
                    'l(x, T): temperature',
                    expression(phi(T) ~ ":" ~ b["0,4"]),
                    expression(phi(T) ~ ":" ~ "temperature"),
                    expression(paste(G[0], "(x', T): ", b["0,3"])),
                    expression(paste(G[0], "(x', T): ", "temperature")),
                    expression(s[0]),
                    expression("G(x', x, T): variance, " ~ nu["0,1"]),
                    "G(x', x, T): variance, load effect",
                    expression(paste(G[0], "(x', T): variance, ") ~ nu["0,3"]),
                    expression(paste(G[0], "(x', T): variance, temperature effect"))
                    )

                    #expression(paste(G[0], "(x', T): Variance, a")),

# Plot lambda elasticity
full_elas$temp = as.factor(full_elas$temp)
full_elas$temp = revalue(full_elas$temp, c("10"="10 C", "16"="16 C", "20"="20 C"))

sens = ggplot(data=full_elas, aes(x=param, y=median)) + geom_point() +
                    geom_errorbar(aes(ymax=upper, ymin=lower), alpha=0.2) +
                    facet_wrap(~temp) + theme_bw() +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=1, size=5)) +
                    xlab("Vital rate (IPM) parameters") +
                    ylab(expression("Elasticity of " *lambda))

sens = sens + scale_x_discrete(labels=param_list_clean[order(param_list)])
sens

ggsave("../results/lower_level_sensitivity.pdf", width=8, height=6)


# Plot Variance to Mean elasticity
full_elas_vm$temp = as.factor(full_elas_vm$temp)
full_elas_vm$temp = revalue(full_elas_vm$temp, c("10"="10 C", "16"="16 C", "20"="20 C"))

vm_plot = ggplot(data=full_elas_vm, aes(x=param, y=median)) + geom_point() +
                    geom_errorbar(aes(ymax=upper, ymin=lower), alpha=0.2) +
                    facet_wrap(~temp) + theme_bw() +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) +
                    xlab("Vital rate (IPM) parameters") +
                    ylab("Elasticity of variance:mean ratio")

vm_plot = vm_plot + scale_x_discrete(labels=param_list_clean[order(param_list)])
vm_plot

ggsave("../results/lower_level_sensitivity_vm.pdf", width=8, height=6)





