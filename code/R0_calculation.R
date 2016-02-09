## Description
## ------------
## Script for calculating R0 for the host-parasite IPM based.  The procedure
## is described completely in Supplementary Material 3.

## This script calculates R0 from the full IPM and the collapsed IPM

## Author: Mark Wilber
## Date: 01/2016

# Load the parameters
source("IPM_functions_for_R.R")
library(ggplot2)
library(plyr)

# Read and set params
params_full = readRDS("../results/IPM_parameters_no_26_linear.rds")

############################
# Parameters for IPM model #
############################

linear = TRUE
min_size = -5
max_size = 18
bins = 50

matrix_vals = set_discretized_values(min_size, max_size, bins)
y = matrix_vals$y # Centers of discretized matrix
h = matrix_vals$h # Delta between bin edges

S_init = 100 # Initial number of susceptible hosts
beta = 9.82 * 10^-4 # from Rach. and Briggs
s_0 = 1 # Survival probability is 1
temps = 4:20

# Arrays for storing results
R0s_approx = array(NA, dim=length(temps))
R0s_exact = array(NA, dim=length(temps))
sbars = array(NA, dim=length(temps))
lbars = array(NA, dim=length(temps))


################################################################
## Step 1: Get R0 based on collapsing the IPM to an SIS model ##
################################################################

for(i in 1:length(temps)){

    desired_temp = temps[i]

    params = set_temp_params(desired_temp, params_full, linear)
    params$prob_inf = 1 - exp(-9.82 * 10^-4) # Fix transmission

    P_matrix = build_full_P(min_size, max_size, bins, params)
    stable_dist = get_stable_dist(P_matrix)
    cond_dist = stable_dist[-1] / sum(stable_dist[-1])

    # Get mean loss function
    mean_loss = sum(sapply(1:bins, function(i) cond_dist[i] * r_x(y[i], params)))
    lbars[i] = mean_loss

    # Get mean survival
    mean_surv = sum(sapply(1:bins, function(i) cond_dist[i] * s_x(y[i], params)))
    sbars[i] = mean_surv

    prob_inf = params$prob_inf # Temperature dependent infection probability
    approx_beta = -log(1 - prob_inf)

    R0 = (beta * (s_0 * S_init)) / (1 - mean_surv * (1 - mean_loss))
    #R0 = (beta * (s_0 * S_init)) + mean_surv * (1 - mean_loss)
    R0s_approx[i] = R0

}

################################################################
## Step 2: Get the exact R0 for the host-parasite IPM based on #
## the matrix algebra given in Klepac and Caswell 2011        ##
################################################################

for(i in 1:length(temps)){

    desired_temp = temps[i]
    params = set_temp_params(desired_temp, params_full, linear)

    # Compute transition probabilities
    G = h * outer(y, y, g_xpx, params=params)
    S = s_x(y, params=params)
    R = r_x(y, params=params)
    U = G %*% diag(S * (1 - R)) # The transition matrix

    # Compute first term in linearization
    StoI = s_0 * h * g0_x(y, params) * beta * S_init
    SI = matrix(rep(StoI, bins), nrow=bins, ncol=bins)

    J = SI + U # Add the two terms to get J

    # Compute R0 using equation 29 in Klepac and Caswell 2011
    F = J - U # The same as SI
    R_mat = SI %*% solve(diag(bins) - U)
    eig_vals = Re(eigen(R_mat)$values)
    R0s_exact[i] = max(eig_vals)

}


####################################
## Step 3: Plot the two R0 results #
####################################

nt = length(temps)
labels = rep(c("Full R0", "Collapsed R0"), c(nt, nt))
R0_data = data.frame(r0=c(R0s_exact, R0s_approx), temp=c(temps, temps),
                        Method=labels)
ggplot(data=R0_data, aes(x=temp, y=r0)) + geom_line(aes(color=Method)) +
        theme_bw() + geom_hline(yintercept=1, linetype="dashed") +
        xlab("Temperature, C") + ylab(expression(R[0]))
ggsave("../results/r0_plot.pdf", width=7, height=5)

