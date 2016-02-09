## Description
## -----------
## Contains some helpful functions for analyzing an IPM model in R. See below
## for the various functions

## Author: Mark Wilber

##############################################################################

library(MASS)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(fields)
library(popbio)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

format_for_heat = function(mat, row_nm, col_nm, var_names){

    # Function converts a 2D matrix for plotting in ggplot

    tmat = mat
    colnames(tmat) = col_nm
    rownames(tmat) = row_nm
    long_mat = melt(tmat)
    colnames(long_mat) = var_names

    return(long_mat)

}

set_discretized_values = function(min_size, max_size, bins){

  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model

  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix

  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)


  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins

  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])

  # Width of cells
  h = y[2] - y[1]

  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
                h=h))

}



get_the_kernel = function(g_xpx, s_x, r_x, bins, y, params, h, plot_it=F,
              save="temp"){

  # Function to get the IPM kernel for the fungal model.  Doesn't include
  # the 0 class

  # Parameters
  # ----------

  # g_xpx : the growth function
  # s_x : the survival function
  # r_x : the loss function
  # bins : number of bins in the discretize model
  # y : midpoints (length(y) = bins)
  # params : list of relevant parameters to pass into vital functions
  # h : width of discretized cell. Used for midpoint evaluation

  # Returns
  # -------
  # : list
  #   P (kernel), G (growth Matrix), R (loss vector), S (survival vector)

  # Make growth matrix
  G = h * outer(y, y, g_xpx, params=params)

  # Make survival matrix
  S = s_x(y, params=params)

  # Make loss matrix
  R = r_x(y, params=params)

  # Make the full kernel
  #P = G
  # Faster without a for loop...
  #for(i in 1:bins) P[,i]=G[,i]* S[i] * (1 - R[i])

  P = G %*% diag(S * (1 - R))

  if(plot_it){

    #growth_plot =

    inf_data = data.frame(size=y, density=dnorm(y, mean=params$clump_mean, sd=params$clump_sd))
    infection_plot = ggplot(inf_data, aes(x=size, y=density)) + geom_line() + labs(title=bquote('Initial infection burden function, '*G[0](x)*'')) +
                      ylab("Probability density") +
                      xlab("Log zoospore load at time t") +
                      ylim(0, 0.45) + theme_bw() + theme(plot.title=element_text(size=12))

      #+ ggtitle(bquote('Initial infection function, '*G[0](x)*''))
    #'
    recov_data = data.frame(size=y, recovery=R)
    recovery_plot = ggplot(recov_data, aes(x=size, y=recovery)) + geom_line() + labs(title="Loss of infection function, l(x)") +
                      ylab("Probability of losing infection") +
                      xlab("Log zoospore load at time t") +
                      ylim(0, 1) + theme_bw() + theme(plot.title=element_text(size=12))
    lts = 6
    myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    #Format and plot G
    long_G = format_for_heat(t(G), y, y, c('size', 'sizeNext', 'values'))
    growth_plot = ggplot(long_G, aes(x=size, y=sizeNext, fill=values)) + geom_tile() + labs(title="Growth function, G(x', x)")
    growth_plot = growth_plot + scale_fill_gradientn(colours = myPalette(100), name="Probability")
    growth_plot = growth_plot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
    growth_plot = growth_plot + theme(legend.text=element_text(size=5),
                                      legend.title=element_text(size=lts),
                                      plot.title=element_text(size=12))

    # Format and plot P
    long_P = format_for_heat(t(P), y, y, c('size', 'sizeNext', 'values'))
    sg_plot = ggplot(long_P, aes(x=size, y=sizeNext, fill=values)) + geom_tile() + labs(title="Survival/growth kernel")
    sg_plot = sg_plot + scale_fill_gradientn(colours = myPalette(100), name="Probability")
    sg_plot = sg_plot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
    sg_plot = sg_plot + geom_abline(intercept=0, slope=1)
    sg_plot = sg_plot + theme(legend.text=element_text(size=5),
                              legend.title=element_text(size=lts),
                              plot.title=element_text(size=12))


    pdf(paste(save, ".pdf", sep=""), width=8, height=6)
    multiplot(growth_plot, sg_plot, infection_plot, recovery_plot, cols=2)
    dev.off()

  }

  return(list(P=P, S=S, R=R, G=G))

}

get_full_P = function(P, row1, col1, min_size, y, plot_it=T){

  # Function takes in the fungal kernel with a zero class and sticks
  # the zero row and column onto it

  # Parameters
  # ----------
  # P : the transition matrix without a zero class
  # row1 : vector specifying the probability of any stage transition to 0
  # col1 : vector specifying the of zero transitioning to any stage
  # min_size : Lower bound of the kernels
  # y : midpoints of each stage, used for plotting

  # Returns
  # -------
  # : matrix
  #   The transition matrix with a 0 class

  full_P = matrix(NA, nrow=dim(P)[1] + 1, ncol=dim(P)[2] + 1)
  full_P[1, ] = row1
  full_P[ , 1] = col1
  full_P[2:dim(full_P)[1], 2:dim(full_P)[2]] = P

  # if(plot_it){

  #   par(mfrow=c(1, 1))

  #   image.plot(c(min_size - 1, y), c(min_size - 1, y), t(full_P),
  #     main="Full transition matrix", xlab="Load at t", ylab="Load at t + 1")

  #   abline(0,  1, lwd=3)

  # }

  return(full_P)

}



absorption_times = function(P){

  # Calculate the mean time to absoprtion (death) given a transition matrix P
  # Calculate the variance as well

  # Parameters
  # ----------
  # P : transition matrix
  #
  # Returns
  # -------
  # list: mean=mean absorption time, var=var absorption time

  I = diag(dim(P)[1])

  # Calculate the fundamental matrix
  N = solve((I - P))

  # Exepected time to absorption starting in transient state x
  mean_absorb_time = colSums(N)
  var_absorb_time = colSums((2*(N %*% N) - N)) - mean_absorb_time^2

  return(list(mean=mean_absorb_time, var=var_absorb_time))

}

get_elasticity_and_sensitivity = function(P, h, y, plot_it=T, save="temp"){

  # Calculate the elasticities and sensitivities of the population growth rate
  # to transition probabilities in the matrix P

  # Parameters
  # ----------
  # P : transition matrix
  # h : cell width
  # y : midpoints of cells, more plotting

  # Returns
  # -------
  # list: elas= elasticity matrix, sens= sensitivity matrix


  lam=Re(eigen(P)$values[1])
  w.eigen = Re(eigen(P)$vectors[,1])
  stable.dist = w.eigen / sum(w.eigen)
  v.eigen = Re(eigen(t(P))$vectors[,1])
  repro.val = v.eigen / v.eigen[1]

  v.dot.w = sum(stable.dist*repro.val) * h

  sens = sensitivity(P)# outer(repro.val,stable.dist) / v.dot.w
  elas = elasticity(P)# matrix(as.vector(sens)*as.vector(P) / lam, nrow=length(y))

  if(plot_it){

    myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

    long_elas = format_for_heat(t(elas), y, y, c('size', 'sizeNext', 'value'))
    long_elas$id = "Elasticity"

    # Just plotting elasticity
    zplot = ggplot(long_elas, aes(x=size, y=sizeNext, fill=value)) + geom_tile()
    zplot = zplot + scale_fill_gradientn(colours = myPalette(100), name="Elasticity")
    zplot = zplot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
    zplot = zplot + facet_wrap(~id)

    ggsave(paste(save, ".pdf", sep=""), height=5, width=7)

  }

  return(list(elas=elas, sens=sens))

}

get_stable_dist = function(P){

  # Get the stable distribution of a P matrix based on Caswell 2001

  w.eigen = Re(eigen(P)$vectors[,1])
  stable_dist = w.eigen / sum(w.eigen)
  return(stable_dist)

}

get_vm_ratio = function(P, y){

  # Get the variance to mean ratio of the stable distribution given the full
  # transition matrix P. y is the step size

  stable_dist = get_stable_dist(P)

  cond_stable_dist = stable_dist[2:length(stable_dist)] / (1 - stable_dist[1])
  cond_mean = sum(cond_stable_dist * y)
  cond_var = sum(cond_stable_dist * y^2) - (cond_mean)^2

  return(cond_var / cond_mean)

}

##############################################################################

g0_x = function(x, params){
  # Define the initial infection function
  return(dnorm(x, mean=params$clump_mean, sd=params$clump_sd))

}

s_x = function(x, params) {

  # Define survival function
  u = exp(params$surv_int + params$surv_size * x)
  return(u / (1 + u))

}


r_x = function(x, params) {

  # Define the loss function, r(x) (also defined as l(x))

  u = exp(params$loss_int + params$loss_size*x + params$loss_temp)
  probs = u / (1 + u)

  # Truncate loss probability. Above this load no loss can occur
  probs[x > params$loss_trunc] = 0
  return(probs)

}


k_xpx = function(xp, x, params){

  # Define the full kernel (for eviction calculations)

  return(g_xpx(xp, x, params) * s_x(x, params) * (1 - r_x(x, params)))

}


g_xpx = function(xp, x, params){

    # Tranistion from x to xp given some parameters

    if(params$linear){ # Linear growth function

      return(dnorm(xp, mean=params$growth_int + params$growth_size * x +
         params$growth_temp,
         sd=sqrt(params$growth_sigma2*exp(2 * params$growth_sigma2_exp * x))))

    } else{ # Non-linear growth function

      x[x <= params$growth_deriv0] = params$growth_deriv0
      return(dnorm(xp, mean=params$growth_int + params$growth_size * x +
         params$growth_temp + params$growth_size_sq * x**2,
         sd=sqrt(params$growth_sigma2*exp(2 * params$growth_sigma2_exp * x))))
    }
  }


set_temp_params = function(desired_temp, params, linear){

  # Specific parameters for each temperature, given a desired temp
  #
  # Parameters
  # ----------
  # desired_temp: int, the desired temperature at which to calculate params
  # params: list, a list-like object with the base params
  # linear: bool, TRUE is you are using a linear growth function. FALSE otherwise
  #
  # Returns
  # --------
  # : Updated parameters

  tparams = data.frame(holder=NA)

  # Two parameters for survival fxn, same for everybody
  tparams$surv_int = params$surv_int
  tparams$surv_size = params$surv_slope

  if(linear){ # Assuming a linear effect of size

    # 5 parameters for growth function
    tparams$growth_int = params$growth_int
    tparams$growth_temp = params$growth_temp * desired_temp
    tparams$growth_size = params$growth_size
    tparams$growth_sigma2 = params$growth_sigma2
    tparams$growth_sigma2_exp = params$growth_sigma2_exp

    # Three parameters for the loss function
    tparams$loss_int = params$loss_int
    tparams$loss_size = params$loss_size
    tparams$loss_temp = params$loss_temp * desired_temp

    # Probability of becoming infected
    u = exp(params$prob_inf_int + params$prob_inf_temp * desired_temp)
    tparams$prob_inf = u / (1 + u)

    # Clumped distribution
    tparams$clump_mean = params$clump_int + params$clump_temp * desired_temp
    tparams$clump_sd = sqrt(params$clump_sigma2*exp(2 * params$clump_sigma2_exp * desired_temp))

    # Class zero survival
    tparams$class_zero_surv = 1
    tparams$linear = TRUE

  } else{ # Assuming a non-linear effect of size

      # Parameters for growth function
    tparams$growth_int = params$growth_int
    tparams$growth_temp = params$growth_temp * desired_temp
    tparams$growth_size = params$growth_size + params$growth_size_temp_int * desired_temp
    tparams$growth_size_sq = params$growth_size_sq
    tparams$growth_sigma2 = params$growth_sigma2
    tparams$growth_sigma2_exp = params$growth_sigma2_exp

    # Find where the derivative is 0
    tparams$growth_deriv0 = (-1 * tparams$growth_size) / (2 * tparams$growth_size_sq)

    # Three parameters for the loss function
    tparams$loss_int = params$loss_int
    tparams$loss_size = params$loss_size
    tparams$loss_temp = params$loss_temp * desired_temp

    # Probability of becoming infected
    u = exp(params$prob_inf_int + params$prob_inf_temp * desired_temp)
    tparams$prob_inf = u / (1 + u)

    # Clumped distribution
    tparams$clump_mean = params$clump_int + params$clump_temp * desired_temp
    tparams$clump_sd = sqrt(params$clump_sigma2*exp(2 * params$clump_sigma2_exp * desired_temp))

    # Class zero survival
    tparams$class_zero_surv = 1
    tparams$linear = FALSE

  }

  return(tparams)
}


### Various transmission functions ###

density_dependent_zoospores = function(n_t, beta, params){

    # Function specifying the probability of infection given the total
    # number zoospores in the population
    #
    # Parameters
    # n_t : Population vector
    # beta : transmission coefficient
    # y : midpoints for IPM (different IPM loads)
    # params : IPM parameters

    # Parameter specifying baseline probability of infection
    a = -1 * log(1 - params$prob_inf)

    # This is funky...need to count total zoospores...un-log transform and then
    # re-transform + 1...what effect is this having?

    total_zoospores_t = log(sum(exp(y) * n_t[2:(bins + 1)]) + 1)

    # Calculate the infection probability
    infection_prob = 1 - exp(-(a + params$beta * total_zoospores_t))

    return(infection_prob)
}

density_dependent_individuals = function(n_t, y, params){

    # Function specifying the probability of infection given the total
    # number zoospores in the population
    #
    # Parameters
    # n_t : Population vector
    # beta : transmission coefficient
    # params : IPM parameters

    # Parameter specifying baseline probability of infection
    a = -1 * log(1 - params$prob_inf)

    total_inf_individuals_t = sum(n_t[2:(bins + 1)])

    # Calculate the infection probability
    infection_prob = 1 - exp(-(a + params$beta * total_inf_individuals_t))

    return(infection_prob)
}

###########################################################################

within_season_simulation = function(n_0, params, bins, bnd, y, h, kernel){

  # Calculates a one step projection given a vector n_0
  #
  # Parameters
  # ----------
  # n_0 : initial population vector
  # params : parameters dataframe
  # bins, bnd, y, h : discretization parameters
  # kernel : density independent IPM kernel

  # Returns
  # -------
  # list of updated pop vector and infection probability

  n_t = as.matrix(n_0, nrow=bins + 1, ncol=1)

  # Calculate the time variant portion of the kernel. This will be updated
  # each time step based on the distribution of zoospores.  Initiall assuming
  # density dependent transmission. The total number of zoospores across all
  # frogs.

  infection_prob = density_dependent_zoospores(n_t, y, params)

  # The zero -> zero transition
  zero_trans = params$class_zero_surv * (1 - infection_prob)

  # The zero to nonzero transitions: Transitioning from zero to infected
  prob_clump = h * dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * infection_prob * prob_clump
  col1 = c(zero_trans, ztnz)

  # Compute the infected to uninfected row. This is just the probability of
  # surviving with load $x$ ($s(x)$) and then losing an infection $r(x)$.  We
  # just need to multiply $s(x) r(x)$

  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)

  # Now we just need to stick these extra columns onto the **P** matrix. The
  # most likely thing that is going to happen is you are going to stay in the
  # uninfected class.

  full_P = get_full_P(P, row1, col1, min_size, y, plot_it=F)

  n_t_new = full_P %*% n_t

  return(list(nt=n_t_new, inf_prob=infection_prob))

}

###########################################################################

build_full_P = function(min_size, max_size, bins, params) {

    # Function that builds the full transition matrix for the host parasite
    # IPM

    # Parameters
    # min_size: float, lower bound of the matrix
    # max_size: float, upper bound of the matrix
    # bins: int, number of bins in the matrix
    # params: list, list of parameters to use to build P

    # Returns
    # full_P : full transition matrix

    # Get the relevant matrix parameters: midpoints, cell edges, and cell width
    matrix_params = set_discretized_values(min_size, max_size, bins)
    bnd = matrix_params$bnd
    y = matrix_params$y
    h = matrix_params$h

    # Calculate and plot eviction.  Eviction is most important when e is large and. Using the script provided by Williams et al. for checked for eviction

    # Get the kernel of the IPM without the zero class
    kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=F)
    P = kernels$P

    # The zero -> zero transition
    zero_trans = params$class_zero_surv * (1 - params$prob_inf)

    # The zero to nonzero transitions: Transitioning to zero to infected
    prob_clump = h * dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
    ztnz = params$class_zero_surv * params$prob_inf * prob_clump
    col1 = c(zero_trans, ztnz)

    # The loss probabilities: Transitioning from non-zero to zero
    nztz = kernels$S * kernels$R
    row1 = c(zero_trans, nztz)

    # Now we just need to stick these extra columns onto the **P** matrix. The
    # most likely thing that is going to happen is you are going to stay in the
    # uninfected class.

    full_P = get_full_P(P, row1, col1, min_size, y, plot_it=F)

    return(full_P)

}

###########################################################################

monte_carlo_params = function(params){
  # For the various vital functions in the IPM, uses the vcov matrix for each
  # set of params to draw new params (assuming multivariate normality under
  # asymptotic likelihood theory) and return a params matrix with the updated
  # params

  # Parameters
  # params : list, properly formatted list of parameters from ipm_vital*.R

  # Returns
  # rand_params : list, list where growth, loss, surv, clump, and prob_inf
  # parameters are random draws

  rand_params = params
  vital_names = c('growth', 'loss', 'surv', 'clump', 'prob_inf')

  for(vital_name in vital_names){

    # Draw params for each vital function
    mean_vect = paste(vital_name, "coef", sep="_")
    vcov_mat = paste(vital_name, "vcov", sep="_")
    #print(paste(mean_vect, vcov_mat))
    vital_draw = mvrnorm(n=1, mu=params[[mean_vect]], Sigma=params[[vcov_mat]])
    vital_draw = as.list(vital_draw)

    # Store params if they exist
    rand_params[[paste(vital_name, "int", sep="_")]] = vital_draw[["(Intercept)"]]
    rand_params[[paste(vital_name, "temp", sep="_")]] = vital_draw[["temp"]]
    rand_params[[paste(vital_name, "size", sep="_")]] = vital_draw[["size"]]
    rand_params[[paste(vital_name, "size_sq", sep="_")]] = vital_draw[["I(size^2)"]]
    rand_params[[paste(vital_name, "size_temp_int", sep="_")]] = vital_draw[["size:temp"]]

  }

  return(rand_params)

}

################################

compute_quartiles = function(est_matrix){
  # Computes the median, lower quartile, and upper quartile of est_matrix
  # column wise

  median = apply(est_matrix, 2, median)
  upper = apply(est_matrix, 2, uq)
  lower = apply(est_matrix, 2, lq)

  return(list(lower=lower, median=median, upper=upper))

}

