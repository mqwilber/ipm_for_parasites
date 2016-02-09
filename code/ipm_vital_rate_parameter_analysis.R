## Description
## -----------
## This script fits the various vital rate functions used in the host-parasite
## IPM described in the manuscript. In short this script fits a,
## variety of statistical models in order to get the parameters for the IPM
## model. This script makes diagnostics plots, compares models, etc. and saves
## the best fit parameter values for each vital rate function.

## See the main manuscript for the model that we are fitting.  In short it
## requires 5 functions
## ---------------------
## 1. A survival function (s(x))
## 2. A growth function (G(x', x))
## 3. A loss of infection function (l(x))
## 4. A transmission function (phi(I(x, t)))
## 5. A initial infection burden function (G_0(x'))

## Notes: This script should be executed with ipm_models/code as the working
## directory

## Author: Mark Wilber

library(AICcmodavg)
library(ggplot2)
library(boot)
library(nlme)
library(lme4)
library(plyr)

# Load in temperature data
temp_data = read.csv("../data/formatted/converted_temp_data.csv")

# Log-transform all Bd loads but zeros so that they can techinically range from
# -inf to inf.  This will be important when we are dealing with eviction
# problems later on

log_trans_data = function(data){

  # Function to log-transform all data Bd-data but not the zeros

  nz_ind = data$size != 0
  data$size[nz_ind] = log(data$size[nz_ind])
  nznext_ind = (data$sizeNext != 0) & !(is.na(data$sizeNext))
  data$sizeNext[nznext_ind] = log(data$sizeNext[nznext_ind])

  return(data)

}

temp_data = log_trans_data(temp_data)

# The time between successive observations on a frog in 3 days.  A few time
# points were longer to let's exclude those

days = 3
temp_data = temp_data[temp_data$time == 3, ]

# A set of data without the temperature 26
temp_data_no_26 = temp_data[temp_data$temp != 26, ]

##############################################################################

# Set up a data.frame to hold the parameters of the vital functions
params_no_26_l=list() # Holds linear parameters
params_no_26_nl=list() # Hold nonlinear parameters

drop_zeros = function(data){

  # Function to remove zeros from size and sizeNext

  data_nozeros = data[data$size != 0, ]

  # Drop individuals that lose infection as well as we are considering this a separate process
  ind = data_nozeros$sizeNext != 0
  ind[is.na(ind)] = TRUE
  data_nozeros = data_nozeros[ind, ]
  return(data_nozeros)

}

temp_data_nozeros = drop_zeros(temp_data)
temp_data_no_26_nozeros = drop_zeros(temp_data_no_26)

# Plot growth kernel
ggplot(data=temp_data_no_26_nozeros, aes(x=size, y=sizeNext)) +
            geom_point(aes(color=as.factor(temp)))

##############################################################################

# Fit the survival function $s(x)$

# Fit the glm and set the parameters. I have transformed size to log(size)

# Because there is not death in 4 or 12, let's only fit the data with 20
# degrees
survival2 = glm(surv ~ size, family=binomial,
            data=temp_data_no_26_nozeros[temp_data_no_26_nozeros$temp == 20, ])

# Fit with all temperature points included
survival3 = glm(surv ~ size, family=binomial, data=temp_data_no_26_nozeros)


vals = seq(-5, 18, len=100)
mod_pred2 = predict.glm(survival2, newdata=data.frame(size=vals), type="response")
mod_pred3 = predict.glm(survival3, newdata=data.frame(size=vals), type="response")

# Plot comparing the two estimated survival functions
comp_data = data.frame(list(size=c(vals, vals), preds=c(mod_pred2, mod_pred3),
              mods=rep(c("Only 20 C", "All temperatures"), c(100, 100))))
ggplot(data=comp_data, aes(x=size, y=preds)) + geom_line(aes(color=mods)) +
              theme_bw() + geom_vline(xintercept=log(10000), linetype=2) +
              xlab("Log zoospore load at time t") +
              ylab("Probability of survival, s(x)") +
              labs(color = "Data set used")

ggsave("../results/compare_survival_fxns.pdf", width=6, height=5)

plot(vals, mod_pred2)
abline(v=log(10000))

# Save the parameters
params_no_26_l$surv_int = coefficients(survival2)[1]
params_no_26_l$surv_slope = coefficients(survival2)[2]
params_no_26_l$surv_vcov = vcov(survival2)
params_no_26_l$surv_coef = coefficients(survival2)

params_no_26_nl$surv_int = coefficients(survival2)[1]
params_no_26_nl$surv_slope = coefficients(survival2)[2]
params_no_26_nl$surv_vcov = vcov(survival2)
params_no_26_nl$surv_coef = coefficients(survival2)

# Plot survival function
rep_vals = c('4 C', '12 C', '20 C')
names(rep_vals) = c(4, 12, 20)

temp_data_no_26_nozeros$temp_fact = revalue(as.factor(temp_data_no_26_nozeros$temp), rep_vals)

ggplot(data=temp_data_no_26_nozeros, aes(x=size, y=surv)) +
  geom_point(aes(color=temp_fact), alpha=0.7) + scale_shape(solid=F) +
  stat_smooth(data=subset(temp_data_no_26_nozeros, temp == 20),
              method=glm, family=binomial, color="black", alpha=0.1) +
  theme_bw() + scale_colour_brewer(palette = "Set2") +
  labs(color = "Temperature, C") +
  xlab("Log zoospore load at time t") +
  ylab("Probability of survival, s(x)") +
  scale_y_continuous(breaks=round(seq(0, 1, len=6), 1)) +
  facet_wrap(~temp_fact, nrow=3) +
  scale_colour_discrete(guide = FALSE)

ggsave("../results/survival_fxn.pdf", width=4, height=5)

ggplot(data=temp_data_no_26_nozeros, aes(x=size, y=surv)) +
  #stat_smooth(method=glm, family=binomial) +
  geom_jitter(position=position_jitter(height=0.1), aes(color=as.factor(temp))) + scale_shape(solid=F) +
  theme_bw() + scale_color_grey() +
  labs(color = "Temperature, C") +
  xlab("Log zoospore load at time t") +
  ylab("Probability of survival, s(x)")

ggsave("../results/survival_fxn_points.pdf", width=6, height=5)

##############################################################################

# Fit the growth function, $G(x', x)$

# Drop the NAs in the dataset
temp_data_nozeros_nona = temp_data_nozeros[!is.na(temp_data_nozeros$sizeNext), ]
temp_data_no_26_nozeros_nona = temp_data_no_26_nozeros[!is.na(temp_data_no_26_nozeros$sizeNext), ]

# Fit two standard linear models with temp as a factor and as a continuous
# variable
growth1 = glm(sizeNext ~ size + as.factor(temp), na.action=na.omit,
                data=temp_data_nozeros_nona)
growth2 = glm(sizeNext ~ size + temp, na.action=na.omit,
                data=temp_data_nozeros_nona)
growth3 = glm(sizeNext ~ size*as.factor(temp), na.action=na.omit,
                data=temp_data_nozeros_nona)
# Compare the two models with ANOVA
anova(growth1, growth2, test='F')
anova(growth1, growth3, test='F')

# Check out a diagnostic plot...looks pretty good
glm.diag.plots(growth1)

# Let's also consider more some more complex models in which we have a
# either a random effect of individual or non-constant error terms

# Same model at growth1
growth1_0 = lm(sizeNext ~ size, data=temp_data_no_26_nozeros)

growth1_1 = gls(sizeNext ~ size + as.factor(temp), na.action=na.omit,
                    data=temp_data_nozeros_nona)

# Random effect of frog
growth1_2 = lme(sizeNext ~ size + as.factor(temp), random= ~1 | individual,
                    na.action=na.omit, data=temp_data_nozeros_nona)

# Variance a function of size and temperature
growth1_3 = gls(sizeNext ~ size + as.factor(temp),
                    weights=varExp(form=~size | as.factor(temp)),
                    na.action=na.omit, data=temp_data_nozeros_nona)

growth1_4 = lme(sizeNext ~ size + as.factor(temp), random= ~1 | individual,
                    weights=varExp(form=~size | as.factor(temp)),
                    na.action=na.omit, data=temp_data_nozeros_nona)

growth1_5 = gls(sizeNext ~ as.factor(temp) + size,
                    weights=varExp(form=~size | as.factor(temp)),
                    na.action=na.omit, data=temp_data_nozeros_nona, method="ML")

growth1_6 = gls(sizeNext ~ as.factor(temp),
                    weights=varExp(form=~size | as.factor(temp)),
                    na.action=na.omit, data=temp_data_nozeros_nona, method="ML")

growth1_7 = gls(sizeNext ~ size,
                    weights=varExp(form=~size | as.factor(temp)),
                    na.action=na.omit, data=temp_data_nozeros_nona, method="ML")

# If I drop temperature 26, is a linear relationship with temp appropriate?
growth1_8 = gls(sizeNext ~ temp + size,
                    weights=varComb(varExp(form=~size), varExp(form=~temp)),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_9 = gls(sizeNext ~ temp + size,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_10 = gls(sizeNext ~ temp + size, na.action=na.omit,
                      data=temp_data_no_26_nozeros_nona, method="ML")

growth1_11 = gls(sizeNext ~ temp,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_12 = gls(sizeNext ~ size,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_13 = lme(sizeNext ~ size + temp, random= ~1 | individual,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_14 = gls(sizeNext ~ size*temp + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_15 = gls(sizeNext ~ size*temp + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_15b = gls(sizeNext ~ size*temp + I(size^2),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_15c = gls(sizeNext ~ size*temp + I(size^2),
                    weights=varExp(form=~size*temp),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_15d = gls(sizeNext ~ size*temp,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_16 = gls(sizeNext ~ as.factor(temp) + size,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_17 = gls(sizeNext ~ temp + size + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")

growth1_18 = gls(sizeNext ~ as.factor(temp) + size + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona, method="ML")


anova(growth1_8, growth1_9, growth1_10, growth1_11, growth1_12, growth1_13,
        growth1_15, growth1_16, growth1_17, growth1_18, growth1_15b, growth1_15c, growth1_15d)
anova(growth1_13, growth1_14) # Random effect of individual needed

# Growth1_15 is the best model when the temp 26 is dropped based on AIC.
# A slightly non-linear effect of size and an interaction between size
# and temperature. Growth1_9 is the best linear model

# The problem now it the following: Because of this non-linearity I think I will
# need to think about the Bd growth function as a piece-wise thing.
# Given a certain temperature, below a certain size let's say 0 (log(1) ZE)
# we are going to assume that the Bd growth kernel is constant.

# Similarly, above a certain size

# Refit Growth1_9 and Growth1_15 with REML
growth1_9 = gls(sizeNext ~ temp + size,
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona,
                    method="REML")

growth1_15 = gls(sizeNext ~ size*temp + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp_data_no_26_nozeros_nona,
                    method="REML")

# What happens if I drop 4 degrees? Is temp still important? Sure is.
temp4 = temp_data_no_26_nozeros[temp_data_no_26_nozeros$temp != 4, ]
growth_no_4_1 = gls(sizeNext ~ as.factor(temp) + size + I(size^2),
                    weights=varExp(form=~size),
                    na.action=na.omit, data=temp4,
                    method="ML")

# There is still a strong positive effect of temperature.

# Make some diagnostic plots of growth1_9
pdf("../results/growth_model_diagnostics_linear.pdf", width=8, height=4)

par(mfrow=c(1, 2))
plot(predict(growth1_9), resid(growth1_9, type="normalized"),
        ylab="Residuals", xlab="Fitted", main="Residual Plot")
lines(lowess(predict(growth1_9), resid(growth1_9, type="normalized")),
          lty=2, col="red")
qqnorm(resid(growth1_9, type="normalized"))
qqline(resid(growth1_9, type="normalized"))

dev.off()


# Make some diagnostic plots of growth1_15
pdf("../results/growth_model_diagnostics_nonlinear.pdf", width=8, height=4)

par(mfrow=c(1, 2))
plot(predict(growth1_15), resid(growth1_15, type="normalized"),
        ylab="Residuals", xlab="Fitted", main="Residual Plot")
lines(lowess(predict(growth1_15), resid(growth1_15, type="normalized")),
          lty=2, col="red")
qqnorm(resid(growth1_15, type="normalized"))
qqline(resid(growth1_15, type="normalized"))

dev.off()

# Reset par
par(mfrow=c(1, 1))

### Predicted plot for different temperatures ###

# Plots for linear growth function
sizes = seq(-1, 12, len=100)
temps = c(4, 12, 20)
temp_holder = numeric()
temp_se_holder = numeric()

for(temp in temps){
    preds = predictSE.gls(growth1_9, newdata=data.frame(size=sizes, temp=temp))
    temp_holder = c(temp_holder, preds$fit)
    temp_se_holder = c(temp_se_holder, preds$se)
}

temps_factor = as.factor(rep(temps, each=100))
temp_plot = data.frame(list(size=rep(sizes, 3), sizeNext=temp_holder,
  temp=temps_factor, upper=temp_holder + temp_se_holder,
  lower=temp_holder - temp_se_holder))

tplot = ggplot(data=temp_data_no_26_nozeros, aes(x=size, y=sizeNext)) +
            geom_point(aes(color=as.factor(temp)), alpha=0.7)
tplot = tplot + geom_line(data=temp_plot, aes(x=size, y=sizeNext, color=temp))
tplot = tplot + geom_ribbon(data=temp_plot, aes(ymin=lower, ymax=upper, color=temp), linetype=0, alpha=0.05)
tplot = tplot + theme_bw()
#tplot = tplot + geom_line(data=temp_plot, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + labs(color="Temperature, C") +
                xlab("Log zoospore load at time t") +
                ylab("Log zoospore load at time t + 1")
#tplot = tplot + scale_colour_brewer(palette = "Set2")
#tplot = tplot + facet_wrap(~temp)
tplot

ggsave("../results/empirical_growth.pdf", width=6, height=5)


tplot = ggplot(data=temp_plot, aes(x=size, y=sizeNext))
tplot = tplot + geom_point(data=temp_data_no_26_nozeros, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + theme_bw()
#tplot = tplot + geom_line(data=temp_plot, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + labs(color="Temperature, C") +
                xlab("Log zoospore load at time t") +
                ylab("Log zoospore load at time t + 1")

ggsave("../results/empirical_growth_points.pdf", width=6, height=5)

# Plots for nonlinear growth function

sizes = seq(-1, 12, len=100)
temps = c(4, 12, 20)
temp_holder = numeric()
temp_se_holder = numeric()

for(temp in temps){
    preds = predictSE.gls(growth1_15, newdata=data.frame(size=sizes, temp=temp))
    temp_holder = c(temp_holder, preds$fit)
    temp_se_holder = c(temp_se_holder, preds$se)
}

temps_factor = rep(temps, each=100)
temp_plot = data.frame(list(size=rep(sizes, 3), sizeNext=temp_holder,
  temp=temps_factor, upper=temp_holder + temp_se_holder,
  lower=temp_holder - temp_se_holder))
tplot = ggplot(data=temp_plot, aes(x=size, y=sizeNext))
tplot = tplot + geom_point(data=temp_data_no_26_nozeros, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + geom_line(aes(color=as.factor(temp)))
tplot = tplot + geom_ribbon(aes(ymin=lower, ymax=upper, color=as.factor(temp)), linetype=0, alpha=0.1)
tplot = tplot + theme_bw()
#tplot = tplot + geom_line(data=temp_plot, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + labs(color="Temperature, C") +
                xlab("Log zoospore load at time t") +
                ylab("Log zoospore load at time t + 1")
tplot

ggsave("../results/empirical_growth_nonlinear.pdf", width=6, height=5)


tplot = ggplot(data=temp_plot, aes(x=size, y=sizeNext))
tplot = tplot + geom_point(data=temp_data_no_26_nozeros, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + theme_bw()
#tplot = tplot + geom_line(data=temp_plot, aes(x=size, y=sizeNext, color=as.factor(temp)))
tplot = tplot + labs(color="Temperature, C") +
                xlab("Log zoospore load at time t") +
                ylab("Log zoospore load at time t + 1")

ggsave("../results/empirical_growth_points_nonlinear.pdf", width=6, height=5)

# Save the parameters of growth1_9; when temperature 26 is dropped

params_no_26_l$growth_int = coefficients(growth1_9)[1]
params_no_26_l$growth_temp = coefficients(growth1_9)[2]
params_no_26_l$growth_size = coefficients(growth1_9)[3]
params_no_26_l$growth_vcov = growth1_9$varBeta # VarCovar matrix for parameters
params_no_26_l$growth_coef = coefficients(growth1_9) # VarCovar matrix for parameters

#  Get the variance structure
params_no_26_l$growth_sigma2 = summary(growth1_9)$sigma^2
params_no_26_l$growth_sigma2_exp = as.numeric(growth1_9$modelStruct$varStruct)[1]

# Save the parameters of growth1_15; when temperature 26 is dropped
params_no_26_nl$growth_int = coefficients(growth1_15)[1]
params_no_26_nl$growth_size = coefficients(growth1_15)[2]
params_no_26_nl$growth_temp = coefficients(growth1_15)[3]
params_no_26_nl$growth_size_sq = coefficients(growth1_15)[4]
params_no_26_nl$growth_size_temp_int = coefficients(growth1_15)[5]
params_no_26_nl$growth_vcov = growth1_15$varBeta # VarCovar matrix for parameters
params_no_26_nl$growth_coef = coefficients(growth1_15) # VarCovar matrix for parameters

#  Get the variance structure
params_no_26_nl$growth_sigma2 = summary(growth1_15)$sigma^2
params_no_26_nl$growth_sigma2_exp = as.numeric(growth1_15$modelStruct$varStruct)[1]


##############################################################################

# Calculate $p_{inf}$ which is the probability of having no infection to
# becoming infected

# This is a function of temperature, similar to growth

# Subset the data to only include the alive individuals
alive_individs = function(data){

  alive = data[data$surv == 1, ]
  zero_data = alive[alive$size == 0, ]

  # Make a new column
  zero_data$transition = 0
  zero_data$transition[zero_data$size == 0 & zero_data$sizeNext != 0] = 1

  return(zero_data)
}

zero_data = alive_individs(temp_data)
zero_data_no_26 = alive_individs(temp_data_no_26)

# Fit a logistic model where the probability of transition is a function of
# temperature

p_infection1_1 = glm(transition ~ temp, family=binomial, data=zero_data)
p_infection1_2 = glm(transition ~ as.factor(temp), family=binomial,
                                  data=zero_data)

p_infection1_3 = glm(transition ~ temp, family=binomial, data=zero_data_no_26)
p_infection1_4 = glm(transition ~ as.factor(temp), family=binomial,
                                  data=zero_data_no_26)

anova(p_infection1_3, p_infection1_4, test="Chisq")

# There seems to be a linear effect of temperature when 26.
# More specifically, there is not compelling evidence to use a more complex model

# Set infection probability as a function of temperature
params_no_26_l$prob_inf_int = coefficients(p_infection1_3)[1]
params_no_26_l$prob_inf_temp = coefficients(p_infection1_3)[2]
params_no_26_l$prob_inf_vcov = vcov(p_infection1_3)
params_no_26_l$prob_inf_coef = coefficients(p_infection1_3)

params_no_26_nl$prob_inf_int = coefficients(p_infection1_3)[1]
params_no_26_nl$prob_inf_temp = coefficients(p_infection1_3)[2]
params_no_26_nl$prob_inf_vcov = vcov(p_infection1_3)
params_no_26_nl$prob_inf_coef = coefficients(p_infection1_3)

##############################################################################

# Calculate the initial infection burden function, $G_0(x')$

alive = temp_data[temp_data$surv == 1, ]
alive_no_26 = temp_data_no_26[temp_data_no_26$surv == 1, ]

transition_fxn = function(data){

  # Find the frogs who went from 0 load to non-zero load
  transition = data[data$size == 0 & data$sizeNext != 0, ]

  ind = (transition$sizeNext > 5)
  print(paste("Dropping", sum(ind), "points"))

  dropped = transition[ind,]

  # Drop values greater than 6 as they are likely PCR error
  transition = transition[!ind, ]

  return(list(transition=transition, dropped=dropped))

}

transition = transition_fxn(alive)$transition
transition_no_26 = transition_fxn(alive_no_26)$transition

# Fit a model that allows variance and mean of clump_distribution to vary by
# temperature
clump_model = gls(sizeNext ~ as.factor(temp), data=transition,
          weights=varIdent(form = ~ 1 | as.factor(temp)))

clump_model2 = gls(sizeNext ~ as.factor(temp), data=transition)

anova(clump_model, clump_model2)
anova(clump_model2)

# Fit a model to the no 26 temp data with variance varying by temperature
clump_model3 = gls(sizeNext ~ temp, data=transition_no_26,
          weights=varExp(form=~temp), method="ML")

clump_model4 = gls(sizeNext ~ temp, data=transition_no_26, method="ML")

clump_model5 = gls(sizeNext ~ as.factor(temp), data=transition_no_26, method="ML")


anova(clump_model3, clump_model4, clump_model5)
anova(clump_model3, clump_model4)

# Refit clump_model3 with REML
clump_model3 = gls(sizeNext ~ temp, data=transition_no_26,
          weights=varExp(form=~temp), method="REML")

# Model diagnostics on from initial infection load

pdf("../results/initial_infection_diagnostics.pdf", width=8, height=4)

par(mfrow=c(1, 2))

plot(predict(clump_model3), resid(clump_model3, type="normalized"),
              xlab="Predicted values", ylab="Normalized residuals",
              main="Residual plot")
qqnorm(resid(clump_model3, type="normalized"))
qqline(resid(clump_model3, type="normalized"))

dev.off()


# Get the parameters for the no 26 models
params_no_26_l$clump_int = coefficients(clump_model3)[1]
params_no_26_l$clump_temp = coefficients(clump_model3)[2]
params_no_26_l$clump_vcov = clump_model3$varBeta
params_no_26_l$clump_coef = coefficients(clump_model3)

# Get the variance relationship
params_no_26_l$clump_sigma2 = summary(clump_model3)$sigma^2
params_no_26_l$clump_sigma2_exp = as.numeric(clump_model3$modelStruct$varStruct)[1]

# Non-linear model
params_no_26_nl$clump_int = coefficients(clump_model3)[1]
params_no_26_nl$clump_temp = coefficients(clump_model3)[2]
params_no_26_nl$clump_vcov = clump_model3$varBeta
params_no_26_nl$clump_coef = coefficients(clump_model3)

# Get the variance relationship
params_no_26_nl$clump_sigma2 = summary(clump_model3)$sigma^2
params_no_26_nl$clump_sigma2_exp = as.numeric(clump_model3$modelStruct$varStruct)[1]

##############################################################################

### Calculate the survival in the zero class, $s(0)$

# Calculate the number of individuals that died after having a load of 0
# Calculate a specific probability for each temperature

zero_data = temp_data[temp_data$size == 0, ]
zero_data_no_26 = temp_data_no_26[temp_data_no_26$size == 0, ]

zero_surv1 = glm(surv ~ as.factor(temp), family=binomial, data=zero_data)

# For now, let's just set the survival probability for no 26 to 1
params_no_26_l$class_zero_surv = 1
params_no_26_nl$class_zero_surv = 1

##############################################################################

### Calculate the loss (recovery) of infection function, $l(x)$

alive_fxn = function(data){
  # Function for preparing data for analysis

  data['loss'] = 0

  # Don't include the zero class
  alive_trun = data[data$size != 0, ]

  alive_trun$loss[alive_trun$sizeNext == 0] = 1

  # Drop cases where loss occurs over 6 ZEs...probably PCR error
  ind = (alive_trun$size > 5 & alive_trun$loss == 1)
  dropped = alive_trun[ind, ]

  print(paste("Dropping", sum(ind), "points"))

  alive_trun = alive_trun[!ind, ]

  return(list(alive_trun=alive_trun, dropped=dropped))

}

alive_trun = alive_fxn(alive)$alive_trun
alive_trun_no_26 = alive_fxn(alive_no_26)$alive_trun


infection_loss1_1 = glm(loss ~ size, data=alive_trun, family=binomial)

infection_loss1_2 = glm(loss ~ size + as.factor(temp), data=alive_trun,
                  family=binomial)

infection_loss1_3 = glm(loss ~ as.factor(temp), data=alive_trun,
                  family=binomial)

infection_loss1_4 = glm(loss ~ size*as.factor(temp), data=alive_trun,
                  family=binomial)


anova(infection_loss1_3, infection_loss1_2, test="Chisq")
anova(infection_loss1_1, infection_loss1_2, test="Chisq")
anova(infection_loss1_2, infection_loss1_4, test="Chisq")

# Losing an infection is relatively independent of infection size...up to
# a certain point

# Fit the loss/recovery without 26

infection_loss1_5 = glm(loss ~ size, data=alive_trun_no_26, family=binomial)

infection_loss1_6 = glm(loss ~ size + temp, data=alive_trun_no_26,
                  family=binomial)

infection_loss1_7 = glm(loss ~ temp, data=alive_trun_no_26,
                  family=binomial)

infection_loss1_8 = glm(loss ~ size*temp, data=alive_trun_no_26,
                  family=binomial)

infection_loss1_9 = glm(loss ~ size + as.factor(temp), data=alive_trun_no_26,
                  family=binomial)

infection_loss1_10 = glm(loss ~ size*temp + I(size^2), data=alive_trun_no_26,
                  family=binomial)

infection_loss1_11 = glm(loss ~ size + temp + I(temp^2), data=alive_trun_no_26,
                  family=binomial)

AIC(infection_loss1_5, infection_loss1_6, infection_loss1_7,
    infection_loss1_8, infection_loss1_9, infection_loss1_10,
    infection_loss1_11)

anova(infection_loss1_6, infection_loss1_11, test="Chisq")

# Recovery depends on both infection size an temperature.
# Model infection_loss1_9 is the best based on AIC, but not based on Chi-squared
# because they are pretty close, I am going to use the linear relationship

# Save the loss parameters
params_no_26_l$loss_int = coefficients(infection_loss1_6)[1]
params_no_26_l$loss_size = coefficients(infection_loss1_6)[2]
params_no_26_l$loss_temp = coefficients(infection_loss1_6)[3]
params_no_26_l$loss_vcov = vcov(infection_loss1_6)
params_no_26_l$loss_coef = coefficients(infection_loss1_6)


# Save the loss parameters for no temp 26
params_no_26_nl$loss_int = coefficients(infection_loss1_6)[1]
params_no_26_nl$loss_size = coefficients(infection_loss1_6)[2]
params_no_26_nl$loss_temp = coefficients(infection_loss1_6)[3]
params_no_26_nl$loss_vcov = vcov(infection_loss1_6)
params_no_26_nl$loss_coef = coefficients(infection_loss1_6)

##############################################################################

# Save all the IPM parameters to use in later analyses

saveRDS(params_no_26_l, "../results/IPM_parameters_no_26_linear.rds")
saveRDS(params_no_26_nl, "../results/IPM_parameters_no_26_nonlinear.rds")
