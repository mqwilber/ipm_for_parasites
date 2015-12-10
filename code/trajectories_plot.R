## Description
## -----------
## Script for making trajectory plot found in the supplementary material. Loads
## in the temperature data and makes a plot of Bd load trajectories for each
## individual across temperatures

## Notes: Assumes that the working directory is ipm_models/code

## Author: Mark Wilber

library(ggplot2)
library(plyr)

# Load in temperature data
temp_dat = read.csv("../data/archival/TempExp.csv")
temp_dat_adult = temp_dat[temp_dat$InitialStage == "adult", ]

# Drop the dead individuals and the 250 days measures
no_dead = temp_dat_adult[temp_dat_adult$ZE != "dead", ]
no_dead = no_dead[no_dead$day < 200, ]
no_dead = no_dead[no_dead$Temperature != 26, ]

# Convert factor to float
no_dead$ZE = as.numeric(levels(no_dead$ZE))[no_dead$ZE]

# Rearrage temperature factors and rename
no_dead$Temperature = factor(no_dead$Temperature, levels=c(4, 12, 20), ordered=T)
no_dead$Temperature = revalue(no_dead$Temperature, c("4"="4 C", "12"="12 C", "20"="20 C"))

split_vect = function(data){strsplit(as.character(str_val), "_")[[1]][2]}

splits = character()
for(str_val in no_dead$individual){
    splits = c(splits, split_vect(str_val))
}

no_dead$ind_group = splits

# Visualize the different trajectories of individuals
ggplot(data=no_dead, aes(x=day, y=log(ZE + 1))) +
                    geom_line(aes(color=ind_group), size=0.1, alpha=0.7) +
                    geom_point(aes(color=ind_group), size=1) + theme_bw() +
                    xlab("Days") + ylab("Log zoospore equivalents + 1") +
                    labs(color="Individual") +
                    facet_wrap(~Temperature, nrow=3, ncol=1)

ggsave("../results/trajectories.pdf", width=6,
    height=7)

