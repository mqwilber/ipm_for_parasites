## Description
## -----------
## This script plots the data generated in `transmission_simulation.R`

## Notes: Script assumes the working directory is ipm_models/code

# Load necessary packages
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plyr)

temps = c(12, 16, 20)
beta_vals = read.csv("../results/transmission_results/beta_vals.csv")
inf_vals = read.csv("../results/transmission_results/inf_vals.csv")

get_transmission_data = function(temps, file_start){
    # Function for converting simulation data to a plotable form

    all_data = list()
    for(i in 1:length(temps)){

        temp = temps[i]
        # Load the data from
        data = read.csv(paste("../results/transmission_results/", file_start, temp, ".csv", sep=""))

        # Format data to use in ggplot heat map
        nameddata = as.matrix(data)[, 2:dim(data)[2]]
        colnames(nameddata) = inf_vals$inf_vals
        rownames(nameddata) = beta_vals$beta_vals

        longdata = melt(nameddata)
        longdata$value = longdata$value
        colnames(longdata) = c("beta", "inf_prob", "values")
        longdata$temperature = temp

        all_data[[i]] = longdata
    }

    plot_data = all_data[[1]]
    for(i in 2:length(temps)){

        plot_data = rbind(plot_data, all_data[[i]])
    }

    # Convert temps to strings for plotting
    temp_strings = sapply(temps, function(x){paste(x, "C")})
    names(temp_strings) = as.character(temps)

    plot_data$temperature = as.factor(plot_data$temperature)
    plot_data$temperature = revalue(plot_data$temperature, temp_strings)

    return(plot_data)

}

# Plot the prevalence data
file_start = "prop_uninfected_"
plot_data = get_transmission_data(temps, file_start)

# Define palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
es = 9

zp1 = ggplot(plot_data,
              aes(x = inf_prob, y = beta, fill = 1 - values))
zp1 = zp1 + geom_tile()
zp1 = zp1 + scale_fill_gradientn(colours = myPalette(100), name="Prevalence")
#zp1 = zp1 + coord_equal()
zp1 = zp1 + theme_bw() +  theme(axis.text.x  = element_text(size=es), axis.text.y  = element_text(size=es))
zp1 = zp1 + facet_wrap(~ temperature, nrow=3, ncol=1)
zp1 = zp1 + ylab(expression("Transmission parameter, " * beta)) + xlab(expression("Environmental infection probability, " * omega))
print(zp1)
ggsave("../results/prevalence.pdf")


# Plot the population data
file_start = "pop_size_post_summer_"
plot_data = get_transmission_data(temps, file_start)

# Define palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

zp1 = ggplot(plot_data,
              aes(x = inf_prob, y = beta, fill = 1 - (values / 100)))
zp1 = zp1 + geom_tile()
zp1 = zp1 + scale_fill_gradientn(colours = myPalette(100), name="Proportion\npopulation\ndecline")
#zp1 = zp1 + coord_equal()
zp1 = zp1 + theme_bw() +  theme(axis.text.x  = element_text(size=es), axis.text.y  = element_text(size=es))
zp1 = zp1 + facet_wrap(~ temperature, nrow=3, ncol=1)
zp1 = zp1 + ylab(expression("Transmission coefficient, " * beta)) + xlab(expression("Environmental infection probability, " * omega))
print(zp1)
ggsave("../results/pop_decline.pdf")


