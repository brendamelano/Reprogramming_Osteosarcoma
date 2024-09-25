library(car)
library(drc)
library(ggplot2)
library(magrittr)
library(dplyr)


# Reading in csv files
OS384_cis_atr <- read.csv("~/Desktop/OS_IC50_analysis/IC50_data/2020_01_04_15day_holiday_384_742_cis_atr_pf/384_15day_cis_atr - 2023_01_03_384_742_CIS_ATR_PF_15day_IC50.csv", 
                          header = F)


# Renaming the row names
rownames(OS384_cis_atr)[1] <- "Cis_1"
rownames(OS384_cis_atr)[2] <- "Cis_2"
rownames(OS384_cis_atr)[3] <- "Cis_3"

rownames(OS384_cis_atr)[4] <- "Atr_1"
rownames(OS384_cis_atr)[5] <- "Atr_2"
rownames(OS384_cis_atr)[6] <- "Atr_3"


# Making a dataframe for the ATR and CIS experiment
OS384_cis <-OS384_cis_atr[1:3,]

OS384_atr <-OS384_cis_atr[4:6,]



#######  ATR analysis    ###########



# Computing the mean of the triplicates
OS384_atr <- OS384_cis_atr[4:6,]
mean_values <- OS384_atr %>% summarise_if(is.numeric, mean)
sd_values <- OS384_atr %>% summarise_if(is.numeric, sd)

OS384_atr <- rbind(OS384_atr, mean_values, sd_values)

# Naming the row with the mean values "mean"
rownames(OS384_atr)[4] <- 'mean'
rownames(OS384_atr)[5] <- 'sd'


# Making a vector with ATR concentrations
concentrations <- c(.5, 0.167, 0.056, 0.0185, 0.00617284, 0.002057613, 0.000685871, 
                    0.0002286237, 0)


# Changing the column names to the concentrations 
names(OS384_atr) <- concentrations


# Creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS384_atr)[,'mean'],
  sd = t(OS384_atr)[,'sd']  # Add the standard deviation
)


# Fitting a dose-response curve
curved_fit <- drm(
  formula = cell_count ~ conc,
  data = dose_response,
  fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
)


# Printing the summary of the model
summary(curved_fit)


# Storing the coefficients from the dose-response curve
coefs <- setNames(
  curved_fit$coefficients,
  c("hill", "min_value", "max_value", "ec_50")
)


# 
ic_50 <- with(
  as.list(coefs),
  exp(
    log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
  )
)


ic_50


# Note the logged x-axis
IC50_curve_384_ATR <- plot(curved_fit, 
                           main = "OS384 ATR IC50", 
                           xlab = "Concentration (uM)", 
                           ylab = "Fluorescent Intensity", 
                           digits = 3)


# changing the name of the 2nd column to predictions
names(IC50_curve_384_ATR)[2] <- 'predictions'


IC50_curve_384_ATR <- ggplot(IC50_curve_384_ATR, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs(title = "OS384 ATR inhibitor IC50", x = "Concentration (uM)", y = "Fluorescent Intensity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        plot.title = element_text(size = 10))  # Change size of title


IC50_curve_384_ATR


# Save the plot as an SVG file
ggsave("~/Desktop/ATR_IC50_OS384.svg", plot = IC50_curve_384_ATR, device = "svg", width = 2.5, height = 2.5)



#############     CIS     #############



# Computing the mean for the cis treatment values
OS384_cis <- rbind(OS384_cis, OS384_cis %>% summarise_if(is.numeric, mean))
rownames(OS384_cis)[4] <- 'mean'


# Making a vector with the concentrations
concentrations <- c(5, 1.67, 0.56, 0.185, 0.0617284, 0.02057613, 0.00685871, 
                    0.002286237, 0)


# changing the column names to the concentrations 
names(OS384_cis) <- concentrations


# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS384_cis)[,'mean']
)


# 
dose_response <- dose_response[-1,]


curved_fit <- drm(
  formula = cell_count ~ conc,
  data = dose_response,
  fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
)


summary(curved_fit)


coefs <- setNames(
  curved_fit$coefficients,
  c("hill", "min_value", "max_value", "ec_50")
)


ic_50 <- with(
  as.list(coefs),
  exp(
    log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
  )
)


# printing the IC50 value
ic_50


# Note the logged x-axis
IC50_curve_384_cis <- plot(curved_fit, main = "OS384 Cisplatin IC50", xlab = "Concentration (uM)", 
     ylab = "Fluorescent Intensity", digits = 3)


# Changing the name of the 2nd column to predictions
names(IC50_curve_384_cis)[2] <- 'predictions'


# plotting the dos response curve with ggplot2
IC50_curve_384_cis_gg <- ggplot(IC50_curve_384_cis, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs(title = "OS384 Cisplatin inhibitor IC50", x = "Concentration (uM)", y = "Fluorescent Intensity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        plot.title = element_text(size = 10))  # Change size of title


# Save the plot as an SVG file
ggsave("~/Desktop/Cis_IC50_OS384.svg", plot = IC50_curve_384_cis_gg, device = "svg", width = 2.5, height = 2.5)

