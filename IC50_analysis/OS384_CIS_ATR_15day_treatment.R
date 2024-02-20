library(car)
library(drc)
library(ggplot2)
library(magrittr)
library(dplyr)


# reading in csv files
OS384_cis_atr <- read.csv("~/Desktop/IC50_data/2020_01_04_15day_holiday_384_742_cis_atr_pf/384_15day_cis_atr - 2023_01_03_384_742_CIS_ATR_PF_15day_IC50.csv", header = F)


# renaming the row names
rownames(OS384_cis_atr)[1] <- "Cis_1"
rownames(OS384_cis_atr)[2] <- "Cis_2"
rownames(OS384_cis_atr)[3] <- "Cis_3"

rownames(OS384_cis_atr)[4] <- "Atr_1"
rownames(OS384_cis_atr)[5] <- "Atr_2"
rownames(OS384_cis_atr)[6] <- "Atr_3"


# making a dataframe for the ATR experiment
OS384_cis <-OS384_cis_atr[1:3,]

OS384_atr <-OS384_cis_atr[4:6,]



###  ATR analysis    ###


# computing the mean of the triplicates
OS384_atr <- rbind(OS384_atr, OS384_atr %>% summarise_if(is.numeric, mean))


# naming the row with the mean values "mean"
rownames(OS384_atr)[4] <- 'mean'


# making a vector with the concentrations
# atr concentrations
concentrations <- c(.5, 0.167, 0.056, 0.0185, 0.00617284, 0.002057613, 0.000685871, 
                    0.0002286237, 0)


# changing the column names to the concentrations 
names(OS384_atr) <- concentrations


# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS384_atr)[,'mean']
)


# plotting the dose response curve with the dose (x-axis) on a log scale
ggplot(
  data = dose_response, aes(x = conc, y = cell_count)) +
  geom_point() +
  geom_line() +
  labs(x = "Drug Concentration", y = "Cell Count") +
  ylim(0, NA) + scale_x_continuous(trans = "log10")


# fitting a dose-response curve
curved_fit <- drm(
  formula = cell_count ~ conc,
  data = dose_response,
  fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
)


# printing the summary of the model
summary(curved_fit)


# storing the coefficients from the dose-response curve
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
# This also actually creates a dataframe that can be used to plot the same data in ggplot2
IC50_curve_384_ATR <- plot(curved_fit, main = "OS384 ATR inhibitor IC50", xlab = "Concentration (uM)", 
     ylab = "Fluorescent Intensity", digits = 3)


# changing the name of the 2nd column to predictions
names(IC50_curve_384_ATR)[2] <- 'predictions'


# plotting the dos response curve with ggplot2
IC50_curve_384_ATR <- ggplot(IC50_curve_384_ATR, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs(title = "OS384 ATR inhibitor IC50", x = "Concentration (uM)", y = "Fluorescent Intensity") +
  theme_bw()


# Save the plot as an SVG file
ggsave("~/Desktop/ATR_IC50_OS384.svg", plot = IC50_curve_384_ATR, device = "svg")


####### cis #######


OS384_cis <- rbind(OS384_cis, OS384_cis %>% summarise_if(is.numeric, mean))
rownames(OS384_cis)[4] <- 'mean'


# making a vector with the concentrations


concentrations <- c(5, 1.67, 0.56, 0.185, 0.0617284, 0.02057613, 0.00685871, 
                    0.002286237, 0)

# changing the column names to the concentrations 
names(OS384_cis) <- concentrations




# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS384_cis)[,'mean']
)


dose_response <- dose_response[-1,]


# plotting the dose response curve with the dose (x-axis) on a log scale
ggplot(
  data = dose_response, aes(x = conc, y = cell_count)) +
  geom_point() +
  geom_line() +
  labs(x = "Drug Concentration", y = "Cell Count") +
  ylim(0, NA) + scale_x_continuous(trans = "log10")



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
# Note the logged x-axis
IC50_curve_384_cis <- plot(curved_fit, main = "OS384 Cisplatin IC50", xlab = "Concentration (uM)", 
     ylab = "Fluorescent Intensity", digits = 3)


# changing the name of the 2nd column to predictions
names(IC50_curve_384_cis)[2] <- 'predictions'


# plotting the dos response curve with ggplot2
IC50_curve_384_cis_gg <- ggplot(IC50_curve_384_cis, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs(title = "OS384 Cisplatin IC50", x = "Concentration (uM)", y = "Fluorescent Intensity") +
  theme_bw()


# Save the plot as an SVG file
ggsave("~/Desktop/Cis_IC50_OS384.svg", plot = IC50_curve_384_cis_gg, device = "svg")

