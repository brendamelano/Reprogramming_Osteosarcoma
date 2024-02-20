library(car)
library(drc)
library(ggplot2)
library(magrittr)
library(dplyr)


# reading in csv files
# removed the first column
OS384_pf <- read.csv("~/Desktop/IC50_data/2020_01_04_15day_holiday_384_742_cis_atr_pf/384_15day_pf - 2023_01_03_384_742_CIS_ATR_PF_15day_IC50.csv", header = F)


# renaming the row names
rownames(OS384_pf)[1] <- "Pf_1"
rownames(OS384_pf)[2] <- "Pf_2"
rownames(OS384_pf)[3] <- "Pf_3"


# computing the mean of the columns
OS384_pf <- rbind(OS384_pf, OS384_pf %>% summarise_if(is.numeric, mean))
rownames(OS384_pf)[4] <- 'mean'


# making a vector with the concentrations
# Pf concentrations
concentrations <- c(2, 0.6666667, 0.2222222, 0.07407407, 0.02469136, 0.008230453, 0.002743484,0.0009144947, 0)


# changing the column names to the concentrations 
names(OS384_pf) <- concentrations


# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS384_pf)[,'mean']
)


# fiting a linear model to the data
fit <- lm(cell_count ~ conc, data = dose_response)


summary(fit)


# making predictions for the different concentrations
predict(fit)


# fitting an IC50 curve to the data
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


ic_50



# Note the logged x-axis
IC50_curve_384_CDK4_6 <- plot(curved_fit, main = "OS384 CDK 4/6 inhibitor IC50", xlab = "Concentration (uM)", 
     ylab = "Fluorescent Intensity", digits = 3)


# changing the name of the 2nd column to predictions
names(IC50_curve_384_CDK4_6)[2] <- 'predictions'


# plotting the dose response curve with ggplot2
IC50_curve_384_CDK4_6_gg <- ggplot(IC50_curve_384_CDK4_6, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs( x = "CDK 4/6 Concentration (uM)", y = "OS384 Fluorescent Intensity") +
  theme_bw()


# Save the plot as an SVG file
ggsave("~/Desktop/PF_IC50_OS384.svg", plot = IC50_curve_384_CDK4_6_gg, device = "svg")
