library(car)
library(drc)
library(ggplot2)
library(magrittr)
library(dplyr)


# reading in csv files
# for some reason only the 1st two rows are read in so I manually add the 3rd row below
OS742_pf <- read.csv("~/Desktop/IC50_data/2020_01_04_15day_holiday_384_742_cis_atr_pf/741_pf_15day - 2023_01_03_384_742_CIS_ATR_PF_15day_IC50.csv", nrow = 2, header = F)


# creating a row for the 3rd row vectors
third_row <- c(161107,1027943,1287906,2369980,1908851,1675216,1618844,1834831,2033248)


# binding the dataframe with the 3rd row vector
OS742_pf <- rbind(OS742_pf, third_row)


# renaming the row names
rownames(OS742_pf)[1] <- "Pf_1"
rownames(OS742_pf)[2] <- "Pf_2"
rownames(OS742_pf)[3] <- "Pf_3"


### pf
# computing the means of the of the fluorescent intensities for the different concentrations
OS742_pf <- rbind(OS742_pf, OS742_pf %>% summarise_if(is.numeric, mean))


# renaming the mean row
rownames(OS742_pf)[4] <- 'mean'


# making a vector with the concentrations
# pf concentration
concentrations <- c(2, 0.6666667, 0.2222222, 0.07407407, 0.02469136, 0.008230453, 0.002743484, 0.0009144947, 0)


# changing the column names to the concentrations 
names(OS742_pf) <- concentrations


# creating a dataframe with the cell counts for the mean counts
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS742_pf)[,'mean']
)


# removing a row
# not sure why this is here
# dose_response <- dose_response[-1,]


# plotting the dose response curve with the dose (x-axis) on a log scale
ggplot(
  data = dose_response, aes(x = conc, y = cell_count)) +
  geom_point() +
  geom_line() +
  labs(x = "Drug Concentration", y = "Cell Count") +
  ylim(0, NA) + scale_x_continuous(trans = "log10")


# 
ggplot(dose_response, aes(x = conc, y = cell_count)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE) +
  ylim(0, NA)



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
plot(curved_fit)


# Note the logged x-axis
# This also actually creates a dataframe that can be used to plot the same data in ggplot2
IC50_curve_742_ATR <- plot(curved_fit, main = "OS742 CDK-4/6 inhibitor IC50", xlab = "Concentration (uM)", 
                           ylab = "Fluorescent Intensity", digits = 3)


# changing the name of the 2nd column to predictions
names(IC50_curve_742_ATR)[2] <- 'predictions'


# plotting the dos response curve with ggplot2
IC50_curve_742_ATR <- ggplot(IC50_curve_742_ATR, aes(x = conc, y = predictions)) +
  geom_line() +  # Add a line plot
  scale_x_log10() +  # Set the x-axis to a logarithmic scale
  labs(title = "OS742 CDK-4/6 inhibitor IC50", x = "Concentration (uM)", y = "Fluorescent Intensity") +
  theme_bw()


# Save the plot as an SVG file
ggsave("~/Desktop/PF_IC50_OS742.svg", plot = IC50_curve_742_ATR, device = "svg")




