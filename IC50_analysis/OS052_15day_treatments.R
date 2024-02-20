library(car)
library(drc)
library(ggplot2)
library(magrittr)
library(dplyr)


######    PF    ########

# reading in csv files
# removed the first column
OS052_pf <- read.csv("~/Desktop/IC50_data/2023_01_31_OS152_052_833_cis_atr_pf_15day.csv", header = F)
OS052_pf <- OS052_pf[3:5,3:11]


# renaming the row names
rownames(OS052_pf)[1] <- "Pf_1"
rownames(OS052_pf)[2] <- "Pf_2"
rownames(OS052_pf)[3] <- "Pf_3"



# computing the mean for the replicates
OS052_pf <- rbind(OS052_pf, OS052_pf %>% summarise_if(is.numeric, mean))


rownames(OS052_pf)[4] <- 'mean'



# making a vector with the concentrations
# pf concentrations
concentrations <- c(2, 0.6666667, 0.2222222, 0.07407407, 0.02469136, 0.008230453, 0.002743484,0.0009144947, 0)



# changing the column names to the concentrations 
names(OS052_pf) <- concentrations




# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS052_pf)[,'mean']
)


dose_response <- dose_response[-1,]


# plotting the dose response curve with the dose (x-axis) on a log scale
ggplot(
  data = dose_response, aes(x = conc, y = cell_count)) +
  geom_point() +
  geom_line() +
  labs(x = "Drug Concentration", y = "Cell Count") +
  ylim(0, NA) + scale_x_continuous(trans = "log10")


# fiting a linear model to the data
fit <- lm(cell_count ~ conc, data = dose_response)


summary(fit)


predict(fit)



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

ic_50


# Note the logged x-axis
plot(curved_fit)


# Save the plot as an SVG file
ggsave("~/Desktop/PF_IC50_052.svg", plot = curved_fit, device = "svg")




######    Cis    ########

# reading in csv files
# removed the first column
OS052_cis <- read.csv("~/Desktop/IC50_data/2023_01_31_OS152_052_833_cis_atr_pf_15day.csv", header = F)
OS052_cis <- OS052_cis[15:17,3:11]


# renaming the row names
rownames(OS052_cis)[1] <- "cis_1"
rownames(OS052_cis)[2] <- "cis_2"
rownames(OS052_cis)[3] <- "cis_3"



# computing the mean for the replicates
OS052_cis <- rbind(OS052_cis, OS052_cis %>% summarise_if(is.numeric, mean))


rownames(OS052_cis)[4] <- 'mean'



# making a vector with the concentrations
# pf concentrations
concentrations <- c(5, 1.67, 0.56, 0.185, 0.0617284, 0.02057613, 0.00685871, 
                    0.002286237, 0)



# changing the column names to the concentrations 
names(OS052_cis) <- concentrations




# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS052_cis)[,'mean']
)


dose_response <- dose_response[-1,]



# fiting a linear model to the data
fit <- lm(cell_count ~ conc, data = dose_response)


summary(fit)


predict(fit)



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

ic_50


# Note the logged x-axis
plot(curved_fit)




######    Atr    ########

# reading in csv files
# removed the first column
OS052_atr <- read.csv("~/Desktop/IC50_data/2023_01_31_OS152_052_833_cis_atr_pf_15day.csv", header = F)
OS052_atr <- OS052_atr[18:20,3:11]


# renaming the row names
rownames(OS052_atr)[1] <- "atr_1"
rownames(OS052_atr)[2] <- "atr_2"
rownames(OS052_atr)[3] <- "atr_3"



# computing the mean for the replicates
OS052_atr <- rbind(OS052_atr, OS052_atr %>% summarise_if(is.numeric, mean))


rownames(OS052_atr)[4] <- 'mean'



# making a vector with the concentrations
# pf concentrations
concentrations <- c(.5, 0.167, 0.056, 0.0185, 0.00617284, 0.002057613, 0.000685871, 
                    0.0002286237, 0)



# changing the column names to the concentrations 
names(OS052_atr) <- concentrations




# creating a dataframe with the cell counts for one observation - need to take average instead of just using one observation
dose_response <- data.frame(
  conc = concentrations,
  cell_count = t(OS052_atr)[,'mean']
)


dose_response <- dose_response[-1,]



# fiting a linear model to the data
fit <- lm(cell_count ~ conc, data = dose_response)


summary(fit)

predict(fit)


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

ic_50


# Note the logged x-axis
plot(curved_fit)



