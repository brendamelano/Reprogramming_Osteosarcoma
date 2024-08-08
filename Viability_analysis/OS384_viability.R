#comparisons of interest
#cis treatment with stat1 KD and ZIC2 KD
#pf treatment with stat1 KD

library(ggplot2)
library(dplyr)
library(dplyr)
library(ggplot2)
library(tidyverse)


### NFE2L3    ######



# Reading in the data
OS384_NFE2L3_viability <- read.csv("~/Desktop/Reprogramming_Osteosarcoma/Viability_analysis/Processed_2024_08_07_OS384_NFE2L3_NR0B1_viability.csv")


# Defining the colors for ctrl and NFE2L3 KD
pastel_colors <- c("ctrl" = "#FFB3BA", "NFE2L3" = "#A6E1FA")


# 
OS384_NFE2L3_viability <- OS384_NFE2L3_viability %>%
  filter(drug %in% c("ATR", "ctrl"), target %in% c("ctrl", "NFE2L3"))


# Computing the mean for each gene/drug pair
summary_df <- OS384_NFE2L3_viability %>%
  group_by(drug, target) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Assuming your data frame OS833_NFE2L3_viability has a column 'value' for the measurement
t_test_results <- t.test(value ~ target,
                         data = filter(OS384_NFE2L3_viability, drug == "ATR"),
                         var.equal = TRUE)  # assuming equal variances for simplicity

p_value <- t_test_results$p.value  # Extracting the p-value from the test results


OS384_NFE2L3 <- ggplot(summary_df, aes(x=drug, y=mean_value, fill=target)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), 
                width=.2,                    # Adjust the width of the error bars
                position=position_dodge(.9)) +
  labs(title="OS384 Atr efficacy with NFE2L3 KD", 
       y="Fluorescent Measurement", x="Drug") +
  
  theme_bw() +  # Use theme_bw as a base theme
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  
  scale_fill_manual(values = pastel_colors) +
  
  # Annotate with the p-value
  annotate("text", x = 1, y = max(summary_df$mean_value + summary_df$sd_value), label = sprintf("p = %.3f", p_value))

print(OS384_NFE2L3)

# Save the plot as SVG to the Desktop
ggsave(filename = "~/Desktop/OS384_NFE2L3.svg", plot = OS384_NFE2L3, width = 4, height = 3, units = "in")


# Calculate the ratios for each gene
ratio_data <- OS833_STAT1_viability %>%
  pivot_wider(names_from = drug, values_from = value) %>%  # Use drug for column names and value for data
  mutate(ratio = cis/ctrl) %>%  # Calculate the ratio
  select(gene, ratio)     


# Plot the ratios
ggplot(ratio_data, aes(x=gene, y=ratio, fill=gene)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Cis/Ctrl Ratio by Gene", 
       y="Cis/Ctrl Ratio", x="Gene") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = pastel_colors)



######    JUN   ######



# Reading in the data
OS384_JUN_viability <- read.csv("~/Desktop/Reprogramming_Osteosarcoma/Viability_analysis/Processed_OS384_JUN_viavbility.xlsx - Sheet1.csv")

#OS833_STAT1_viability <- OS833_STAT1_viability %>% filter(!(drug %in% c("atr", "pf")))

pastel_colors <- c("ctrl" = "#FFB3BA", "JUN" = "#A6E1FA")

# OS384_JUN_viability <- OS384_JUN_viability %>%
#   filter(drug %in% c("cis", "ctrl"), gene %in% c("ctrl", "stat1"))


# Computing the mean for each gene/drug pair
summary_df <- OS384_JUN_viability %>%
  group_by(drug, gene) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )


OS384_JUN <- ggplot(summary_df, aes(x=drug, y=mean_value, fill=gene)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), 
                width=.2,                    # Adjust the width of the error bars
                position=position_dodge(.9)) +
  labs(title="In vitro Cis effiacy with STAT1 KD", 
       y="Fluorescent Measurement", x="Gene") +
  theme_bw() +  # Use theme_bw as a base theme
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  scale_fill_manual(values = pastel_colors)  # Use pastel colors

print(OS384_JUN)




# Calculate the ratios for each gene
ratio_data <- OS384_JUN_viability %>%
  pivot_wider(names_from = drug, values_from = value) %>%  # Use drug for column names and value for data
  mutate(ratio = Cis/DMSO) %>%  # Calculate the ratio
  select(gene, ratio)     


# Plot the ratios
ggplot(ratio_data, aes(x=gene, y=ratio, fill=gene)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Cis/Ctrl Ratio by Gene", 
       y="Cis/Ctrl Ratio", x="Gene") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = pastel_colors)




#### scaling by dividing JUN cells by WT cells


# Assuming you have a control measurement
# Compute the scaling factor
scaling_factor <- mean(OS384_JUN_viability$value[OS384_JUN_viability$gene == "ctrl"]) / 
  mean(OS384_JUN_viability$value[OS384_JUN_viability$gene == "JUN"])


# Apply the scaling factor to the JUN cell measurements
OS384_JUN_viability$value_scaled <- ifelse(OS384_JUN_viability$gene == "JUN",
                                           OS384_JUN_viability$value * scaling_factor,
                                           OS384_JUN_viability$value)

# Compute the mean for each gene/drug pair with the scaled values
summary_df <- OS384_JUN_viability %>%
  group_by(drug, gene) %>%
  summarize(
    mean_value = mean(value_scaled, na.rm = TRUE),
    sd_value = sd(value_scaled, na.rm = TRUE)
  )

# Plotting
OS384_JUN <- ggplot(summary_df, aes(x=drug, y=mean_value, fill=gene)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), 
                width=.2,                    
                position=position_dodge(.9)) +
  labs(title="In vitro Cis effiacy with STAT1 KD", 
       y="Fluorescent Measurement", x="Gene") +
  theme_bw() +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = pastel_colors)

print(OS384_JUN)




