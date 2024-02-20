#comparisons of interest
#cis treatment with stat1 KD and ZIC2 KD
#pf treatment with stat1 KD

library(ggplot2)
library(dplyr)
library(dplyr)
library(ggplot2)
library(tidyverse)

######    STAT1   ######



# Reading in the data
OS833_STAT1_viability <- read.csv("~/Desktop/in_vitro_viability_analysis/Processed_833_STAT1_viability.csv")

#OS833_STAT1_viability <- OS833_STAT1_viability %>% filter(!(drug %in% c("atr", "pf")))

pastel_colors <- c("ctrl" = "#FFB3BA", "stat1" = "#A6E1FA")

OS833_STAT1_viability <- OS833_STAT1_viability %>%
  filter(drug %in% c("cis", "ctrl"), gene %in% c("ctrl", "stat1"))

# Computing the mean for each gene/drug pair
summary_df <- OS833_STAT1_viability %>%
  group_by(drug, gene) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# 
OS833_STAT1 <- ggplot(summary_df, aes(x=drug, y=mean_value, fill=gene)) +
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

print(OS833_STAT1)




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



######    ZIC2    #######
OS833_ZIC2_viability <- read.csv("~/Desktop/in_vitro_viability_analysis/Processed_833_ZIC2_viability.csv")

pastel_colors <- c("ctrl" = "#FFB3BA", "ZIC2" = "#A6E1FA")


OS833_ZIC2_viability <- OS833_ZIC2_viability %>%
  filter(drug %in% c("cis", "ctrl"), gene %in% c("ctrl", "ZIC2"))

summary_df_ZIC2 <- OS833_ZIC2_viability %>%
  group_by(drug, gene) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )


# 
# 
OS833_ZIC2 <- ggplot(summary_df_ZIC2, aes(x=drug, y=mean_value, fill=gene)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Fluorescent Measurement by Gene and Drug", 
       y="Fluorescent Measurement", x="Gene") +
  geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), 
                width=.2,                    # Adjust the width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +  # Use theme_bw as a base theme
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  scale_fill_manual(values = pastel_colors)  # Use pastel colors


print(OS833_ZIC2)
