#comparisons of interest
#cis treatment with stat1 KD and ZIC2 KD
#pf treatment with stat1 KD

library(ggplot2)
library(dplyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(broom)

######    NFE2L3   ######



# Reading in the data
OS833_NFE2L3_viability <- read.csv("~/Desktop/Reprogramming_Osteosarcoma/Viability_analysis/Processed_2024_04_23_833_NFE2L3_viability_15_day.csv")
OS833_NFE2L3_viability$drug[25:36] <- "CDK4/6"

# Define the plotting function
plot_viability_by_target <- function(data, targets = NULL, output_dir = "~/Desktop/") {
  
  # If targets are NULL, include all targets excluding control
  if (is.null(targets)) {
    targets <- setdiff(unique(data$target), "ctrl")
  }
  
  # Define colors for targets (ensure all your targets are included)
  pastel_colors <- c("ctrl" = "#FFB3BA", 
                     "nfe2l3" = "#A6E1FA",
                     "NFE2L3" = "#A6E1FA",
                     "NR0B1" = "#A1C299")
  
  # Loop over each target
  for (target_name in targets) {
    # Filter data for current target and control
    filtered_data <- data %>% filter(target %in% c("ctrl", target_name))
    
    # List of drugs that have both ctrl and target_name
    drugs_with_both <- filtered_data %>%
      group_by(drug) %>%
      filter(n_distinct(target) == 2) %>%
      distinct(drug) %>%
      pull(drug)
    
    # Further filter data to include only drugs that have both groups
    filtered_data <- filtered_data %>% filter(drug %in% drugs_with_both)
    
    # Check if there are any drugs left to plot
    if (nrow(filtered_data) == 0) {
      warning(paste("No data available for target:", target_name))
      next
    }
    
    # Compute summary statistics
    summary_df <- filtered_data %>%
      group_by(drug, target) %>%
      summarize(
        mean_value = mean(value, na.rm = TRUE),
        sd_value = sd(value, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Compute p-values for each drug
    p_values <- filtered_data %>%
      group_by(drug) %>%
      summarize(
        p_value = {
          # Perform t-test
          t_test_result <- tryCatch({
            t.test(value ~ target, data = pick(everything()), var.equal = TRUE)
          }, error = function(e) {
            warning(paste("t-test failed for drug:", unique(pick(everything())$drug), "-", e$message))
            return(NA_real_)
          })
          t_test_result$p.value
        },
        .groups = 'drop'
      )
    
    # Add p-values to summary_df
    summary_df <- left_join(summary_df, p_values, by = "drug")
    
    # Prepare data for p-value labels
    p_label_df <- summary_df %>%
      group_by(drug) %>%
      summarize(
        p_value = ifelse(length(na.omit(p_value)) > 0, na.omit(p_value)[1], NA_real_),
        y_position = max(mean_value + sd_value, na.rm = TRUE),
        .groups = 'drop'
      )
    
    
    # Create the plot
    p <- ggplot(summary_df, aes(x = drug, y = mean_value, fill = target)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                    width = 0.2,
                    position = position_dodge(0.9)) +
      labs(title = paste("OS833 Viability -", target_name),
           y = "Fluorescent Measurement", x = "Drug") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12)  # Change the font size here
      ) +
      scale_fill_manual(values = pastel_colors) +
      geom_text(data = p_label_df,
                aes(x = drug,
                    y = y_position + 0.05 * max(summary_df$mean_value + summary_df$sd_value, na.rm = TRUE),
                    label = ifelse(!is.na(p_value),
                                   sprintf("p = %.3f", p_value),
                                   "N/A")),
                inherit.aes = FALSE,
                vjust = 0)
    
    # Print the plot
    print(p)
    
    # Save the plot
    output_filename <- paste0(output_dir, "Viability_Plot_", target_name, ".svg")
    ggsave(filename = output_filename, plot = p, width = 4.5, height = 2.5, units = "in")
  }
}


names(OS833_NFE2L3_viability)[2] <- "target"
OS833_NFE2L3_viability$target[OS833_NFE2L3_viability$target == 'nfe2l3'] <- 'NFE2L3'


# Assuming your data frame is named OS384_NFE2L3_viability
plot_viability_by_target(data = OS833_NFE2L3_viability, output_dir = "~/Desktop/")



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
