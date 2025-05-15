

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


# Defining the plotting function
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
      labs(title = paste("OS186 Viability -", target_name),
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


# 
plot_viability_by_target(data = OS833_NFE2L3_viability, output_dir = "~/Desktop/")



######    OS186 NFE2L3 viability analysis   ######



# Reading in the data
OS186_NFE2L3_viability <- read.csv("~/Desktop/Reprogramming_Osteosarcoma/Viability_analysis/Processed_2025_03_19_OS186_NFE2L3.csv")

OS186_NFE2L3_viability$Drug[37:48] <- "CDK4/6"

names(OS186_NFE2L3_viability)[1] <- "target"

names(OS186_NFE2L3_viability)[2] <- "drug"

names(OS186_NFE2L3_viability)[3] <- "value"


# 
plot_viability_by_target(data = OS186_NFE2L3_viability, output_dir = "~/Desktop/")







