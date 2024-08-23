#library(VennDiagram)
library(stringdist)
#library(DescTools)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(stringr)
library(ggpubr)
library(stats)
library(dplyr)
library(tidyr)
library(purrr)
library(mgcv) # GLMGAM regression
library(grid)
library(png)



############      PROCESSING SAMPLES FOR CTRL D13      #####################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/13_384_ctrl13_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/14_384_ctrl13_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/15_384_ctrl13_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_ctrl_13_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# 
OS384_ctrl_13_merged <- process_and_filter_barcodes(OS384_ctrl_13_merged, "ctrl_13", time_0_barcodes)


# Performing cpm scaling with the function
OS384_ctrl_13_scaled <- cpm_scaling(OS384_ctrl_13_merged)

OS384_ctrl13_log_scaled <- log_scales_and_means(OS384_ctrl_13_scaled, "ctrl_13")


# # Perform regression analysis
# model <- lm(barcode_count_ctrl_13_2_log ~ barcode_count_ctrl_13_3_log, data = OS384ctrl13_log_scaled)
# 
# 
# # Extract r-squared value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot
# OS384_D13_replicate <- ggplot(OS384ctrl13_log_scaled, aes(barcode_count_ctrl_13_2_log, barcode_count_ctrl_13_3_log)) +
#   geom_point() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   xlab("Log Transformed Barcode Count - Replicate 1") +
#   ylab("Log Transformed Barcode Count - Replicate 2") +
#   ggtitle("OS384 Barcode Count Correlation") +
#   geom_text(x = min(OS384ctrl13_log_scaled$barcode_count_ctrl_13_3_log),
#             y = max(OS384ctrl13_log_scaled$barcode_count_ctrl_13_2_log),
#             label = paste("R-squared =", round(r_squared, 2), "\n"),
#             hjust = 0, vjust = 1, parse = TRUE)
# 
# 
# # Save the plot as an SVG file
# ggsave("~/Desktop/OS384_D13_replicate.svg", plot = OS384_D13_replicate, device = "svg")



###############     ATR BARCODE ANALYSIS       ##################



# Reading in the ATR gDNA barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/19_384_atr_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/20_384_atr_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/21_384_atr_3.fastq.gz.out2.txt')


# Applying the process file function to the file paths
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_atr_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Example usage:
OS384_atr_merged <- process_and_filter_barcodes(OS384_atr_merged, "atr", time_0_barcodes)


# Performing cpm scaling with the function
OS384_atr_scaled <- cpm_scaling(OS384_atr_merged)


# Merging the control and treated counts
OS384_atr_ctrl13 <-  merge(OS384_ctrl_13_scaled, OS384_atr_scaled, by='barcode')


# Creating a dataframe for the barcodes not in the white list
#test_sample_extra <- OS384_atr_ctrl13 %>% filter(!(barcode %in% time_0_barcodes))
#test_sample <- OS384_atr_ctrl13 %>% filter((barcode %in% time_0_barcodes))


## hamming distance section ##
# # Writing a for loop to add values from the extra values to the whitelist
# for (barcode in 1:nrow(test_sample_extra)) {
#   barcode_seq <- test_sample_extra$barcode[barcode]
#   for (barcode_white in 1:nrow(test_sample)){
#     barcode_white_seq <- test_sample$barcode[barcode_white]
#     if (StrDist(barcode_seq, barcode_white_seq, method = "hamming") == 1){
#       test_sample[barcode_white, 2:4] <- test_sample[barcode_white, 2:4] + test_sample_extra[barcode, 2:4]
#     }
#   }
# }


# Performing cpm scaling with the function
#OS384_atr_scaled <- cpm_scaling(OS384_atr_final)


# Log transforming the scaled values
OS384_atr_log_scaled <- compute_log2_scaled(OS384_atr_ctrl13)


# Computing the mean per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_atr_cpm = rowMeans(select(., c("barcode_count_atr_1_scaled", 
                                                 "barcode_count_atr_2_scaled", 
                                                 "barcode_count_atr_3_scaled"))))


# Computing the mean cpm per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_ctrl13_cpm = rowMeans(select(., c("barcode_count_ctrl_13_1_scaled", 
                                                        "barcode_count_ctrl_13_2_scaled", 
                                                        "barcode_count_ctrl_13_3_scaled"))))

# Computing the mean per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_atr_log = rowMeans(select(., c("barcode_count_atr_1_scaled_log", 
                                                     "barcode_count_atr_2_scaled_log", 
                                                     "barcode_count_atr_3_scaled_log"))))


# Computing the mean cpm per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_ctrl13_log = rowMeans(select(., c("barcode_count_ctrl_13_1_scaled_log", 
                                                        "barcode_count_ctrl_13_2_scaled_log", 
                                                        "barcode_count_ctrl_13_3_scaled_log"))))


OS384_atr_final <- OS384_atr_log_scaled


### PLOTTING THE REPLICATES
#Perform regression analysis
model <- lm(barcode_count_384_atr_1_scaled_log ~ barcode_count_384_atr_2_scaled_log, data = OS384_atr_final)


# Extract r-squared and p-value
r_squared <- summary(model)$r.squared


# Create the ggplot
OS384_atr_replicate <- ggplot(OS384_atr_final, aes(barcode_count_384_atr_2_scaled_log, barcode_count_384_atr_1_scaled_log)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 7),  
        plot.title = element_text(size = 7)) +  
  xlab("Replicate 1") +
  ylab("Replicate 2") +
  ggtitle("OS384 Count Correlation ATR-i") +
  annotate("text",
           x = min(OS384_atr_final$barcode_count_384_atr_2_scaled_log),
           y = max(OS384_atr_final$barcode_count_384_atr_1_scaled_log),
           label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
           hjust = 0, vjust = 1, size = 3)

OS384_atr_replicate


# Save the plot as an SVG file
ggsave("~/Desktop/OS384_atr_replicate.png", plot = OS384_atr_replicate,  device = "png", width = 2.2, height = 2.2)




###########    PF (PFIZER) BARCODES    ##################



# Reading in the PF barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/22_384_pf_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/23_384_pf_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/24_384_pf_3.fastq.gz.out2.txt')


# Combining all the text files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_pf_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


OS384_pf_merged <- process_and_filter_barcodes(OS384_pf_merged, "pf", time_0_barcodes)


# Performing cpm scaling with the function
OS384_pf_scaled <- cpm_scaling(OS384_pf_merged)


OS384_pf_log_scaled <- log_scales_and_means(OS384_pf_scaled, "pf")


# Merging the control and treated counts
OS384_pf_final <-  merge(OS384_ctrl13_log_scaled, OS384_pf_log_scaled, by='barcode')




###############     CIS BARCODES    ############



# Reading in the Cisplatin treated counts
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/16_384_cis_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/17_384_cis_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/18_384_cis_3.fastq.gz.out2.txt')


# Processing all the txt files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_cis_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


OS384_cis_merged <- process_and_filter_barcodes(OS384_cis_merged, "cis", time_0_barcodes)


# Performing cpm scaling with the function
OS384_cis_scaled <- cpm_scaling(OS384_cis_merged)


# 
OS384_cis_log_scaled <- compute_log2_scaled(OS384_cis_scaled)


# computing the mean per barcode for the merged dataframe
# Make sure that this does not have to be the mean
OS384_cis_log_scaled <- OS384_cis_log_scaled %>%
  mutate(barcode_mean_cis_log = rowMeans(select(., c("barcode_count_cis_1_scaled_log",
                                                     "barcode_count_cis_2_scaled_log",
                                                     "barcode_count_cis_3_scaled_log"))))


OS384_cis_log_scaled <- OS384_cis_log_scaled %>% 
  mutate(barcode_mean_cis_cpm = rowMeans(select(., c("barcode_count_cis_1_scaled", 
                                                     "barcode_count_cis_2_scaled", 
                                                     "barcode_count_cis_3_scaled"))))


# Merging the control and treated counts
OS384_cis_final <-  merge(OS384_ctrl13_log_scaled, OS384_cis_log_scaled, by='barcode')



########      QC for all samples    ##########



# Subset the dataframe to exclude columns with 'log', 'scaled', 'mean'
raw_counts_df_pf <- OS384_pf_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))


raw_counts_df_atr <- OS384_atr_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))

# Creating cis df to keep only raw values
raw_counts_df_cis <- OS384_cis_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))


# Merging the raw counts from the atr and pf treatment
raw_counts_df <- left_join(raw_counts_df_atr, raw_counts_df_pf)


# Assigning the na values to 0
raw_counts_df[is.na(raw_counts_df)] <- 0


# 
raw_counts_df <- raw_counts_df[,-1]


# Compute the sum of each column
sums <- colSums(raw_counts_df)


# Creating a df to keep the sums of all the samples
sums_df <- data.frame(
  sample_type = names(sums),
  total_sum = sums
)


# Create the plot
plot <- ggplot(sums_df, aes(x = sample_type, y = total_sum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "OS384 Count sums", x = "Sample", y = "Total Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate X labels for readability
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank()) # Remove minor gridlines


# Print the plot
print(plot)

ggsave("~/Desktop/OS384_total_counts.svg", plot, device = "svg", width = 3.3, height = 4)


# pltting mean and median
data_long <- pivot_longer(raw_counts_df, 
                          cols = starts_with("barcode_count"),
                          names_to = "condition",
                          values_to = "counts")


# Cap data at 1,000
data_long <- data_long %>%
  mutate(counts = ifelse(counts > 1000, 1000, counts))


# Calculate mean, median, and std
summary_data <- data_long %>%
  group_by(condition) %>%
  summarise(mean = mean(counts, na.rm = TRUE),
            median = median(counts, na.rm = TRUE),
            std = sd(counts, na.rm = TRUE)) %>%
  mutate(lower = mean - std, upper = mean + std)


# Merge the summary back to the original long data for plotting
data_long_summary <- merge(data_long, summary_data, by = "condition")
# Plotting
# Adjust x-axis labels
data_long_summary$condition <- data_long_summary$condition %>%
  str_replace("barcode_count_", "") %>% # Remove "barcode_count_"
  str_replace_all("_", " ") %>% # Replace underscores with spaces
  str_replace("ctrl 13", "Ctrl Day-13") # Change "ctrl 13" to "Ctrl Day-13"



# Plot
p <- ggplot(data_long_summary, aes(x = condition, y = counts)) +
  geom_jitter_rast(aes(color = "Data Points"), width = 0.2, height = 0, alpha = 0.5) + # Rasterized data points
  geom_errorbar(aes(ymin = lower, ymax = upper, x = condition), width = 0.2) + # Error bars
  geom_point(aes(y = mean, color = "Median"), size = 3) + # Median points
  scale_color_manual("", values = c("Median" = "red")) +
  theme_bw(base_size = 8) + # Set base font size to 8
  labs(title = "OS384 Count Distribution",
       x = "Condition",
       y = "Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # X-axis text size
        axis.text.y = element_text(size = 8), # Y-axis text size
        axis.title = element_text(size = 8), # Axis title size
        plot.title = element_text(size = 8), # Plot title size
        legend.position = "right", # Move legend to the right
        legend.text = element_text(size = 8), # Legend text size
        legend.title = element_text(size = 8)) # Legend title size

p

# Save as SVG with rasterized points
ggsave("~/Desktop/OS384_count_distribution.svg", plot = p, width = 3, height = 2.5)



