library(dplyr)
library(ggplot2)
library(tidyr)
library(VennDiagram)
library(ggpubr)
library(tidyverse)
library(statmod)
library(bioseq)



###############     384 TIME 0 BARCODES       ##################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/10_384_ctrl0_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/11_384_ctrl0_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/12_384_ctrl0_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Performing cpm scaling
#OS384_ctrl_0_scaled <- cpm_scaling(OS384_ctrl_0_merged)


# Compute the sum of each column
sums <- colSums(OS384_ctrl_0_merged[2:4])


# Dividing the sum by 1 million to get the scaling factor for cpm
scaling_factor <- sums / 1000000


# Renaming the merged dataframe before scaling it
df <- OS384_ctrl_0_merged


# Divide each value in the data frame by its corresponding scaling factor
for (i in 2:ncol(df)) {
  df[, i] <- df[, i] / scaling_factor[i - 1]
}


# Renaming the columns
names(df)[2:4] <- c("barcode_count_ctrl_0_1_scaled", "barcode_count_ctrl_0_2_scaled", "barcode_count_ctrl_0_3_scaled")


# Merging the counts with the scaled counts
OS384ctrl0_merged <- merge(OS384_ctrl_0_merged, df, by = "V1")

  
# Computing logs of cpm values
OS384ctrl0_log_scaled <- OS384ctrl0_merged %>% mutate(barcode_count_ctrl_0_1_log = log2(barcode_count_ctrl_0_1_scaled))
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_2_log = log2(barcode_count_ctrl_0_2_scaled))
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_3_log = log2(barcode_count_ctrl_0_3_scaled))


# Computing the mean per barcode for the merged dataframe
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% 
  mutate(barcode_log_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_log", 
                                                        "barcode_count_ctrl_0_2_log", 
                                                        "barcode_count_ctrl_0_3_log"))))


# Computing the mean of the cpm scaled values
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))

# Ordering based on cpm means
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled[order(-OS384ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]

OS384ctrl0_log_scaled$Index <- seq_along(OS384ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


# # Use ggplot to create the plot
# p <- ggplot(OS384ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
#   geom_line() + # Draw lines
#   geom_point() + # Add points
#   scale_y_log10() + # Log scale for Y axis
#   labs(title = "OS384 LT Ranked Barcodes", y = "Mean CPM", x = "Ranked LT Barcodes") + # Add titles and labels
#   theme_bw() + # Use a minimal theme for a cleaner look
#   theme(
#     panel.grid.major = element_blank(), # Remove major grid lines
#     panel.grid.minor = element_blank(), # Remove minor grid lines
#     plot.title = element_text(size = 9, face = "bold"), # Title font size
#     axis.title.x = element_text(size = 8), # X-axis title font size
#     axis.title.y = element_text(size = 8), # Y-axis title font size
#     axis.text.x = element_text(size = 8), # X-axis tick label font size
#     axis.text.y = element_text(size = 8) # Y-axis tick label font size
#   )
# 
# 
# # Saving svg file
# ggsave("~/Desktop/OS384_ranked_barcode_LT.svg", plot = p, device = "svg", width = 2.2, height = 2.2, units = "in")


# filter barcodes to only keep those that have counts above 2 (first identified the elbow) by plotting the 
# counts in order
OS384ctrl0_filtered <- OS384ctrl0_log_scaled %>% filter(barcode_log_mean_ctrl_0 > 2)


# Renaming the barcode column
names(OS384ctrl0_filtered)[1] <- "barcode"


# making the list of barcodes for OS384 time 0 for a whitelist
time_0_barcodes <- OS384ctrl0_filtered$barcode


# writing the csv to the single cell folder
#write.csv(time_0_barcodes, "~/Desktop/OS384_time0_barcodes.csv")


### Performing regression for r^2 value of replicates ##


# Perform regression analysis
model <- lm(barcode_count_ctrl_0_3_log ~ barcode_count_ctrl_0_2_log, data = OS384ctrl0_log_scaled)


# Extract r-squared and p-value
r_squared <- summary(model)$r.squared


# Create the ggplot for replicate correlation in D0 control
first_two_replicates_384_D0 <- ggplot(OS384ctrl0_filtered, aes(barcode_count_ctrl_0_3_log, barcode_count_ctrl_0_2_log)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Log Transformed Barcode Count - Replicate 1") +
  ylab("Log Transformed Barcode Count - Replicate 2") +
  ggtitle("OS384 Barcode Count Correlation Time 0") +
  geom_text(x = min(OS384ctrl0_log_scaled$barcode_count_ctrl_0_1_log),
            y = max(OS384ctrl0_log_scaled$barcode_count_ctrl_0_2_log),
            label = paste("R-squared =", round(r_squared, 2)),
            hjust = 0, vjust = 1, parse = TRUE)


# Save the plot as an SVG file
ggsave("~/Desktop/first_two_replicates_384_D0.svg", plot = first_two_replicates_384_D0, device = "svg")



###   Filtering barcodes for trajectory visualization   ###


# getting the 10% quantiles
#quantile(OS384_ctrl_0_unique$barcode_count_ctrl_0, probs = seq(0, 1, 0.1))


# fitering the barcodes for those with a high count for trajectory analysis
# OS384_ctrl_0_trajectory <- OS384ctrl0_filtered %>% filter(barcode_mean_ctrl_0 > 10)
# 
# 
# # converting the barcodes to dna object in order to get the reverse complement
# OS384_barcodes_trajectory <- dna(OS384_ctrl_0_trajectory$barcode)
# 
# 
# # getting the reverse complement of the top barcodes
# rc_384_trajectory_barcodes <- seq_complement(seq_reverse(OS384_barcodes_trajectory))
# 
# 
# # creating a dataframe of the reverse complement of the barcodes that should be used to study trajectories
# OS384_trajectory_barcodes <- as.data.frame(rc_384_trajectory_barcodes)
# 
# 
# # writing the csv to the single cell folder
# write.csv(OS384_trajectory_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")




#############      742 TIME 0 BARCODES      ################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/1_742_ctrl0_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/2__742_ctrl0_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/3__742_ctrl0_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS742_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Renaming the count columns
# Make sure that this renaming is not needed for the subsequent function
#names(OS742_ctrl_0_merged)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Scaling the data based on cpm
OS742_ctrl_0_scaled <- cpm_scaling(OS742_ctrl_0_merged)


# Renaming the count columns
names(OS742_ctrl_0_scaled)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Renaming the columns with scaled values
names(OS742_ctrl_0_scaled)[5:7] <- c("barcode_count_ctrl_0_1_scaled", "barcode_count_ctrl_0_2_scaled", "barcode_count_ctrl_0_3_scaled")


# Computing logs of cpm values
OS742ctrl0_log_scaled <- OS742_ctrl_0_scaled %>% mutate(barcode_count_ctrl_0_1_log = log2(barcode_count_ctrl_0_1_scaled))
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_2_log = log2(barcode_count_ctrl_0_2_scaled))
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_3_log = log2(barcode_count_ctrl_0_3_scaled))


# Computing the mean log value per barcode for the merged dataframe
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% 
  mutate(barcode_log_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_log", 
                                                        "barcode_count_ctrl_0_2_log", 
                                                        "barcode_count_ctrl_0_3_log"))))


# Computing the mean log value per barcode for the merged dataframe
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))



OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled[order(-OS742ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]

OS742ctrl0_log_scaled$Index <- seq_along(OS742ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


# Use ggplot to create the plot
p <- ggplot(OS742ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
  geom_line() + # Draw lines
  geom_point() + # Add points
  scale_y_log10() + # Log scale for Y axis
  #geom_vline(xintercept = 988, color = "red") + 
  labs(title = "OS742 Lineage Tracing ranked barcode plot", y = "Mean CPM", x = "Ranked LT Barcodes") + # Add titles and labels
  theme_bw() + # Use a minimal theme for a cleaner look
  theme(panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) # Remove minor grid lines

# Saving as svg
ggsave("~/Desktop/OS742_ranked_barcode_LT.svg", plot = p, device = "svg", width = 4, height = 4, units = "in")



# filter barcodes to only keep those that have counts above 2 (first identified the elbow) by plotting the 
# counts in order
OS742ctrl0_filtered <- OS742ctrl0_log_scaled %>% filter(barcode_log_mean_ctrl_0 > 2)


# Making the list of barcodes for OS384 time 0 for a whitelist
OS052_time_0_barcodes <- OS742ctrl0_filtered$barcode


## performing regression for r^2 value of replicates ##


# Perform regression analysis
model <- lm(barcode_count_ctrl_0_3_log ~ barcode_count_ctrl_0_2_log, data = OS742ctrl0_filtered)


# Extract r-squared and p-value
r_squared <- summary(model)$r.squared


# Create the ggplot for replicate correlation in D0 control
first_two_replicates_742_D0 <- ggplot(OS742ctrl0_filtered, aes(barcode_count_ctrl_0_3_log, barcode_count_ctrl_0_2_log)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Log Barcode Count - Replicate 1") +
  ylab("Log Barcode Count - Replicate 2") +
  ggtitle("OS742 Barcode Count Correlation") +
  geom_text(x = min(OS742ctrl0_filtered$barcode_count_ctrl_0_1_log),
            y = max(OS742ctrl0_filtered$barcode_count_ctrl_0_2_log),
            label = paste("R-squared =", round(r_squared, 2)),
            hjust = 0, vjust = 1, parse = TRUE)


# Save the plot as an SVG file
ggsave("~/Desktop/first_two_replicates_742_D0.svg", plot = first_two_replicates_742_D0, device = "svg", width = 4, height = 4, units = "in")




###   Filtering barcodes for trajectory visualization   ###


# getting the 10% quantiles
#quantile(OS384_ctrl_0_unique$barcode_count_ctrl_0, probs = seq(0, 1, 0.1))


# fitering the barcodes for those with a high count for trajectory analysis
OS384_ctrl_0_trajectory <- OS384ctrl0_filtered %>% filter(barcode_mean_ctrl_0 > 10)


# converting the barcodes to dna object in order to get the reverse complement
OS384_barcodes_trajectory <- dna(OS384_ctrl_0_trajectory$barcode)


# getting the reverse complement of the top barcodes
rc_384_trajectory_barcodes <- seq_complement(seq_reverse(OS384_barcodes_trajectory))


# creating a dataframe of the reverse complement of the barcodes that should be used to study trajectories
OS384_trajectory_barcodes <- as.data.frame(rc_384_trajectory_barcodes)


# writing the csv to the single cell folder
write.csv(OS384_trajectory_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")


seq_complement(seq_reverse(dna('CATGGCATGTATGAAAAC')))



#############      OS052 TIME 0 BARCODES      ################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_1_S28_L003_R1_001.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_2_S29_L003_R1_001.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_3_S30_L003_R1_001.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS052_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Scaling the data based on cpm
OS052_ctrl_0_scaled <- cpm_scaling(OS052_ctrl_0_merged)


# Renaming the count columns
names(OS052_ctrl_0_scaled)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Renaming the columns with scaled values
names(OS052_ctrl_0_scaled)[5:7] <- c("barcode_count_ctrl_0_1_scaled", "barcode_count_ctrl_0_2_scaled", "barcode_count_ctrl_0_3_scaled")


# Computing logs of cpm values
OS052ctrl0_log_scaled <- OS052_ctrl_0_scaled %>% mutate(barcode_count_ctrl_0_1_log = log2(barcode_count_ctrl_0_1_scaled))
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_2_log = log2(barcode_count_ctrl_0_2_scaled))
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>% mutate(barcode_count_ctrl_0_3_log = log2(barcode_count_ctrl_0_3_scaled))


# Computing the mean log value per barcode for the merged dataframe
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>% 
  mutate(barcode_log_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_log", 
                                                        "barcode_count_ctrl_0_2_log", 
                                                        "barcode_count_ctrl_0_3_log"))))


# Computing the mean of the cpm scaled values
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))


# Ordering based on the average cpm per barcode
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled[order(-OS052ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]


# Creating an index for the order
OS052ctrl0_log_scaled$Index <- seq_along(OS052ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


# # Use ggplot to create the ranked barcode plot
# ggplot(OS052ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
#   geom_line() + # Draw lines
#   geom_point() + # Add points
#   scale_y_log10() + # Log scale for Y axis
#   geom_vline(xintercept = 900, color = "red") + # Vertical line at x = 900
#   labs(title = "OS052 Lineage Tracing ranked barcode plot", y = "mean cpm", x = "") + # Add titles and labels
#   theme_bw() + # Use a minimal theme for a cleaner look
#   theme(panel.grid.major = element_blank(), # Remove major grid lines
#   panel.grid.minor = element_blank()) # Remove minor grid lines


# Computing standard deviation of the barcodes
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>%
  rowwise() %>%
  mutate(StdDev = sd(c(barcode_count_ctrl_0_1, barcode_count_ctrl_0_2, barcode_count_ctrl_0_3), na.rm = TRUE)) %>%
  ungroup()


# Defining a threshold for the standard deviation cutoff
threshold <- 8000


# Filtering based on the stdev threshhold
filtered_df <- OS052ctrl0_log_scaled %>%
  filter(StdDev <= threshold)


## Plotting replicates

# # Perform regression analysis
# model <- lm(barcode_count_ctrl_0_3_log ~ barcode_count_ctrl_0_2_log, data = filtered_df)
# 
# 
# # Extract r-squared and p-value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot for replicate correlation in D0 control
# first_two_replicates_052_D0 <- ggplot(filtered_df, aes(barcode_count_ctrl_0_3_log, barcode_count_ctrl_0_2_log)) +
#   geom_point() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   xlab("Log Transformed Barcode Count - Replicate 1") +
#   ylab("Log Transformed Barcode Count - Replicate 2") +
#   ggtitle("OS052 Barcode Count Correlation") +
#   geom_text(x = min(OS052ctrl0_log_scaled$barcode_count_ctrl_0_1_log),
#             y = max(OS052ctrl0_log_scaled$barcode_count_ctrl_0_2_log),
#             label = paste("R-squared =", round(r_squared, 2)),
#             hjust = 0, vjust = 1, parse = TRUE)
# 
# first_two_replicates_052_D0

# ggsave("~/Desktop/first_two_replicates_052_D0.svg", 
#        plot = first_two_replicates_052_D0, 
#        device = "svg", 
#        width = 4,  # Width in inches
#        height = 4, # Height in inches
#        dpi = 300)  # DPI, optional for SVG

# filter barcodes to only keep those that have counts above 2 (first identified the elbow) by plotting the counts in order
OS052ctrl0_filtered <- filtered_df %>% filter(barcode_cpm_mean_ctrl_0 > 2)


# Making the list of barcodes for OS384 time 0 for a whitelist
OS052_time_0_barcodes <- OS052ctrl0_filtered$barcode





###   Filtering barcodes for trajectory visualization   ###


# getting the 10% quantiles
#quantile(OS384_ctrl_0_unique$barcode_count_ctrl_0, probs = seq(0, 1, 0.1))


# fitering the barcodes for those with a high count for trajectory analysis
OS384_ctrl_0_trajectory <- OS384ctrl0_filtered %>% filter(barcode_mean_ctrl_0 > 10)


# converting the barcodes to dna object in order to get the reverse complement
OS384_barcodes_trajectory <- dna(OS384_ctrl_0_trajectory$barcode)


# getting the reverse complement of the top barcodes
rc_384_trajectory_barcodes <- seq_complement(seq_reverse(OS384_barcodes_trajectory))


# creating a dataframe of the reverse complement of the barcodes that should be used to study trajectories
OS384_trajectory_barcodes <- as.data.frame(rc_384_trajectory_barcodes)


# writing the csv to the single cell folder
#write.csv(OS384_trajectory_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")



