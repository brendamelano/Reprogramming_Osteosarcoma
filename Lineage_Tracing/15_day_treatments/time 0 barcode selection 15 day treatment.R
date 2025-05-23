library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(statmod)
library(ggpubr)
library(bioseq)
library(dplyr)
library(tidyr)




#############      OS052 TIME 0 BARCODES      ################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_1_S1_L001_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_2_S2_L001_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_0_-_3_S3_L001_R1_001.fastq.gz.out1.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS052_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Renaming the count columns
names(OS052_ctrl_0_merged)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Computing standard deviation of the barcodes
OS052_ctrl_0_merged <- OS052_ctrl_0_merged %>%
  rowwise() %>%
  mutate(StdDev = sd(c(barcode_count_ctrl_0_1, barcode_count_ctrl_0_2, barcode_count_ctrl_0_3), na.rm = TRUE)) %>%
  ungroup()


# Defining a threshold for the standard deviation cutoff
# Retry analysis with more stringent cutoff
threshold <- 1500


# Filtering based on the stdev threshhold
filtered_df <- OS052_ctrl_0_merged %>%
  dplyr::filter(StdDev <= threshold)


# Removing the stdev column so that it is not scaled
filtered_df <- filtered_df[,-c(5)]


# Scaling the data based on cpm
OS052_ctrl_0_scaled <- cpm_scaling(filtered_df)


# Computing the logs of the cpm values
OS052ctrl0_log_scaled <- compute_log2_scaled(OS052_ctrl_0_scaled)


# Computing the mean of the cpm scaled values
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))


# Filtering barcodes with low counts to keep likely barcodes
OS052ctrl0_filtered <- OS052ctrl0_log_scaled %>% dplyr::filter(barcode_cpm_mean_ctrl_0 > 90)


# Renaming the barcode column
names(OS052ctrl0_filtered)[1] <- 'barcode'


# Making the list of barcodes for OS384 time 0 for a whitelist
OS052_time_0_barcodes <- OS052ctrl0_filtered$barcode

## Prepping data for visualization

# Ordering based on the average cpm per barcode
OS052ctrl0_log_scaled <- OS052ctrl0_log_scaled[order(-OS052ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]


# Creating an index for the order
OS052ctrl0_log_scaled$Index <- seq_along(OS052ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


#write.csv(OS052ctrl0_log_scaled, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/processed_counts/OS052_T0_counts.csv", row.names = TRUE)



## Plotting replicates

# 
# #Perform regression analysis
# model <- lm(barcode_count_ctrl_0_3_scaled_log ~ barcode_count_ctrl_0_2_scaled_log, data = OS052ctrl0_log_scaled)
# 
# 
# # Extract r-squared and p-value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot for replicate correlation in D0 control
# first_two_replicates_052_D0 <- ggplot(OS052ctrl0_log_scaled, aes(barcode_count_ctrl_0_3_scaled_log, barcode_count_ctrl_0_2_scaled_log)) +
#   geom_point_rast() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 7),
#         plot.title = element_text(size = 7)) +
#   xlab("Replicate 1") +
#   ylab("Replicate 2") +
#   ggtitle("OS052 Count Correlation Day-0") +
#   annotate("text",
#            x = min(OS052ctrl0_log_scaled$barcode_count_ctrl_0_2_scaled_log),
#            y = max(OS052ctrl0_log_scaled$barcode_count_ctrl_0_3_scaled_log),
#            label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
#            hjust = 0, vjust = 1, size = 3)
# 
# 
# first_two_replicates_052_D0
# 
# 
# ggsave("~/Desktop/first_two_replicates_052_D0.svg",
#        plot = first_two_replicates_052_D0,
#        device = "svg",
#        width = 2.2,  # Width in inches
#        height = 2.2, # Height in inches
#        dpi = 300)  # DPI, optional for SVG




# Ranked barcode plot with rasterized points
# plot <- ggplot(OS052ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
#   geom_line() + # Draw lines
#   rasterise(geom_point(), dpi = 300) + # Add rasterized points
#   scale_y_log10() + # Log scale for Y axis
#   labs(title = "OS052 LT ranked barcode plot", y = "mean cpm", x = "Ranked LT barcodes") + # Add titles and labels
#   theme_bw() + # Use a minimal theme for a cleaner look
#   theme(panel.grid.major = element_blank(), # Remove major grid lines
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 10),
#         plot.title = element_text(size = 10)) # Remove minor grid lines
# 
# 
# # Save the plot as SVG with specified dimensions
# ggsave("~/Desktop/OS052_lineage_tracing_plot.svg", plot = plot, width = 2.5, height = 2.5, units = "in")



# ###   Filtering barcodes for trajectory visualization   ###
# 
# 
# # getting the 10% quantiles
# #quantile(OS384_ctrl_0_unique$barcode_count_ctrl_0, probs = seq(0, 1, 0.1))
# 
# 
# # fitering the barcodes for those with a high count for trajectory analysis
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
# #write.csv(OS384_trajectory_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")
# 
# 



###############     OS384 TIME 0 BARCODES       ##################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/10_384_ctrl0_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/11_384_ctrl0_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/12_384_ctrl0_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Renaming the count columns
names(OS384_ctrl_0_merged)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Performing cpm scaling
OS384_ctrl_0_scaled <- cpm_scaling(OS384_ctrl_0_merged)


# Computing the logs of the cpm values
OS384ctrl0_log_scaled <- compute_log2_scaled(OS384_ctrl_0_scaled)


# Computing the mean per barcode for the merged dataframe
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% 
  mutate(barcode_log_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled_log", 
                                                        "barcode_count_ctrl_0_2_scaled_log", 
                                                        "barcode_count_ctrl_0_3_scaled_log"))))

# Computing the mean of the cpm scaled values
OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))


# Ordering based on cpm means
#OS384ctrl0_log_scaled <- OS384ctrl0_log_scaled[order(-OS384ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]


OS384ctrl0_log_scaled$Index <- seq_along(OS384ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


# Use ggplot to create the Ranked barcode plot
# p <- ggplot(OS384ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
#   geom_line() + # Draw lines
#   rasterise(geom_point(), dpi = 300) + # Add rasterized points
#   scale_y_log10() + # Log scale for Y axis
#   labs(title = "OS384 LT ranked barcode plot", y = "mean cpm", x = "Ranked LT barcodes") + # Add titles and labels
#   theme_bw() + # Use a minimal theme for a cleaner look
#   theme(panel.grid.major = element_blank(), # Remove major grid lines
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 10),
#         plot.title = element_text(size = 10)) # Remove minor grid lines
# 
# p
# # 
# # 
# # # Save the plot as SVG with specified dimensions
# ggsave("~/Desktop/OS384_lineage_tracing_plot.svg", plot = p, width = 2.5, height = 2.5, units = "in")


# filter barcodes to only keep those that have counts above 2 (first identified the elbow) by plotting the counts in order
OS384ctrl0_filtered <- OS384ctrl0_log_scaled %>% filter(barcode_log_mean_ctrl_0 > 2)


#write.csv(OS384ctrl0_filtered, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/OS384_T0_processed_counts.csv", row.names = TRUE)


# Renaming the barcode column
names(OS384ctrl0_filtered)[1] <- "barcode"


# making the list of barcodes for OS384 time 0 for a whitelist
time_0_barcodes <- OS384ctrl0_filtered$barcode


### Performing regression for r^2 value of replicates ##


# Perform regression analysis
# model <- lm(barcode_count_ctrl_0_3_log ~ barcode_count_ctrl_0_2_log, data = OS384ctrl0_log_scaled)
# 
# 
# # Extract r-squared and p-value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot for replicate correlation in D0 control
# first_two_replicates_384_D0 <- ggplot(OS384ctrl0_log_scaled, aes(barcode_count_ctrl_0_3_log, barcode_count_ctrl_0_2_log)) +
#   geom_point_rast() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 7),  
#         plot.title = element_text(size = 7)) +  
#   xlab("Replicate 1") +
#   ylab("Replicate 2") +
#   ggtitle("OS384 Count Correlation Day-0") +
#   annotate("text",
#            x = min(OS384ctrl0_log_scaled$barcode_count_ctrl_0_1_log),
#            y = max(OS384ctrl0_log_scaled$barcode_count_ctrl_0_2_log),
#            label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
#            hjust = 0, vjust = 1, size = 3)


# Save the plot as an SVG file
#ggsave("~/Desktop/first_two_replicates_384_D0.png", plot = first_two_replicates_384_D0, device = "png", width = 2.2, height = 2.2)



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




#############      OS742 TIME 0 BARCODES      ################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/1_742_ctrl0_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/2__742_ctrl0_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/3__742_ctrl0_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS742_ctrl_0_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)



# Scaling the data based on cpm
OS742_ctrl_0_scaled <- cpm_scaling(OS742_ctrl_0_merged)


# Renaming the count columns
names(OS742_ctrl_0_scaled)[2:4] <- c("barcode_count_ctrl_0_1", "barcode_count_ctrl_0_2", "barcode_count_ctrl_0_3")


# Renaming the columns with scaled values
names(OS742_ctrl_0_scaled)[5:7] <- c("barcode_count_ctrl_0_1_scaled", "barcode_count_ctrl_0_2_scaled", "barcode_count_ctrl_0_3_scaled")


# Computing the logs of the cpm values
OS742ctrl0_log_scaled <- compute_log2_scaled(OS742_ctrl_0_scaled)


# Computing the mean log value per barcode for the merged dataframe
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% 
  mutate(barcode_log_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled_log", 
                                                        "barcode_count_ctrl_0_2_scaled_log", 
                                                        "barcode_count_ctrl_0_3_scaled_log"))))


# Computing the mean log value per barcode for the merged dataframe
OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled %>% 
  mutate(barcode_cpm_mean_ctrl_0 = rowMeans(select(., c("barcode_count_ctrl_0_1_scaled", 
                                                        "barcode_count_ctrl_0_2_scaled", 
                                                        "barcode_count_ctrl_0_3_scaled"))))



OS742ctrl0_log_scaled <- OS742ctrl0_log_scaled[order(-OS742ctrl0_log_scaled$barcode_cpm_mean_ctrl_0), ]


OS742ctrl0_log_scaled$Index <- seq_along(OS742ctrl0_log_scaled$barcode_cpm_mean_ctrl_0)


# # Use ggplot to create the plot
# p <- ggplot(OS742ctrl0_log_scaled, aes(x = Index, y = barcode_cpm_mean_ctrl_0)) +
#   geom_line() + # Draw lines
#   rasterise(geom_point(), dpi = 300) + # Add rasterized points
#   scale_y_log10() + # Log scale for Y axis
#   labs(title = "OS742 LT ranked barcode plot", y = "mean cpm", x = "Ranked LT barcodes") + # Add titles and labels
#   theme_bw() + # Use a minimal theme for a cleaner look
#   theme(panel.grid.major = element_blank(), # Remove major grid lines
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 10),
#         plot.title = element_text(size = 10)) # Remove minor grid lines
# 
# 
# # Save the plot as SVG with specified dimensions
# ggsave("~/Desktop/OS742_lineage_tracing_plot.svg", plot = p, width = 2.5, height = 2.5, units = "in")



# filter barcodes to only keep those that have counts above 2 (first identified the elbow) by plotting the 
# counts in order
OS742ctrl0_filtered <- OS742ctrl0_log_scaled %>% filter(barcode_log_mean_ctrl_0 > 2)


#write.csv(OS742ctrl0_filtered, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS742/OS742_T0_processed_counts.csv", row.names = TRUE)


# Making the list of barcodes for OS384 time 0 for a whitelist
OS742_time_0_barcodes <- OS742ctrl0_filtered$V1


## performing regression for r^2 value of replicates ##

# 
# # Perform regression analysis
# model <- lm(barcode_count_ctrl_0_1_scaled_log ~ barcode_count_ctrl_0_2_scaled_log, data = OS742ctrl0_filtered)
# 
# 
# # Extract r-squared and p-value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot for replicate correlation in D0 control
# first_two_replicates_742_D0 <- ggplot(OS742ctrl0_filtered, aes(barcode_count_ctrl_0_1_scaled_log, barcode_count_ctrl_0_2_scaled_log)) +
#   geom_point_rast() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 7),  
#         plot.title = element_text(size = 7)) +  
#   xlab("Replicate 1") +
#   ylab("Replicate 2") +
#   ggtitle("OS742 Count Correlation Day-0") +
#   annotate("text",
#            x = min(OS742ctrl0_filtered$barcode_count_ctrl_0_1_scaled_log),
#            y = max(OS742ctrl0_filtered$barcode_count_ctrl_0_2_scaled_log),
#            label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
#            hjust = 0, vjust = 1, size = 3)
# 
# first_two_replicates_742_D0
# 
# # Save the plot as an SVG file
# ggsave("~/Desktop/first_two_replicates_742_D0.svg", plot = first_two_replicates_742_D0, device = "svg", width = 2.2, height = 2.2, units = "in")
# 


