library(stats)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(stringdist)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(ggrastr)
library(stringr)
library(DESeq2)
library(tidyverse)
library(mgcv) # GLMGAM regression
library(purrr)
library(DescTools)
library(ggplot2)
library(grid)
library(png)



############      PROCESSING SAMPLES FOR CTRL D13      #####################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/4__742_ctrl13_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/5_742_ctrl13_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/6_742_ctrl13_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS742_ctrl_13_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)

# Example usage:
OS742_ctrl_13_merged <- process_and_filter_barcodes(OS742_ctrl_13_merged, "ctrl_13", OS742_time_0_barcodes)


# Performing cpm scaling with the function
OS742_ctrl_13_scaled <- cpm_scaling(OS742_ctrl_13_merged)


## PLOTTING THE REPLICATES


# # Perform regression analysis
# model <- lm(barcode_count_ctrl_13_2_log ~ barcode_count_ctrl_13_3_log, data = OS742_ctrl13_scaled)
# 
# 
# # Extract r-squared value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot
# OS742_D13_replicate <- ggplot(OS742_ctrl13_scaled, aes(barcode_count_ctrl_13_2_log, barcode_count_ctrl_13_3_log)) +
#   geom_point() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlab("Log Barcode Count - Replicate 1") +
#   ylab("Log Barcode Count - Replicate 2") +
#   ggtitle("OS742 Barcode Count Correlation") +
#   geom_text(x = min(OS742_ctrl13_scaled$barcode_count_ctrl_13_3_log),
#             y = max(OS742_ctrl13_scaled$barcode_count_ctrl_13_2_log),
#             label = paste("R-squared =", round(r_squared, 2), "\n"),
#             hjust = 0, vjust = 1, parse = TRUE)
# 
# ggsave("~/Desktop/OS742_ranked_barcode_LT.svg", plot = OS742_D13_replicate, device = "svg", width = 4, height = 4, units = "in")


#
# # Create a plot of just the data points (without axes, labels, and titles)
# plot_without_axes <- ggplot(OS384ctrl13_log_scaled, aes(barcode_count_ctrl_13_2_log, barcode_count_ctrl_13_3_log)) +
#   geom_point() +
#   theme_void()  # This removes all axes, text, etc.
# 
# 
# # Save the data points as a PNG
# png(filename = "data_points.png", width = 800, height = 600)
# print(plot_without_axes)
# dev.off()
# 
# # Read the PNG back in
# img <- rasterGrob(readPNG("data_points.png"), interpolate=TRUE)
# 
# # Overlay the axes, labels, and titles on the PNG
# OS384_D13_replicate <- ggplot(OS384ctrl13_log_scaled, aes(barcode_count_ctrl_13_2_log, barcode_count_ctrl_13_3_log)) +
#   annotation_custom(img, xmin = min(OS384ctrl13_log_scaled$barcode_count_ctrl_13_2_log), 
#                     xmax = max(OS384ctrl13_log_scaled$barcode_count_ctrl_13_2_log), 
#                     ymin = min(OS384ctrl13_log_scaled$barcode_count_ctrl_13_3_log), 
#                     ymax = max(OS384ctrl13_log_scaled$barcode_count_ctrl_13_3_log)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   xlab("Log Transformed Barcode Count - Replicate 1") +
#   ylab("Log Transformed Barcode Count - Replicate 2") +
#   ggtitle("OS384 Barcode Count Correlation") +
#   geom_text(x = min(OS384ctrl13_log_scaled$barcode_count_ctrl_13_3_log),
#             y = max(OS384ctrl13_log_scaled$barcode_count_ctrl_13_2_log),
#             label = paste("R-squared =", round(r_squared, 2)),
#             hjust = 0, vjust = 1, parse = TRUE)



# Save the plot as an SVG file
# ggsave("~/Desktop/OS384_D13_replicate.svg", plot = OS384_D13_replicate, device = "svg")



###########    PF (PFIZER) CDK 2/4/6 BARCODES    ##################



# Reading in the PF barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/7_742_pf_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/8_742_pf_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/9_742_pf_3.fastq.gz.out2.txt')


# Combining all of the barcode dataframes
result_list <- lapply(file_paths, process_file)


# # Merging the data frames by 'V1'
OS742_pf_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


#
OS742_pf_merged <- process_and_filter_barcodes(OS742_pf_merged, "pf", OS742_time_0_barcodes)


# Performing cpm scaling with the function
OS742_pf_scaled <- cpm_scaling(OS742_pf_merged)


# Merging the control and treated counts
OS742_pf_ctrl13 <-  merge(OS742_ctrl_13_scaled, OS742_pf_scaled, by='barcode')


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


# 
OS742_pf_ctrl13 <- compute_log2_scaled(OS742_pf_ctrl13)



# Computing the mean per barcode for the merged dataframe
OS742_pf_ctrl13 <- OS742_pf_ctrl13 %>% 
  mutate(barcode_mean_pf_cpm = rowMeans(select(., c("barcode_count_pf_1_scaled", 
                                                     "barcode_count_pf_2_scaled", 
                                                     "barcode_count_pf_3_scaled"))))


# Computing the mean cpm per barcode for the merged dataframe
OS742_pf_ctrl13 <- OS742_pf_ctrl13 %>% 
  mutate(barcode_mean_ctrl13_cpm = rowMeans(select(., c("barcode_count_ctrl_13_1_scaled", 
                                                        "barcode_count_ctrl_13_2_scaled", 
                                                        "barcode_count_ctrl_13_3_scaled"))))


OS742_pf_final <- OS742_pf_ctrl13


## PLOTTING THE REPLICATES

# 
# # Perform regression analysis
# model <- lm(barcode_count_pf_1_scaled_log ~ barcode_count_pf_2_scaled_log, data = OS742_pf_final)
# 
# 
# # Extract r-squared value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot
# OS742_pf_replicate <- ggplot(OS742_pf_final, aes(barcode_count_pf_1_scaled_log, barcode_count_pf_2_scaled_log)) +
#   geom_point_rast() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 7),  
#         plot.title = element_text(size = 7)) + 
#   xlab("Replicate 1") +
#   ylab("Replicate 2") +
#   ggtitle("OS742 Count Correlation CDK-4/6 i ") +
#   annotate("text",
#            x = min(OS742ctrl0_filtered$barcode_count_ctrl_0_1_scaled_log),
#            y = max(OS742ctrl0_filtered$barcode_count_ctrl_0_2_scaled_log),
#            label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
#            hjust = 0, vjust = 1, size = 3)
# 
# 
# OS742_pf_replicate
# 
# ggsave("~/Desktop/OS742_pf_correlation.svg", plot = OS742_pf_replicate, device = "svg", width = 2.2, height = 2.2, units = "in")
# 


## running the function

#test_sample <- OS384_pf


# running the function
#OS384_pf <- test_barcode_analysis(test_sample, time_0_barcodes)

## Finishing the run of the function

# 
# # taking the difference of the log transformed values
# # need to add the barcode_mean_ctrl_13 (add the log suffix too)
# pf_diff_merged <- OS384_pf_final %>% mutate(difference_pf_log2 = barcode_log_mean_ctrl_13 - barcode_mean_pf_log)
# 
# 
# # scaling the data    
# pf_diff_merged <- pf_diff_merged %>% mutate(log2_diff_zscore_pf = scale(difference_pf_log2))
# 
# 
# # merging atr_diff_ordered with ctrl day 0 barcode counts by barcode
# pf_diff_merged <- merge(pf_diff_merged, OS384ctrl0_log_scaled, by = "barcode")
# 
# 
# # ordering based on the difference of z scores
# pf_diff_ordered <- pf_diff_merged[order(pf_diff_merged$difference_pf_log2, decreasing = T), ]
# 
# 
# # Plotting z score v counts at time 0
# ggplot(pf_diff_ordered, aes(x = barcode_mean_ctrl_0 , y = log2_diff_zscore_pf )) + 
#   geom_point() +
#   scale_x_continuous(trans='log10') + 
#   geom_line(y=1.96) + 
#   geom_line(y = -1.96) + 
#   theme_minimal() + 
#   ylab('Z-score of log2 difference') + 
#   xlab('log transformed counts - Day 0') +
#   ggtitle('Barcode selection with CDK 4/6 inhibitor')




##########    IDENTIFYING THE OVERLAPPING BARCODE DROPOUTS   #############




# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(deplted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)


# Writing the dropout barcodes to the single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/depleted_barcodes_OS384_inVivo_LT.csv")



enriched_barcodes <- unique(c(cis_barcodes, pf_barcodes, atr_barcodes))
enriched_barcodes <- dna(enriched_barcodes)
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)
write.csv(rc_enriched_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/enriched_barcodes_OS384_inVivo_LT.csv")



# Reading in the fastq files for read 1 and read 2 to combine the reads from the scRNAseq data
read1_file <- read.delim("/Users/bmelano/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/384-in-vivo_S1_L001_R1_001.fastq")
read2_file <- read.delim("/Users/bmelano/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/384-in-vivo_S1_L001_R2_001.fastq")



#######   Visualizing overlapping barcodes for the different samples and checking significance using hypergeometric tets   ########


# # finding the barcodes that overlapped in selected barcodes for all samples
# venn.diagram(
#   x = list(cis_barcodes_all, pf_barcodes_all, atr_barcodes_all),
#   category.names = c("Cisplatin" , "CDK 4/6 i", "ATR i"),
#   filename = '~/Desktop/overlapping_barcodes.svg',
#   output=TRUE,
#   height = 2150,
#   width = 2150,
#   main = "OS384 Barcode Selection Overlap"
# )

venn.diagram(
  x = list(cis_barcodes_all, pf_barcodes_all, atr_barcodes_all),
  category.names = c("Cisplatin" , "CDK 4/6 i", "ATR i"),
  filename = '~/Desktop/overlapping_barcodes.png',
  output=TRUE,
  height = 2150,
  width = 2150,
  main = "OS384 Barcode Selection Overlap"
)



# carrying out hypergeometric test for overlapping barcodes of all drugs
1 - phyper(q= 15,m = 824,n = 88,k = 215, lower.tail = T, log.p = F)


## checking correlation of z-score differences for different drugs

#merge all data frames in list
df_list_merged <- merge(z_score_diff_filtered_cis, z_score_diff_filtered_atr, by = 'barcode')
df_list_merged <- merge(df_list_merged, z_score_diff_filtered_pf)


# filtering the merged dataframe to only keep the columns that have the difference 
difference_4_drugs <- df_list_merged[,c(1, grep("zscore", names(df_list_merged)))]


# creating a scatter plot for the different drugs to see if the same barcodes are enriched or depleted with different drugs
ggplot(difference_4_drugs, aes(x = log2_diff__zscore_cis, y = log2_diff__zscore_pf)) + 
  geom_point() 
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')

  
### creating a list of the barcodes for dropouts based on log2 difference


# pivoting the data
dropout_barcodes_pivot <- difference_4_drugs %>%
  pivot_longer(colnames(difference_4_drugs)[2:5]) %>%
  as.data.frame()


# plotting the barcodes that have the highest dropouts in all 4 samples
ggplot(dropout_barcodes_pivot, aes(x = value)) + 
  geom_histogram(bins = 200) +
  facet_wrap(~ name) +
  scale_x_continuous(trans='log10')


# creatina a list of the comparisons
my_comparisons <- list( c("diference_atr, diference_cis"), c("diference_atr, diference_dox"), c("diference_atr, diference_pf") )


#
dropout_barcodes_pivot_filtered <- dropout_barcodes_pivot %>% filter(name %in% c("log2_diff_zscore_atr", "log2_diff_zscore_cis"))


# desktop
ggplot(dropout_barcodes_pivot_filtered, aes(x = name, y = value)) + 
  geom_boxplot() +
  scale_y_continuous(trans='log10')  +
  stat_compare_means() +
  stat_compare_means() 


ggplot(dropout_barcodes_pivot, aes(x = name, y = value)) + 
  geom_boxplot() +
  scale_y_continuous(trans='log10')  +
  stat_compare_means() +
  stat_compare_means() 



# creating a dataframe without the barcodes
diff_only <- difference_4_drugs[,-1]


######  QC analysis for all treatments together  ###



# Subset the dataframe to exclude columns with 'log', 'scaled', 'mean' from the atr dataset
raw_counts_df <- OS742_pf_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))



# Setting the na values to 0
raw_counts_df[is.na(raw_counts_df)] <- 0


raw_counts_df <- raw_counts_df[,-1]


# Compute the sum of each column
sums <- colSums(raw_counts_df)


sums_df <- data.frame(
  sample_type = names(sums),
  total_sum = sums
)


# Modify the sample_type to remove "barcode_count_" and underscores
sums_df$sample_type <- gsub("barcode_count_", "", sums_df$sample_type)
sums_df$sample_type <- gsub("ctrl_13", "Ctrl-D13", sums_df$sample_type)  # Capitalize the "c" in "ctrl"
sums_df$sample_type <- gsub("_", " ", sums_df$sample_type)


# Replace 'pf' with 'CDK-4/6 i' and 'atr' with 'ATR i'
sums_df$sample_type <- gsub("pf", "CDK-4/6 i", sums_df$sample_type)

# Order the factor levels so that control samples appear first
sums_df$sample_type <- factor(sums_df$sample_type, levels = sort(unique(sums_df$sample_type), decreasing = TRUE))

# Create the plot with font sizes set to 8
# Create the plot with a log scale on the y-axis and font sizes set to 8
plot <- ggplot(sums_df, aes(x = sample_type, y = total_sum)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +  # Set base font size to 8
  labs(title = "OS742 Count sums", x = "Sample", y = "Total Counts (log scale)") +
  scale_y_log10() +  # Set y-axis to log scale
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate X labels for readability, font size 8
        axis.text.y = element_text(size = 10), # Y-axis text size 8
        axis.title.x = element_text(size = 11), # X-axis title size 8
        axis.title.y = element_text(size = 11), # Y-axis title size 8
        plot.title = element_text(size = 11),   # Title font size 8
        panel.grid.major = element_blank(),    # Remove major gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

# Print the plot
print(plot)


ggsave("~/Desktop/OS052_total_counts.svg", plot, device = "svg", width = 2.5, height = 3)



### Plotting mean and median


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
  str_replace("ctrl 13", "Ctrl Day-13") %>% # Change "ctrl 13" to "Ctrl Day-13"
  str_replace("pf", "CDK-4/6 i") %>% # Change "pf" to "CDK-4/6 i"
  str_replace("atr", "Atr i") %>% # Change "atr" to "Atr i"
  str_to_title()  # Capitalize the first letter of each word


# Plot
p <- ggplot(data_long_summary, aes(x = condition, y = counts)) +
  geom_jitter_rast(aes(color = "Data Points"), width = 0.2, height = 0, alpha = 0.5) + # Rasterized data points
  geom_errorbar(aes(ymin = lower, ymax = upper, x = condition), width = 0.2) + # Error bars
  geom_point(aes(y = median, color = "Median"), size = 1) + # Median points
  scale_color_manual("", values = c("Median" = "red")) +
  theme_bw(base_size = 8) + # Set base font size to 8
  labs(title = "OS742 Count Distribution",
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
ggsave("~/Desktop/OS742_count_distribution.svg", plot = p, width = 2.5, height = 2.5)






