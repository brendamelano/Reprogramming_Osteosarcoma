library(VennDiagram)
library(stringdist)
library(tidyverse)
library(ggplot2)
library(tidyverse)
#library(DESeq2)
#library(DescTools)
library(stats)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(stringr)
library(tidyr)
library(ggpubr)
library(mgcv) # GLMGAM regression
library(purrr)
library(grid)
library(png)



############      PROCESSING SAMPLES FOR CTRL D13      #####################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_13_-_1_S4_L001_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_13_-_2_S5_L001_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_ctrl_13_-_3_S6_L001_R1_001.fastq.gz.out1.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS052_ctrl_13_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Filtering barcodes based on the time 0 barcodes
OS052_ctrl_13_merged <- process_and_filter_barcodes(OS052_ctrl_13_merged, "ctrl_13", OS052_time_0_barcodes)


# Performing cpm scaling with the function
OS052_ctrl_13_scaled <- cpm_scaling(OS052_ctrl_13_merged)


OS052_ctrl_13_scaled <- log_scales_and_means(OS052_ctrl_13_scaled, "ctrl_13")


############   ATR ANALYSIS     ##############



# Reading in the ATR gDNA barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_atr_1_S154_L008_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_atr_2_S155_L008_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_atr_3_S156_L008_R1_001.fastq.gz.out1.txt')


# Applying the process file function to the file paths
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS052_atr_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Filtering based on time 0 barcodes
OS052_atr_merged <- process_and_filter_barcodes(OS052_atr_merged, "atr", OS052_time_0_barcodes)


# Performing cpm scaling with the function
OS052_atr_scaled <- cpm_scaling(OS052_atr_merged)


# Creating a dataframe for the barcodes not in the white list
#test_sample_extra <- OS052_atr_ctrl13 %>% dplyr::filter(!(barcode %in% OS052_time_0_barcodes))
#test_sample <- OS052_atr_ctrl13 %>% dplyr::filter((barcode %in% OS052_time_0_barcodes))


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


# Merging the control and treated counts
OS052_atr_ctrl13 <-  merge(OS052_ctrl_13_scaled, OS052_atr_scaled, by='barcode')



# Computing logs of cpm values
OS052ctrl13_log_scaled <- OS052_atr_ctrl13 %>% mutate(barcode_count_ctrl_13_1_log = log2(barcode_count_ctrl_13_1_scaled))
OS052ctrl13_log_scaled <- OS052ctrl13_log_scaled %>% mutate(barcode_count_ctrl_13_2_log = log2(barcode_count_ctrl_13_2_scaled))
OS052ctrl13_log_scaled <- OS052ctrl13_log_scaled %>% mutate(barcode_count_ctrl_13_3_log = log2(barcode_count_ctrl_13_3_scaled))



# Computing logs of cpm values
OS052_atr_log_scaled <- OS052ctrl13_log_scaled %>% mutate(barcode_count_atr_1_log = log2(barcode_count_atr_1_scaled))
OS052_atr_log_scaled <- OS052_atr_log_scaled %>% mutate(barcode_count_atr_2_log = log2(barcode_count_atr_2_scaled))
OS052_atr_log_scaled <- OS052_atr_log_scaled %>% mutate(barcode_count_atr_3_log = log2(barcode_count_atr_3_scaled))


# Computing the mean control log per barcode for the merged dataframe
OS052ctrl13_log_scaled <- OS052ctrl13_log_scaled %>%
   mutate(barcode_log_mean_ctrl_13 = rowMeans(select(., c("barcode_count_ctrl_13_1_log",
                                                          "barcode_count_ctrl_13_2_log",
                                                          "barcode_count_ctrl_13_3_log"))))




## PLOTTING THE REPLICATES






# Computing the mean per barcode for the merged dataframe
OS052_atr_log_scaled <- OS052_atr_log_scaled %>% 
  mutate(barcode_log_mean_atr = rowMeans(select(., c("barcode_count_atr_1_log", 
                                                 "barcode_count_atr_2_log", 
                                                 "barcode_count_atr_3_log"))))


# Computing the mean per barcode for the merged dataframe
OS052_atr_log_scaled <- OS052_atr_log_scaled %>% 
  mutate(barcode_cpm_mean_atr = rowMeans(select(., c("barcode_count_atr_1_scaled", 
                                                     "barcode_count_atr_2_scaled", 
                                                     "barcode_count_atr_3_scaled"))))


# Computing the mean per barcode for the merged dataframe
OS052_atr_log_scaled <- OS052_atr_log_scaled %>% 
  mutate(barcode_mean_ctrl13_cpm = rowMeans(select(., c("barcode_count_ctrl_13_1_scaled", 
                                                     "barcode_count_ctrl_13_2_scaled", 
                                                     "barcode_count_ctrl_13_3_scaled"))))


# Reassigning the test sample
OS052_atr_final <- OS052_atr_log_scaled



### PLOTTING THE REPLICATES TO CHECK PRECISION

# Perform regression analysis
# model <- lm(barcode_count_atr_1_log ~ barcode_count_atr_3_log, data = OS052_atr_final)
# 
# 
# # Extract r-squared and p-value
# r_squared <- summary(model)$r.squared
# 
# 
# # Create the ggplot
# OS052_atr_replicate <- ggplot(OS052_atr_final, aes(barcode_count_atr_3_log, barcode_count_atr_1_log)) +
#   geom_point_rast() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 7),
#         plot.title = element_text(size = 7)) +
#   xlab("Replicate 1") +
#   ylab("Replicate 2") +
#   ggtitle("OS052 Count Correlation ATR") +
#   annotate("text",
#            x = min(OS052ctrl0_log_scaled$barcode_count_ctrl_0_2_scaled_log),
#            y = max(OS052ctrl0_log_scaled$barcode_count_ctrl_0_3_scaled_log),
#            label = as.expression(bquote(R^2 == .(round(r_squared, 2)))),
#            hjust = 0, vjust = 1, size = 3)
# 
# 
# 
# 
# OS052_atr_replicate
# 
# 
# ggsave("~/Desktop/OS052_atr_replicate.svg",
#        plot = OS052_atr_replicate,
#        device = "svg",
#        width = 2.2,  # Width in inches
#        height = 2.2, # Height in inches
#        dpi = 300)  # DPI, optional for SVG




###########    PF (PFIZER) BARCODES    ##################



# Reading in the PF barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_pf_1_S157_L008_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_pf_2_S158_L008_R1_001.fastq.gz.out1.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/text_files/052_pf_3_S159_L008_R1_001.fastq.gz.out1.txt')


# Combining the counts into a single file
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS052_pf_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Filtering barcodes based on whitelist
OS052_pf_merged <- process_and_filter_barcodes(OS052_pf_merged, "pf", OS052_time_0_barcodes)


# Scaling the data based on cpm
OS052_pf_scaled <- cpm_scaling(OS052_pf_merged)


# Merging the control and treated counts
OS052_pf_ctrl13 <-  merge(OS052_ctrl_13_scaled, OS052_pf_scaled, by='barcode')


# Computing logs and means of barcode counts
OS052_pf_log_scaled <- log_scales_and_means(OS052_pf_ctrl13, "pf")


# Reassigning the test_sample
OS052_pf_final <- OS052_pf_log_scaled



######  QC analysis for all treatments together  ###



# Subset the dataframe to exclude columns with 'log', 'scaled', 'mean'
raw_counts_df_pf <- OS052_pf_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))


# Subset the dataframe to exclude columns with 'log', 'scaled', 'mean' from the atr dataset
raw_counts_df_atr <- OS052_atr_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))


# Merging the raw counts from the atr and pf treatment
raw_counts_df <- left_join(raw_counts_df_atr, raw_counts_df_pf)

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
sums_df$sample_type <- gsub("atr", "ATR i", sums_df$sample_type)

# Order the factor levels so that control samples appear first
sums_df$sample_type <- factor(sums_df$sample_type, levels = sort(unique(sums_df$sample_type), decreasing = TRUE))

# Create the plot with font sizes set to 8
# Create the plot with a log scale on the y-axis and font sizes set to 8
plot <- ggplot(sums_df, aes(x = sample_type, y = total_sum)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8) +  # Set base font size to 8
  labs(title = "OS052 Count sums", x = "Sample", y = "Total Counts (log scale)") +
  scale_y_log10() +  # Set y-axis to log scale
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # Rotate X labels for readability, font size 8
        axis.text.y = element_text(size = 8), # Y-axis text size 8
        axis.title.x = element_text(size = 8), # X-axis title size 8
        axis.title.y = element_text(size = 8), # Y-axis title size 8
        plot.title = element_text(size = 8),   # Title font size 8
        panel.grid.major = element_blank(),    # Remove major gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

# Print the plot
print(plot)


ggsave("~/Desktop/OS052_total_counts.svg", plot, device = "svg", width = 3, height = 3)



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
  labs(title = "OS052 Count Distribution",
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
ggsave("~/Desktop/OS052_count_distribution.svg", plot = p, width = 3, height = 2.5)




# #######   Visualizing overlapping barcodes for the different samples and checking significance using hypergeometric tests   ########
# 

library(VennDiagram)


depleted_pf <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/depleted_LT_barcodes_pf_OS052_LT.txt", header = F)
depleted_pf <- depleted_pf[,1]


depleted_atr <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/depleted_LT_barcodes_atr_OS052_LT.txt", header = F)
depleted_atr <- depleted_atr[,1]


# Specify the path to save the file on your Desktop
output_file <- "~/Desktop/barcode_overlap_venn.pdf"  # Modify path if necessary

# Set PDF graphic device to save the plot correctly
pdf(file = output_file, width = 7, height = 7)

# Generate the Venn diagram without saving to file
venn.plot <- venn.diagram(
  x = list("CDK-4/6 i depleted" = depleted_pf, "ATR i depleted" = depleted_atr),
  filename = NULL,  # Do not save to file directly
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  main = "OS052 depleted barcode overlap",
  main.cex = 2
)

library(grid)


# Draw the Venn diagram
grid.draw(venn.plot)

# Close the graphics device
dev.off()


# carrying out hypergeometric test for overlapping barcodes of all drugs
1 - phyper(q= 15,m = 824,n = 88,k = 215, lower.tail = T, log.p = F)
# 
OS052_pf_final_length <- as.double(nrow(OS052_pf_final))
OS052_atr_final_length <- as.double(nrow(OS052_atr_final))
N <- OS052_pf_final_length + OS052_atr_final_length

# 2. Calculate overlaps and non-overlaps
x <- as.double(length(intersect(depleted_pf, depleted_atr)))
b <- as.double(length(setdiff(depleted_pf, depleted_atr)))    # Depleted only in PF
c <- as.double(length(setdiff(depleted_atr, depleted_pf)))    # Depleted only in ATR
d <- N - x - b - c                                # Not depleted in either

# Create the contingency table
contingency_table <- matrix(c(x, b, c, d), nrow = 2,
                            dimnames = list(
                              "Depleted in ATR" = c("Yes", "No"),
                              "Depleted in PF" = c("Yes", "No")
                            ))
print("Contingency Table:")
print(contingency_table)

# Perform Fisher's Exact Test
test_result <- fisher.test(contingency_table, alternative = "greater")
print("Fisher's Exact Test Result:")
print(test_result)




# 
# ### creating a list of the barcodes for dropouts based on log2 difference
# 
# 
# # pivoting the data
# dropout_barcodes_pivot <- difference_4_drugs %>%
#   pivot_longer(colnames(difference_4_drugs)[2:5]) %>%
#   as.data.frame()
# 
# 
# # plotting the barcodes that have the highest dropouts in all 4 samples
# ggplot(dropout_barcodes_pivot, aes(x = value)) + 
#   geom_histogram(bins = 200) +
#   facet_wrap(~ name) +
#   scale_x_continuous(trans='log10')
# 
# 
# # creating a list of the comparisons
# my_comparisons <- list( c("diference_atr, diference_cis"), c("diference_atr, diference_dox"), c("diference_atr, diference_pf") )
# 
# 
# #
# dropout_barcodes_pivot_filtered <- dropout_barcodes_pivot %>% filter(name %in% c("log2_diff_zscore_atr", "log2_diff_zscore_cis"))
# 
# 
# # desktop
# ggplot(dropout_barcodes_pivot_filtered, aes(x = name, y = value)) + 
#   geom_boxplot() +
#   scale_y_continuous(trans='log10')  +
#   stat_compare_means() +
#   stat_compare_means() 
# 
# 
# ggplot(dropout_barcodes_pivot, aes(x = name, y = value)) + 
#   geom_boxplot() +
#   scale_y_continuous(trans='log10')  +
#   stat_compare_means() +
#   stat_compare_means() 
# 
# 
# 
# # creating a dataframe without the barcodes
# diff_only <- difference_4_drugs[,-1]
# 
# 
# 
# 
# 
# 
