library(stringdist)
library(tidyverse)
library(ggplot2)
library(ggpubr)
#library(VennDiagram)
library(stats)
library(dplyr)
library(tidyr)
library(mgcv) # GLMGAM regression
library(purrr)
library(grid)
library(png)
#library(DescTools)



############      PROCESSING SAMPLES FOR CTRL D13      #####################



# Creating the file paths to read in
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/13_384_ctrl13_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/14_384_ctrl13_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/15_384_ctrl13_3.fastq.gz.out2.txt')


# Applying the file_paths function to read in and process the files
result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_ctrl_13_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Renaming the first column to barcode
names(OS384_ctrl_13_merged)[1] <- 'barcode'


# Changing the column names for the scaled counts
names(OS384_ctrl_13_merged)[2:4] <- c("barcode_count_ctrl_13_1", 
                                     "barcode_count_ctrl_13_2", 
                                     "barcode_count_ctrl_13_3")



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


# Renaming the first column to barcode
names(OS384_atr_merged)[1] <- 'barcode'


# Merging the control and treated counts
OS384_atr_ctrl13 <-  merge(OS384_ctrl_13_merged, OS384_atr_merged, by='barcode')


# Creating a dataframe for the barcodes not in the white list
test_sample_extra <- OS384_atr_ctrl13 %>% filter(!(barcode %in% time_0_barcodes))
test_sample <- OS384_atr_ctrl13 %>% filter((barcode %in% time_0_barcodes))


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


# Reassigning the test sample
OS384_atr_final <- test_sample


# Performing cpm scaling with the function
OS384_atr_scaled <- cpm_scaling(OS384_atr_final)


# 
OS384_atr_log_scaled <- compute_log2_scaled(OS384_atr_scaled)


# Computing the mean per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_atr_cpm = rowMeans(select(., c("barcode_count_384_atr_1_scaled", 
                                                 "barcode_count_384_atr_2_scaled", 
                                                 "barcode_count_384_atr_3_scaled"))))


# Computing the mean cpm per barcode for the merged dataframe
OS384_atr_log_scaled <- OS384_atr_log_scaled %>% 
  mutate(barcode_mean_ctrl13_cpm = rowMeans(select(., c("barcode_count_ctrl_13_1_scaled", 
                                                        "barcode_count_ctrl_13_2_scaled", 
                                                        "barcode_count_ctrl_13_3_scaled"))))




###

#OS384_atr_final <- test_barcode_analysis(test_sample = OS384_atr, time_0_barcodes, ctrl_sample = OS384ctrl_6_final)

## Finishing the run of the function

## Identifying enriched and depleted barcodes by taking the difference of log and z-score of the differences


# Taking the difference of the log transformed values
atr_diff_merged <- OS384_atr_final %>% mutate(difference_atr_log2 = barcode_log_mean_ctrl_13 - barcode_mean_atr)


### PLOTTING THE REPLICATES
#Perform regression analysis
model <- lm(barcode_count_384_atr_1_scaled_log ~ barcode_count_384_atr_3_scaled_log, data = OS384_atr_log_scaled)


# Extract r-squared and p-value
r_squared <- summary(model)$r.squared


# Create the ggplot
OS384_atr_replicate <- ggplot(OS384_atr_log_scaled, aes(barcode_count_384_atr_3_scaled_log, barcode_count_384_atr_1_scaled_log)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),  # Change the font size for axis titles
        plot.title = element_text(size = 9)) +  # Change the font size for the main title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Log Transformed Barcode Count - Replicate 1") +
  ylab("Log Transformed Barcode Count - Replicate 2") +
  ggtitle("OS384 ATR Barcode Count Correlation") +
  geom_text(x = min(OS384_atr_log_scaled$barcode_count_384_atr_3_scaled_log),
            y = max(OS384_atr_log_scaled$barcode_count_384_atr_1_scaled_log),
            label = paste("R-squared =", round(r_squared, 2), "\n"),
            hjust = 0, vjust = 1, parse = TRUE)


# Save the plot as an SVG file
ggsave("~/Desktop/OS384_atr_replicate.png", plot = OS384_atr_replicate,  device = "png", width = 3.2, height = 3.2)


# Plot the data
ggplot(filtered_data, aes(x = barcode)) + 
  geom_line(aes(y = barcode_log_mean_ctrl_13, color = "Control")) +
  geom_line(aes(y = barcode_mean_atr, color = "Treatment")) +
  labs(title = "Barcode Counts",
       y = "Count",
       color = "Condition") +
  theme_minimal()



###########    PF (PFIZER) BARCODES    ##################



# Reading in the PF barcodes
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/22_384_pf_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/23_384_pf_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/24_384_pf_3.fastq.gz.out2.txt')


result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_pf_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Changing the barcode column name
colnames(OS384_pf_merged)[1] <- "barcode"


# Performing cpm scaling with the function
OS384_pf_scaled <- cpm_scaling(OS384_pf_merged)


# 
OS384_pf_log_scaled <- compute_log2_scaled(OS384_pf_scaled)



# Chaning the name of the V1 column
names(OS384_pf_log_scaled)[1] <- "barcode"



# Computing the mean per barcode for the merged dataframe
# Make sure that this does not have to be the mean
OS384_pf_log_scaled <- OS384_pf_log_scaled %>% 
  mutate(barcode_mean_pf_log = rowMeans(select(., c("barcode_count_384_pf_1_scaled_log", 
                                                    "barcode_count_384_pf_2_scaled_log", 
                                                    "barcode_count_384_pf_3_scaled_log"))))


# Filtering out the barcodes based on T0 barcodes
OS384_pf_final <- OS384_pf_log_scaled %>% filter(barcode %in% time_0_barcodes)



#test_sample <- OS384_pf


# running the function
#OS384_pf <- test_barcode_analysis(test_sample, time_0_barcodes)


# merging atr_diff_ordered with ctrl day 0 barcode counts by barcode
pf_diff_ordered <- merge(pf_diff_ordered, OS384_ctrl_0_unique, by = "barcode")


# merging atr_diff_ordered with ctrl day 0 barcode counts
pf_diff_ordered <- merge(pf_diff_ordered, OS384_ctrl_0_unique, by = "barcode")


# plotting 
ggplot(pf_diff_ordered, aes(x = barcode_count_ctrl_0.y , y = log2_diff_zscore_pf )) + 
  geom_point() +
  scale_x_continuous(trans='log10')



###############     CIS BARCODES    ############



# Reading in the Cisplatin treated counts
file_paths <- c('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/16_384_cis_1.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/17_384_cis_2.fastq.gz.out2.txt',
                '~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_74_counts/18_384_cis_3.fastq.gz.out2.txt')


result_list <- lapply(file_paths, process_file)


# Merging the data frames by 'V1'
OS384_cis_merged <- Reduce(function(x, y) merge(x, y, by = "V1"), result_list)


# Performing cpm scaling with the function
OS384_cis_scaled <- cpm_scaling(OS384_cis_merged)


# Renaming the columns with the raw values
names(OS384_cis_scaled)[2:4] <- c("barcode_count_cis_1", "barcode_count_cis_2", "barcode_count_cis_3")


# Changing the colun names for the scaled counts
names(OS384_cis_scaled)[5:7] <- c("barcode_count_cis_1_scaled", "barcode_count_cis_2_scaled", "barcode_count_cis_3_scaled")


# Computing logs of cpm values
OS384_cis_log_scaled <- OS384_cis_scaled %>% mutate(barcode_count_cis_1_log = log2(barcode_count_cis_1_scaled))
OS384_cis_log_scaled <- OS384_cis_log_scaled %>% mutate(barcode_count_cis_2_log = log2(barcode_count_cis_2_scaled))
OS384_cis_log_scaled <- OS384_cis_log_scaled %>% mutate(barcode_count_cis_3_log = log2(barcode_count_cis_3_scaled))


# computing the mean per barcode for the merged dataframe
# Make sure that this does not have to be the mean
OS384_cis_log_scaled <- OS384_cis_log_scaled %>% 
  mutate(barcode_mean_cis_log = rowMeans(select(., c("barcode_count_cis_1_log", 
                                                     "barcode_count_cis_2_log", 
                                                     "barcode_count_cis_3_log"))))


OS384_cis_log_scaled <- OS384_cis_log_scaled %>% 
  mutate(barcode_mean_cis_cpm = rowMeans(select(., c("barcode_count_cis_1_scaled", 
                                                     "barcode_count_cis_2_scaled", 
                                                     "barcode_count_cis_3_scaled"))))

names(OS384_cis_log_scaled)[1] <- "barcode"

# filtering out the barcodes based on T0 barcodes
test_sample <- OS384_cis_log_scaled %>% filter(barcode %in% time_0_barcodes)
ctrl_sample <- OS384_cis_log_scaled %>% filter(barcode %in% time_0_barcodes)


# Merging the control and treated counts
test_sample <-  merge(ctrl_sample, test_sample, by='barcode')


# Reassigning the test_sample
OS384_cis_final <- test_sample


# Taking the difference of the log transformed values
cis_diff_merged <- OS384_cis_final %>% mutate(difference_cis_log2 = barcode_log_mean_ctrl_13 - barcode_mean_cis_log)


# Scaling the data    
cis_diff_merged <- cis_diff_merged %>% mutate(log2_diff_zscore_cis = scale(difference_cis_log2))


# merging atr_diff_ordered with ctrl day 0 barcode counts by barcode
cis_diff_merged <- merge(cis_diff_merged, OS384ctrl0_log_scaled, by = "barcode")


## performing logistic regression to get the p-values

# filtering out the barcodes based on T0 barcodes
test_sample <- OS384_cis_log_scaled %>% filter(barcode %in% time_0_barcodes)
ctrl_sample <- OS384ctrl13_log_scaled %>% filter(barcode %in% time_0_barcodes)


# 
OS384_cis_logReg <- OS384_cis_merged[,-(2:4)]
names(OS384_cis_logReg)[2:4] <- c("count_scaled_1", "count_scaled_2", "count_scaled_3")
OS384_cis_logReg$treated <- 1


# control D13 samples
OS384_ctrl13_logReg <- OS384_ctrl13_scaled[,-(2:4)]
names(OS384_ctrl13_logReg)[2:4] <- c("count_scaled_1", "count_scaled_2", "count_scaled_3")
OS384_ctrl13_logReg$treated <- 0


# binding the cis and ctrl dataframe
OS384_logReg_df <- rbind(OS384_cis_logReg, OS384_ctrl13_logReg)



########      QC for all samples    ##########



# Subset the dataframe to exclude columns with 'log', 'scaled', 'mean'
raw_counts_df_pf <- OS384_pf_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))


raw_counts_df_atr <- OS384_atr_final %>%
  select(-contains("log"), -contains("scaled"), 
         -contains("mean"), -contains("StdDev"), -contains("Index"))

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
ggplot(data_long_summary, aes(x = condition, y = counts)) +
  geom_jitter(aes(color = "Data Points"), width = 0.2, height = 0, alpha = 0.5) + # Actual data points
  geom_errorbar(aes(ymin = lower, ymax = upper, x = condition), width = 0.2) + # Error bars
  geom_point(aes(y = mean, color = "Mean"), size = 3) + # Mean points
  geom_line(aes(y = median, group = condition, color = "Median"), size = 1) + # Optional: Line for median
  scale_color_manual("", values = c("Data Points" = "black", "Mean" = "red", "Median" = "blue")) +
  theme_bw() + # Using theme_bw
  labs(title = "Summary Statistics of Barcode Counts with Data Points",
       x = "Condition",
       y = "Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") # Improve readability





##########    IDENTIFYING THE OVERLAPPING BARCODE DROPOUTS   #############



# Concatenating the barcode lists and keeping the unique barcodes
deplted_barcodes <- unique(c(cis_barcodes, pf_barcodes, atr_barcodes))


# changing the barcode type to DNA in order to later take the reverse complement
deplted_barcodes <- dna(deplted_barcodes)


# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(deplted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)


# Writing the dropout barcodes to thw single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/depleted_barcodes_OS384_inVivo_LT.csv")



enriched_barcodes <- unique(c(cis_barcodes, pf_barcodes, atr_barcodes))
enriched_barcodes <- dna(enriched_barcodes)
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)
write.csv(rc_enriched_barcodes, "~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/enriched_barcodes_OS384_inVivo_LT.csv")



# Reading in the fastq files for read 1 and read 2 to combine the reads from the scRNAseq data
read1_file <- read.delim("/Users/bmelano/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/384-in-vivo_S1_L001_R1_001.fastq")
read2_file <- read.delim("/Users/bmelano/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/384-in-vivo_S1_L001_R2_001.fastq")


# finding the barcodes that overlapped in selected barcodes for all samples
venn.diagram(
  x = list(cis_barcodes, pf_barcodes, atr_barcodes),
  category.names = c("Cisplatin" , "CDK 4/6 i", "ATR i"),
  filename = '~/Desktop/overlapping_barcodes.svg',
  output=TRUE,
  height = 2150,
  width = 2150,
  main = "OS384 Barcode Selection Overlap"
)


# Making a vector of the barcodes for all the treatments
cis_barcodes_all <- OS384_cis_final$barcode
pf_barcodes_all <- OS384_pf_final$barcode
atr_barcodes_all <- OS384_atr_final$barcode


# finding the barcodes that overlapped in selected barcodes for all samples
venn.diagram(
  x = list(cis_barcodes_all, pf_barcodes_all, atr_barcodes_all),
  category.names = c("Cisplatin" , "CDK 4/6 i", "ATR i"),
  filename = '~/Desktop/overlapping_barcodes.svg',
  output=TRUE,
  height = 2150,
  width = 2150,
  main = "OS384 Barcode Selection Overlap"
)

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
