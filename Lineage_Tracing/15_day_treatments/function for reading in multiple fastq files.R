library(dplyr)
library(ggplot2)
library(tidyr)


# Function to read in and process the files
process_file <- function(file_path) {
  
  # extracting the basename from the file path
  file_name <- basename(file_path)
  
  # extracing the sequence after the number prefix
  file_name <- substring(file_name, 4)
  
  # removing the suffix from the base names
  column_name <- paste0("barcode_count_", gsub(".fastq.gz.out2.txt", "", file_name))
  
  # reading in the files in the file path
  data <- read.delim(file_path, header = FALSE)
  
  # counting the number of each barcode
  data_summary <- data %>%
    group_by(V1) %>%
    summarize(!!column_name := n())
  
  # returning the data summary
  return(data_summary)
}


# Creating a function to compute scaled values of merged dataframe
cpm_scaling <- function(merged_df){
  
  # renaming the scaled dataframe
  colnames(merged_df)[1] = "barcode"
  
  # computing the sum of the count columns
  sums <- colSums(merged_df[2:4])
  
  scaling_factor <- sums / 1000000
  
  df <- merged_df
  
  for (i in 2:ncol(df)) {
    df[, i] <- df[, i] / scaling_factor[i - 1]
  }
  
  # merging the original dataframe with the scaled df
  merged_df <- merge(merged_df, df, by = "barcode")
  
  # returning the merged dataframe
  return(merged_df)
}


# Function for creating log values and average log and cpm values
creating_logs <- function(cpm_scaled_data){
  
  
}


# Cell ID columns should be in the first column labeled V1 
ctrl_barcode_analysis <- function(ctrl_sample, time_0_barcodes){
  
  # changing the barcode to a categorical variable
  ctrl_sample[['barcode']] <- as.factor(ctrl_sample[['barcode']])
  
  
  # counting the unique barcode counts for the control condition
  ctrl_sample <- ctrl_sample %>%
    # grouping the samples by barcode
    group_by(barcode) %>%
    # summarizing the 
    summarize("barcode_count_ctrl_15" = n())
  
  # filtering out the barcodes based on T0 barcodes
  ctrl_sample <- ctrl_sample %>% filter(barcode %in% time_0_barcodes)
  
  # returning the ctrl dataframe
  return(ctrl_sample)
  
}


# function to prep the data for the test condition
# Cell ID columns should be in the first column labeled V1 
test_barcode_analysis <- function(test_sample, time_0_barcodes, ctrl_sample){
  
  # changing the barcode to a categorical variable
  test_sample[['barcode']] <- as.factor(test_sample[['barcode']])
  
  # counting the unique barcode counts for the control condition
  test_sample <- test_sample %>%
    group_by(barcode) %>%
    summarize("barcode_count_test" = n())
  
  # filtering out the barcodes based on T0 barcodes
  test_sample <- test_sample %>% filter(barcode %in% time_0_barcodes)
  
  # merging the control and treated counts
  test_sample <-  merge(ctrl_sample, test_sample, by='barcode')
  
  # Log2 scaling the test counts
  test_sample <- test_sample %>% mutate(test_log2 = log2(barcode_count_test))
  
  # Log2 scaling the ctrl 6 counts
  test_sample <- test_sample %>% mutate(ctrl_15_log2 = log2(barcode_count_ctrl_15))
  
  # returning the test sample dataframe
  return(test_sample)
  
}


# computing the difference of the z scores
paste0(test_condition, "_diff") <- paste0(test_condition, "_diff") %>% mutate(paste0(test_condition, "_diff_zscores") = paste0("barcode_count_", ctrl_condition) - paste0("barcode_count_", test_condition))


# ordering based on the difference of z scores
paste0(test_condition, "_diff_ordered") <- paste0(test_condition, "_diff")[order(paste0(test_condition, "_diff")$paste0("difference_", test_condition, "_zscores"), decreasing = T), ]


# filtering the z-scores to keep the top dropouts
z_score_diff_filtered_atr <- atr_diff_ordered[1:100,]




  


