library(ggplot2)
library(dplyr)
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
cpm_scaling <- function(merged_df) {
  
  # Exclude the first column and compute the sum of the remaining columns
  sums <- colSums(merged_df[, -1])
  
  # Divide sums by 1 million
  scaling_factor <- sums / 1000000
  
  # Create a copy of the original dataframe to hold the scaled values
  df <- merged_df
  
  # Apply the scaling factor to all columns except the first
  for (i in 2:ncol(df)) {
    df[, i] <- df[, i] / scaling_factor[i - 1]
  }
  
  # Rename the scaled columns
  scaled_column_names <- paste0(names(df)[-1], "_scaled")
  names(df)[-1] <- scaled_column_names
  
  # Merge the original dataframe with the scaled one
  merged_df <- cbind(merged_df, df[, -1])
  
  # Return the merged dataframe
  return(merged_df)
}


## Log scaling function
compute_log2_scaled <- function(df) {
  
  # Find all column names that contain '_scaled'
  scaled_cols <- grep("_scaled$", names(df), value = TRUE)
  
  # Create log2 columns for each scaled column
  for (col in scaled_cols) {
    log_col <- paste0(col, "_log")
    df[[log_col]] <- log2(df[[col]])
  }
  
  return(df)
}



# Function fror chisquared analysis
compute_chisq_test <- function(df, barcode_col, test_cols, ctrl_cols) {
  # Initialize an empty dataframe to store the results
  p_values <- data.frame(barcode = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Iterate through each unique barcode
  for (barcode_id in unique(df[[barcode_col]])) {
    
    # Extract counts for test and control groups
    barcode_test <- sum(sapply(test_cols, function(col) sum(df[[col]][df[[barcode_col]] == barcode_id])))
    
    barcode_ctrl <- sum(sapply(ctrl_cols, function(col) sum(df[[col]][df[[barcode_col]] == barcode_id])))
    
    ctrl_total <- sum(sapply(ctrl_cols, function(col) sum(df[[col]][df[[barcode_col]] != barcode_id])))
    
    test_total <- sum(sapply(test_cols, function(col) sum(df[[col]][df[[barcode_col]] != barcode_id])))
    
    # Create data frame for chi-squared test
    counts_barcode <- c(barcode_ctrl, barcode_test)
    counts_total <- c(ctrl_total, test_total)
    barcode_df <- data.frame(counts_barcode, counts_total)
    
    # Perform chi-squared test and extract p-value
    chi_sq <- chisq.test(barcode_df)$statistic
    p_value <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
    
    # Append barcode id and p-value to results dataframe
    p_values <- rbind(p_values, data.frame(barcode = barcode_id, p_value = p_value))
  }
  
  return(p_values)
}

# 
# # Cell ID columns should be in the first column labeled V1 
# ctrl_barcode_analysis <- function(ctrl_sample, time_0_barcodes){
#   
#   # changing the barcode to a categorical variable
#   ctrl_sample[['barcode']] <- as.factor(ctrl_sample[['barcode']])
#   
#   
#   # counting the unique barcode counts for the control condition
#   ctrl_sample <- ctrl_sample %>%
#     # grouping the samples by barcode
#     group_by(barcode) %>%
#     # summarizing the 
#     summarize("barcode_count_ctrl_15" = n())
#   
#   # filtering out the barcodes based on T0 barcodes
#   ctrl_sample <- ctrl_sample %>% filter(barcode %in% time_0_barcodes)
#   
#   # returning the ctrl dataframe
#   return(ctrl_sample)
#   
# }
# 
# 
# # function to prep the data for the test condition
# # Cell ID columns should be in the first column labeled V1 
# test_barcode_analysis <- function(test_sample, time_0_barcodes, ctrl_sample){
#   
#   # changing the barcode to a categorical variable
#   test_sample[['barcode']] <- as.factor(test_sample[['barcode']])
#   
#   # counting the unique barcode counts for the control condition
#   test_sample <- test_sample %>%
#     group_by(barcode) %>%
#     summarize("barcode_count_test" = n())
#   
#   # filtering out the barcodes based on T0 barcodes
#   test_sample <- test_sample %>% filter(barcode %in% time_0_barcodes)
#   
#   # merging the control and treated counts
#   test_sample <-  merge(ctrl_sample, test_sample, by='barcode')
#   
#   # Log2 scaling the test counts
#   test_sample <- test_sample %>% mutate(test_log2 = log2(barcode_count_test))
#   
#   # Log2 scaling the ctrl 6 counts
#   test_sample <- test_sample %>% mutate(ctrl_15_log2 = log2(barcode_count_ctrl_15))
#   
#   # returning the test sample dataframe
#   return(test_sample)
#   
# }
# 
# 
# 
#   
# 
# 
