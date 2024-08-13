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

process_and_filter_barcodes <- function(input_df, sample_name, time_0_barcodes) {
  
  # Renaming the first column to 'barcode'
  names(input_df)[1] <- 'barcode'
  
  # Changing the column names for the scaled counts based on the sample name
  for (i in 1:3) {
    names(input_df)[i + 1] <- paste0("barcode_count_", sample_name, "_", i)
  }
  
  # Filtering based on the time 0 barcodes
  filtered_df <- input_df %>% dplyr::filter(barcode %in% time_0_barcodes)
  
  return(filtered_df)
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
