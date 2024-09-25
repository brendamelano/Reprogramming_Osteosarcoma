library(ggplot2)
library(dplyr)
library(tidyr)
library(dplyr)
library(stringdist)
library(data.table)


# Function to read in and process the files
process_file <- function(file_path) {
  
  # extracting the basename from the file path
  file_name <- basename(file_path)
  
  # extracting the sequence after the number prefix
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


# Filtering based on time 0 barcodes
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

process_and_filter_barcodes <- function(input_df, sample_name, time_0_barcodes) {

  # Renaming the first column to 'barcode'
  names(input_df)[1] <- 'barcode'
  
  # Changing the column names for the scaled counts based on the sample name
  for (i in 1:3) {
    names(input_df)[i + 1] <- paste0("barcode_count_", sample_name, "_", i)
  }
  
  # Identify barcodes in time_0_barcodes and not in time_0_barcodes
  in_whitelist_df <- input_df %>% filter(barcode %in% time_0_barcodes)
  not_in_whitelist_df <- input_df %>% filter(!(barcode %in% time_0_barcodes))
  
  # Get the barcodes
  barcodes_in_whitelist <- as.character(in_whitelist_df$barcode)
  barcodes_not_in_whitelist <- as.character(not_in_whitelist_df$barcode)
  
  # Define the alphabet
  alphabet <- c("A", "C", "G", "T", "N")
  
  # Function to generate barcodes at Hamming distance one
  generate_hamming_distance_one_barcodes <- function(barcode, alphabet) {
    barcode_chars <- unlist(strsplit(barcode, split=""))
    positions <- seq_along(barcode_chars)
    hamming_barcodes <- c()
    for (pos in positions) {
      for (letter in alphabet) {
        if (letter != barcode_chars[pos]) {
          new_barcode_chars <- barcode_chars
          new_barcode_chars[pos] <- letter
          new_barcode <- paste0(new_barcode_chars, collapse="")
          hamming_barcodes <- c(hamming_barcodes, new_barcode)
        }
      }
    }
    return(hamming_barcodes)
  }
  
  # Create a mapping from barcodes at Hamming distance one to the original barcode
  hamming_mapping <- data.frame(generated_barcode=character(), original_barcode=character(), stringsAsFactors=FALSE)
  
  for (barcode in barcodes_in_whitelist) {
    hamming_barcodes <- generate_hamming_distance_one_barcodes(barcode, alphabet)
    hamming_mapping <- rbind(hamming_mapping, data.frame(generated_barcode=hamming_barcodes, original_barcode=barcode, stringsAsFactors=FALSE))
  }
  
  # Proceed only if there are barcodes not in the whitelist and the hamming mapping is not empty
  if (nrow(not_in_whitelist_df) > 0 && nrow(hamming_mapping) > 0) {
    # Merge non-whitelist barcodes with the Hamming mapping
    merged_df <- merge(not_in_whitelist_df, hamming_mapping, by.x="barcode", by.y="generated_barcode")
    
    if (nrow(merged_df) > 0) {
      # Define counts columns
      counts_columns <- names(not_in_whitelist_df)[2:4]
      
      # Sum the counts for each original_barcode
      additional_counts <- merged_df %>% group_by(original_barcode) %>% summarise(across(all_of(counts_columns), sum, .names = "{.col}_additional"))
      
      # Merge additional_counts into in_whitelist_df
      in_whitelist_df <- in_whitelist_df %>% left_join(additional_counts, by=c("barcode"="original_barcode"))
      
      # Replace NA in additional counts with 0 and sum the counts
      for (col in counts_columns) {
        additional_col <- paste0(col, "_additional")
        if (!(additional_col %in% names(in_whitelist_df))) {
          # If the additional counts column doesn't exist, create it with zeros
          in_whitelist_df[[additional_col]] <- rep(0, nrow(in_whitelist_df))
        } else {
          in_whitelist_df[[additional_col]][is.na(in_whitelist_df[[additional_col]])] <- 0
        }
        in_whitelist_df[[col]] <- in_whitelist_df[[col]] + in_whitelist_df[[additional_col]]
      }
      
      # Remove the additional counts columns
      in_whitelist_df <- in_whitelist_df[, !(names(in_whitelist_df) %in% paste0(counts_columns, "_additional"))]
    }
  }
  
  return(in_whitelist_df)
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



log_scales_and_means <- function(data, group_prefix) {
  # Compute log2 scaled values
  scaled_log_data <- compute_log2_scaled(data)
  
  # Dynamically create column names based on the group prefix
  log_columns <- paste0("barcode_count_", group_prefix, "_", 1:3, "_scaled_log")
  scaled_columns <- paste0("barcode_count_", group_prefix, "_", 1:3, "_scaled")
  
  # Compute the mean per barcode for the log-scaled columns and add the group prefix to the variable name
  scaled_log_data <- scaled_log_data %>%
    mutate(!!paste0("barcode_mean_", group_prefix, "_log") := rowMeans(select(., all_of(log_columns))))
  
  # Compute the mean per barcode for the original scaled columns and add the group prefix to the variable name
  scaled_log_data <- scaled_log_data %>%
    mutate(!!paste0("barcode_mean_", group_prefix, "_cpm") := rowMeans(select(., all_of(scaled_columns))))
  
  return(scaled_log_data)
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


barcode_volcano_plot <- function(data, sample_name, drug, sig_level = 0.05, fc_cutoff = 1) {
  
  # Calculate the number of enriched and depleted barcodes
  enriched_count <- sum(data$logFC > fc_cutoff & data$adjusted_p_value < sig_level)
  depleted_count <- sum(data$logFC < -fc_cutoff & data$adjusted_p_value < sig_level)
  
  # Create the volcano plot
  volcano_plot <- ggplot(data, aes(x = logFC, y = -log10(adjusted_p_value))) +
    geom_point_rast(size = 0.5, aes(color = ifelse(adjusted_p_value < sig_level & (logFC > fc_cutoff | logFC < -fc_cutoff), "red", "black")), show.legend = FALSE) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = paste(sample_name, drug), x = "logFC", y = "-log10(p-value)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 8)) + 
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(sig_level), linetype = "dashed", color = "gray") +
    ylim(0, 30) +
    
    # Add annotations for enriched and depleted counts
    annotate("text", x = -6.3, y = 28, label = paste("Depleted:", depleted_count), hjust = 0, size = 2.7) +
    annotate("text", x = 6.3, y = 28, label = paste("Enriched:", enriched_count), hjust = 1, size = 2.7)
  
  return(volcano_plot)
}

