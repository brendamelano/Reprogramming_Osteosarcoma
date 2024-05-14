


########    OS384     ##############



# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS384_atr_log_scaled)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_384_atr_1", "barcode_count_384_atr_2", "barcode_count_384_atr_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_atr_log_scaled, barcode_col, test_cols, ctrl_cols)


OS384_atr_log_scaled <- merge(OS384_atr_log_scaled, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_atr_log_scaled$logFC <- log2(OS384_atr_log_scaled$barcode_mean_atr_cpm / OS384_atr_log_scaled$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- ggplot(OS384_atr_log_scaled, aes(x=logFC, y=-log10(p_value))) +
  geom_point(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 Atr-i treatment", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8.5)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) 


# Save the plot
ggsave("~/Desktop/OS384_atr_volcano_plot.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")


## Idenitfying depleted and enriched barcodes

# Identifying the depleted and enriched barcodes
# filtering based on log fold change to identify positive fold change
depleted_filtered_atr <- OS384_atr_log_scaled %>% filter(logFC < -1 & p_value < 0.05)


# making a vector of the top dropout barcodes for all the treatments
depleted_barcodes <- depleted_filtered_atr$barcode


# changing the barcode type to DNA in order to later take the reverse complement
depleted_barcodes <- dna(depleted_barcodes)

# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)


# Writing the dropout barcodes to the single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_atr_depleted_barcodes_LT.csv")


########    OS742     ##############


# Initialize empty dataframe to store results
p_values_df <- data.frame(barcode = character(), p_value = numeric())

names(OS742_pf_log_scaled)
names(OS742_pf_log_scaled)[1] <- 'barcode'

# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_42_pf_1", "barcode_count_42_pf_2", "barcode_count_42_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS742_pf_log_scaled, barcode_col, test_cols, ctrl_cols)


# reassigned this to cis_diff_merged in the forloop above
OS742_pf_final <- merge(OS742_pf_log_scaled, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS742_pf_final$logFC <- log2(OS742_pf_final$barcode_mean_pf_cpm / OS742_pf_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs

sig_level <- 0.05
fc_cutoff <- 1


# Create the volcano plot
volcano_plot <- ggplot(OS742_pf_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point(size=0.5, aes(color=ifelse(p_value<sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS742 CDK 4/6-i treatment", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) +
  ylim(0, 8) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  xlim(-2.5, 2.5)


# Save the plot
ggsave("~/Desktop/OS742_cdk_volcano_plot.svg", plot = volcano_plot, width = 2.4, height = 2.4, units = "in")


## ENRICHED AND DEPLETED BARCODES




# Identifying the depleted and enriched barcodes
# filtering based on log fold change to identify positive fold change
enriched_filtered_pf <- OS742_pf_final %>% filter(logFC > 1 & p_value < 0.05)
depleted_filtered_pf <- OS742_pf_final %>% filter(logFC < -1 & p_value < 0.05)


# making a vector of the top dropout barcodes for all the treatments
depleted_barcodes <- depleted_filtered_pf$barcode
enriched_barcodes <- enriched_filtered_pf$barcode


# changing the barcode type to DNA in order to later take the reverse complement
depleted_barcodes <- dna(depleted_barcodes)


# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)

# Writing the dropout barcodes to the single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/Osteo_Lineage_Tracing_Analysis/15_day_treatments/depleted_barcodes_OS742_inVivo_LT.csv")



enriched_barcodes <- unique(c(cis_barcodes, pf_barcodes, atr_barcodes))

enriched_barcodes <- dna(enriched_barcodes)
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)
write.csv(rc_enriched_barcodes, "~/Desktop/Osteo_Lineage_Tracing_Analysis/15_day_treatments/enriched_barcodes_OS742_inVivo_LT.csv")



######     OS052    ########



# function for performing chi-squared analysis
perform_chi_squared_analysis <- function(data, drug, barcode_colname, test_colnames, ctrl_colnames) {
  
  # Initialize an empty data frame to store results
  p_values <- data.frame(barcode = character(), p_value = numeric())
  
  # Loop through unique barcode IDs
  for (barcode_id in unique(data[[barcode_colname]])) {
    # Extract counts for test and control groups dynamically based on column names
    barcode_test <- sum(data[data[[barcode_colname]] == barcode_id, test_colnames])
    barcode_ctrl <- sum(data[data[[barcode_colname]] == barcode_id, ctrl_colnames])
    
    # Calculate totals for control and test excluding the current barcode
    ctrl_total <- sum(data[data[[barcode_colname]] != barcode_id, ctrl_colnames])
    test_total <- sum(data[data[[barcode_colname]] != barcode_id, test_colnames])
    
    # Create data frame for chi-squared test
    counts_barcode <- c(barcode_ctrl, barcode_test)
    counts_total <- c(ctrl_total, test_total)
    barcode_df <- data.frame(counts_barcode, counts_total)
    
    # Perform chi-squared test and extract p-value
    chi_sq_test_result <- chisq.test(barcode_df)
    chi_sq <- chi_sq_test_result$statistic
    p_value <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
    
    # Append barcode ID and p-value to results data frame
    p_values <- rbind(p_values, data.frame(barcode = barcode_id, p_value = p_value))
  }
  
  return(p_values)
}


atr_p_values <- perform_chi_squared_analysis(
  data = OS052_atr_final,
  drug = "atr",
  barcode_colname = "barcode",
  test_colnames = c("barcode_count_atr_1", "barcode_count_atr_2", "barcode_count_atr_3"),
  ctrl_colnames = c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")
)


# reassigned this to cis_diff_merged in the forloop above
OS052_atr_final <- merge(OS052_atr_final, atr_p_values, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS052_atr_final$logFC <- log2(OS052_atr_final$barcode_cpm_mean_atr / OS052_atr_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


# Create the volcano plot
ggplot(OS052_atr_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point(size=0.5, aes(color=ifelse(p_value<sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS052 ATR inhibitor Treated Barcode Fold change", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) +
  #ylim(0, 8) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray")# +
  #xlim(-2.5, 2.5)



# Identifying the depleted and enriched barcodes
# filtering based on log fold change to identify positive fold change
depleted_filtered_atr <- OS052_atr_final %>% filter(logFC < -1 & p_value < 0.05)


# making a vector of the top dropout barcodes for all the treatments
depleted_barcodes <- depleted_filtered_atr$barcode
#enriched_barcodes <- depleted_filtered_atr$barcode


# changing the barcode type to DNA in order to later take the reverse complement
depleted_barcodes <- dna(depleted_barcodes)

rc_depleted_barcodes <- depleted_barcodes

# getting the reverse complement of the top barcodes
#rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)


# Writing the dropout barcodes to the single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/depleted_barcodes_OS052_LT.csv")



