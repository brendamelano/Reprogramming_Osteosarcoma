library(ggrastr)
library(Biostrings)

########    OS384  ATR-i   ##############



# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS384_atr_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_atr_1", "barcode_count_atr_2", "barcode_count_atr_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_atr_final, barcode_col, test_cols, ctrl_cols)

p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")

# Merging the p-values to the final df
OS384_atr_final <- merge(OS384_atr_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_atr_final$logFC <- log2(OS384_atr_final$barcode_mean_atr_cpm / OS384_atr_final$barcode_mean_ctrl_13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1



# Calculate the number of enriched and depleted barcodes
enriched_count <- sum(OS384_atr_final$logFC > 1 & OS384_atr_final$p_value < sig_level)
depleted_count <- sum(OS384_atr_final$logFC < -1 & OS384_atr_final$p_value < sig_level)


# Plot with annotations
volcano_plot <- ggplot(OS384_atr_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 Atr inhibitor", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) +
  # Add annotations for enriched and depleted counts
  annotate("text", x = -6.3, y = 28, label = paste("Depleted:", depleted_count), hjust = 0, size = 2.7) +
  annotate("text", x = 6.3, y = 28, label = paste("Enriched:", enriched_count), hjust = 1, size = 2.7)


volcano_plot

# Save the plot
ggsave("~/Desktop/OS384_atr_volcano_plot.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")



#### CISPLATIN ##



# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS384_cis_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_cis_1", "barcode_count_cis_2", "barcode_count_cis_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_cis_final, barcode_col, test_cols, ctrl_cols)

p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")

OS384_cis_final <- merge(OS384_cis_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_cis_final$logFC <- log2(OS384_cis_final$barcode_mean_cis_cpm / OS384_cis_final$barcode_mean_ctrl_13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


# Calculate the number of enriched and depleted barcodes
enriched_count <- sum(OS384_cis_final$logFC > 1 & OS384_cis_final$p_value < sig_level)
depleted_count <- sum(OS384_cis_final$logFC < -1 & OS384_cis_final$p_value < sig_level)

# Plot with annotations
volcano_plot <- ggplot(OS384_cis_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 Cis treatment", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) +
  # Add annotations for enriched and depleted counts
  annotate("text", x = -6.3, y = 28, label = paste("Depleted:", depleted_count), hjust = 0, size = 2.7) +
  annotate("text", x = 6.3, y = 28, label = paste("Enriched:", enriched_count), hjust = 1, size = 2.7)

volcano_plot

# Save the plot
ggsave("~/Desktop/OS384_cis_volcano_plot.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")



## PF ##



# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS384_pf_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_pf_1", "barcode_count_pf_2", "barcode_count_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_pf_final, barcode_col, test_cols, ctrl_cols)

p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")

OS384_pf_final <- merge(OS384_pf_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_pf_final$logFC <- log2(OS384_pf_final$barcode_mean_pf_cpm / OS384_pf_final$barcode_mean_ctrl_13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1

# Calculate the number of enriched and depleted barcodes
enriched_count <- sum(OS384_pf_final$logFC > 1 & OS384_pf_final$p_value < sig_level)
depleted_count <- sum(OS384_pf_final$logFC < -1 & OS384_pf_final$p_value < sig_level)

volcano_plot <- ggplot(OS384_pf_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 CDK-4/6 inhibitor", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) +
  # Add annotations for enriched and depleted counts
  annotate("text", x = -6.3, y = 28, label = paste("Depleted:", depleted_count), hjust = 0, size = 2.7) +
  annotate("text", x = 6.3, y = 28, label = paste("Enriched:", enriched_count), hjust = 1, size = 2.7)

volcano_plot

# Save the plot
ggsave("~/Desktop/OS384_pf_volcano_plot.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")



#########     IDENTIFYING ENRICHED AND DEPLETED BARCODES    ##



# Identifying the depleted and enriched barcodes
depleted_filtered_atr <- OS384_atr_final %>% filter(logFC < -1 & p_value < 0.05)
depleted_filtered_pf <- OS384_pf_final %>% filter(logFC < -1 & p_value < 0.05)
depleted_filtered_cis <- OS384_cis_final %>% filter(logFC < -1 & p_value < 0.05)


# making a vector of the top dropout barcodes for all the treatments
depleted_barcodes_atr <- depleted_filtered_atr$barcode
depleted_barcodes_pf <- depleted_filtered_pf$barcode
depleted_barcodes_cis <- depleted_filtered_cis$barcode


depleted_barcodes <- c(depleted_barcodes_atr)


# changing the barcode type to DNA in order to later take the reverse complement
depleted_barcodes <- dna(depleted_barcodes)


# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)


# Writing the dropout barcodes to the single cell analysis folder
write.csv(rc_depleted_barcodes, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384_atr_depleted_barcodes_LT.csv")


write.table(rc_depleted_barcodes$rc_depleted_barcodes, file = "~/Desktop/OS384_atr_depleted_barcodes_LT.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)



enriched_depleted_barcodes(
  data = OS384_atr_final,
  depleted_logFC_threshold = -1,
  enriched_logFC_threshold = 1,
  p_value_threshold = 0.05,
  output_file_depleted = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_atr_OS384_LT.txt",
  output_file_enriched = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/enriched_LT_barcodes_atr_OS384_LT.txt"
)

enriched_depleted_barcodes(
  data = OS384_pf_final,
  depleted_logFC_threshold = -1,
  enriched_logFC_threshold = 1,
  p_value_threshold = 0.05,
  output_file_depleted = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_pf_OS384_LT.txt",
  output_file_enriched = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/enriched_LT_barcodes_pf_OS384_LT.txt"
)


enriched_depleted_barcodes(
  data = OS384_cis_final,
  depleted_logFC_threshold = -1,
  enriched_logFC_threshold = 1,
  p_value_threshold = 0.05,
  output_file_depleted = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_cis_OS384_LT.txt",
  output_file_enriched = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/enriched_LT_barcodes_cis_OS384_LT.txt"
)


### Overlapping barcode analysis

library(VennDiagram)
library(grid)


# Read in depleted data frames
depleted_pf <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_pf_OS384_LT.txt", header = FALSE)
depleted_pf <- depleted_pf[,1]

depleted_atr <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_atr_OS384_LT.txt", header = FALSE)
depleted_atr <- depleted_atr[,1]

depleted_cis <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS384/depleted_LT_barcodes_cis_OS384_LT.txt", header = FALSE)
depleted_cis <- depleted_cis[,1]


# Extract the Barcode columns (assuming the column is named "Barcode")
total_pf <- OS384_pf_final$barcode
total_atr <- OS384_atr_final$barcode
total_cis <- OS384_cis_final$barcode

# Create a named list of depleted data frames
depleted_list <- list(
  "CDK-4/6 i depleted" = depleted_pf,
  "ATR i depleted" = depleted_atr,
  "Cisplatin depleted" = depleted_cis
)

# Create a named list of total barcodes
total_list <- list(
  "CDK-4/6 i total" = total_pf,
  "ATR i total" = total_atr,
  "Cisplatin total" = total_cis
)

# Get all combinations of the data frames (pairs)
depleted_pairs <- combn(names(depleted_list), 2, simplify = FALSE)

# Initialize vectors to store p-values and pair identifiers
p_values <- numeric(length(depleted_pairs))
pair_names <- character(length(depleted_pairs))

# Counter for indexing
counter <- 1

# Loop through each pair and perform statistical tests
for (pair in depleted_pairs) {
  # Extract the names and data for the pair
  name1 <- pair[1]
  name2 <- pair[2]
  data1 <- depleted_list[[name1]]
  data2 <- depleted_list[[name2]]
  
  # Extract total barcodes for each treatment
  total1_name <- gsub("depleted", "total", name1)
  total2_name <- gsub("depleted", "total", name2)
  total1 <- total_list[[total1_name]]
  total2 <- total_list[[total2_name]]
  
  # Calculate the total number of unique barcodes across both treatments
  total_barcodes <- union(total1, total2)
  N <- length(total_barcodes)
  
  # Calculate overlaps and counts
  x <- as.numeric(length(intersect(data1, data2)))  # Overlap of depleted barcodes
  b <- as.numeric(length(setdiff(data1, data2)))    # Depleted only in Treatment 1
  c <- as.numeric(length(setdiff(data2, data1)))    # Depleted only in Treatment 2
  d <- as.numeric(N - x - b - c)                    # Not depleted in either treatment
  
  # Create contingency table
  contingency_table <- matrix(c(x, b, c, d), nrow = 2,
                              dimnames = list(
                                c("Depleted in Both", "Not Depleted in Both"),
                                c(name1, name2)
                              ))
  
  # Perform Fisher's Exact Test
  test_result <- fisher.test(contingency_table, alternative = "greater")
  p_value <- test_result$p.value
  
  # Store the p-value and pair name
  p_values[counter] <- p_value
  pair_names[counter] <- paste(name1, "vs", name2, sep = " ")
  
  # Increment counter
  counter <- counter + 1
}

# Apply Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(p_values, method = "BH")


# Create a data frame with results
results <- data.frame(
  Pair = pair_names,
  Raw_P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

# Print the results
print(results)

# Reset counter
counter <- 1

# Loop through each pair again to create Venn diagrams with adjusted p-values
for (pair in depleted_pairs) {
  # Extract the names and data for the pair
  name1 <- pair[1]
  name2 <- pair[2]
  data1 <- depleted_list[[name1]]
  data2 <- depleted_list[[name2]]
  
  # Get the adjusted p-value for this pair
  adjusted_p_value <- adjusted_p_values[counter]
  
  # Generate the Venn diagram with adjusted p-value in the title
  venn.plot <- venn.diagram(
    x = list(data1, data2),
    category.names = c(name1, name2),
    filename = NULL,  # Do not save to file directly
    fill = c("skyblue", "pink"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = c(-20, 20),
    main = paste("OS384 Depleted Barcode Overlap:\n", name1, "vs", name2, "\nAdjusted p-value =", signif(adjusted_p_value, 3)),
    main.cex = 1.2
  )
  
  # Define the output file path for this pair
  # Replace spaces and slashes in names to create valid filenames
  file_name_safe <- function(name) {
    gsub("[ /]", "_", name)
  }
  output_file <- paste0("~/Desktop/", file_name_safe(name1), "_vs_", file_name_safe(name2), "_venn.pdf")
  
  # Save the Venn diagram to a PDF file
  pdf(file = output_file, width = 7, height = 7)
  grid.draw(venn.plot)
  dev.off()
  
  # Increment counter
  counter <- counter + 1
}



########    OS742     ##############



# Initialize empty dataframe to store results
p_values_df <- data.frame(barcode = character(), p_value = numeric())


names(OS742_pf_final)[1] <- 'barcode'

# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_pf_1", "barcode_count_pf_2", "barcode_count_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS742_pf_final, barcode_col, test_cols, ctrl_cols)

p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")

# reassigned this to cis_diff_merged in the forloop above
OS742_pf_final <- merge(OS742_pf_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS742_pf_final$logFC <- log2(OS742_pf_final$barcode_mean_pf_cpm / OS742_pf_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1

volcano_plot <- barcode_volcano_plot(OS742_pf_final, sample_name = "OS742", drug = "CDK-4/6 inhibitor", sig_level = 0.05, fc_cutoff = 1)
print(volcano_plot)


volcano_plot

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
enriched_barcodes <- dna(enriched_barcodes)


# getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)


# Writing the dropout barcodes to the single cell analysis folder

write.table(rc_enriched_barcodes$rc_enriched_barcodes, file = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/15_day_treatments/enriched_barcodes_OS742_PF_LT.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rc_depleted_barcodes$rc_depleted_barcodes, file = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/15_day_treatments/depleted_barcodes_OS742_PF_LT.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)




######     OS052    ########


## ATR ##

# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS052_atr_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_atr_1", "barcode_count_atr_2", "barcode_count_atr_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS052_atr_final, barcode_col, test_cols, ctrl_cols)


p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")


# Reassigned this to cis_diff_merged in the forloop above
OS052_atr_final <- merge(OS052_atr_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
# test difference of logs as well
OS052_atr_final$logFC <- log2(OS052_atr_final$barcode_cpm_mean_atr / OS052_atr_final$barcode_mean_ctrl13_cpm)


write.csv(OS052_atr_final, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/processed_counts/OS052_atr_final.csv", row.names = TRUE)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- barcode_volcano_plot(OS052_atr_final, sample_name = "OS052", drug = "Atr inhibitor", sig_level = 0.05, fc_cutoff = 1)
print(volcano_plot)


# Save the plot
ggsave("~/Desktop/OS052_atr_volcano_plot.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")


##### PF    ###


# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS052_pf_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_pf_1", "barcode_count_pf_2", "barcode_count_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS052_pf_final, barcode_col, test_cols, ctrl_cols)


# Adjusting the p-values based on BH correction
p_values_df$adjusted_p_value <- p.adjust(p_values_df$p_value, method = "BH")


# reassigned this to cis_diff_merged in the forloop above
OS052_pf_final <- merge(OS052_pf_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
# test difference of logs as well
OS052_pf_final$logFC <- log2(OS052_pf_final$barcode_mean_pf_cpm / OS052_pf_final$barcode_mean_ctrl_13_cpm)


write.csv(OS052_pf_final, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_gDNA_barcodes/processed_counts/OS052_pf_final.csv", row.names = TRUE)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- barcode_volcano_plot(OS052_pf_final, sample_name = "OS052", drug = "CDK-4/6 inhibitor", sig_level = 0.05, fc_cutoff = 1)
print(volcano_plot)


# Saving the PF plot
ggsave("~/Desktop/OS052_PF_volcano.svg", plot = volcano_plot, width = 2.2, height = 2.2, units = "in")



##





enriched_depleted_barcodes(
  data = OS052_pf_final,
  depleted_logFC_threshold = -1,
  enriched_logFC_threshold = 1,
  p_value_threshold = 0.05,
  output_file_depleted = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/depleted_LT_barcodes_pf_OS052_LT.txt",
  output_file_enriched = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/enriched_LT_barcodes_pf_OS052_LT.txt"
)






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



