
library(ggrastr)

########    OS384     ##############



# Initialize empty dataframe to store results
p_values <- data.frame(barcode = character(), p_value = numeric())


# Printing names to define test and ctrl cols
names(OS384_atr_final)


# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_384_atr_1", "barcode_count_384_atr_2", "barcode_count_384_atr_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_atr_final, barcode_col, test_cols, ctrl_cols)


OS384_atr_final <- merge(OS384_atr_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_atr_final$logFC <- log2(OS384_atr_final$barcode_mean_atr_cpm / OS384_atr_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- ggplot(OS384_atr_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 ATR inhibitor", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8.5)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) 


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
test_cols <- c("barcode_count_384_cis_1", "barcode_count_384_cis_2", "barcode_count_384_cis_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_cis_final, barcode_col, test_cols, ctrl_cols)


OS384_cis_final <- merge(OS384_cis_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_cis_final$logFC <- log2(OS384_cis_final$barcode_mean_cis_cpm / OS384_cis_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- ggplot(OS384_cis_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point(size=0.5, aes(color=ifelse(p_value < sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS384 Cis treatment", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8.5)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30) 

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
test_cols <- c("barcode_count_384_pf_1", "barcode_count_384_pf_2", "barcode_count_384_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS384_pf_final, barcode_col, test_cols, ctrl_cols)


OS384_pf_final <- merge(OS384_pf_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
OS384_pf_final$logFC <- log2(OS384_pf_final$barcode_mean_pf_cpm / OS384_pf_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


volcano_plot <- ggplot(OS384_pf_final, aes(x=logFC, y=-log10(p_value))) +
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



#### IDENTIFYING ENRICHED AND DEPLETED BARCODES ##



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


########    OS742     ##############


# Initialize empty dataframe to store results
p_values_df <- data.frame(barcode = character(), p_value = numeric())

names(OS742_pf_final)
names(OS742_pf_final)[1] <- 'barcode'

# Define the column names
barcode_col <- "barcode"
test_cols <- c("barcode_count_42_pf_1", "barcode_count_42_pf_2", "barcode_count_42_pf_3")
ctrl_cols <- c("barcode_count_ctrl_13_1", "barcode_count_ctrl_13_2", "barcode_count_ctrl_13_3")


# Compute chi-squared tests
p_values_df <- compute_chisq_test(OS742_pf_final, barcode_col, test_cols, ctrl_cols)


# reassigned this to cis_diff_merged in the forloop above
OS742_pf_final <- merge(OS742_pf_final, p_values_df, by = 'barcode')


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

# Exporting enriched barcodes
enriched_barcodes <- unique(c(cis_barcodes, pf_barcodes, atr_barcodes))

enriched_barcodes <- dna(enriched_barcodes)
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)
write.csv(rc_enriched_barcodes, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/scRNAseq_LT_analysis/OS742/enriched_barcodes_OS742_inVivo_LT.csv")



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


# Reassigned this to cis_diff_merged in the forloop above
OS052_atr_final <- merge(OS052_atr_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
# test difference of logs as well
OS052_atr_final$logFC <- log2(OS052_atr_final$barcode_cpm_mean_atr / OS052_atr_final$barcode_mean_ctrl13_cpm)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1


# Create the volcano plot
volcano_plot <- ggplot(OS052_atr_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value<sig_level & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS052 ATR inhibitor", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8.5)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30)


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


# reassigned this to cis_diff_merged in the forloop above
OS052_pf_final <- merge(OS052_pf_final, p_values_df, by = 'barcode')


# computing log fold change for different samples
# try computing with log values to see if the values in the middle with high p-values change
# test difference of logs as well
OS052_pf_final$logFC <- log2(OS052_pf_final$barcode_cpm_mean_pf / OS052_pf_final$barcode_cpm_mean_ctrl_13)


# Set significance level and fold change cutoffs
sig_level <- 0.05
fc_cutoff <- 1

# Create the volcano plot
OS052_PF_volcano <- ggplot(OS052_pf_final, aes(x=logFC, y=-log10(p_value))) +
  geom_point_rast(size=0.5, aes(color=ifelse(p_value< 0.05 & (logFC > 1 | logFC < -1), "red", "black")), show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  labs(title="OS052 CDK-4/6 inhibitor", x="logFC", y="-log10(p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8.5)) + 
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed", color="gray") +
  ylim(0, 30)


# Saving the PF plot
ggsave("~/Desktop/OS052_PF_volcano.svg", plot = OS052_PF_volcano, width = 2.2, height = 2.2, units = "in")



##



# Identifying the depleted and enriched barcodes
# filtering based on log fold change to identify positive fold change
depleted_filtered_atr <- OS052_atr_final %>% filter(logFC < -2.5 & p_value < 0.05)
enriched_filtered_atr <- OS052_atr_final %>% dplyr::filter(logFC > 1 & p_value < 0.05)
depleted_filtered_pf <- OS052_pf_final %>% filter(logFC < -1 & p_value < 0.05)



# making a vector of the top dropout barcodes for all the treatments
depleted_barcodes_atr <- depleted_filtered_atr$barcode
enriched_filtered_atr <- enriched_filtered_atr$barcode
depleted_barcodes_pf <- depleted_filtered_pf$barcode
depleted_barcodes_both <- c(depleted_barcodes_atr, depleted_barcodes_pf)


# changing the barcode type to DNA in order to later take the reverse complement
enriched_barcodes <- dna(enriched_filtered_atr)


# Getting the reverse complement of the top barcodes
rc_depleted_barcodes <- seq_complement(seq_reverse(depleted_barcodes))
rc_enriched_barcodes <- seq_complement(seq_reverse(enriched_barcodes))


# Creating a dataframe of the depleted barcodes
rc_depleted_barcodes <- as.data.frame(rc_depleted_barcodes)
rc_enriched_barcodes <- as.data.frame(rc_enriched_barcodes)


# Writing the dropout barcodes to the single cell analysis folder as txt file to be uploaded onto wynton
write.table(rc_enriched_barcodes$rc_enriched_barcodes, file = "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/enriched_LT_barcodes_atr_OS052_LT.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)



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



