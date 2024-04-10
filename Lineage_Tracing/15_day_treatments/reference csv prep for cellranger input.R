


# reading in sample csv
sample_csv <- read.csv('~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/LT_whitelist_barcodes_sample.csv')


# Taking the reverse complement of the white list barcodes
OS052_time_0_barcodes_rc <- dna(OS052_time_0_barcodes)

# Taking the reverse compliment of the time 0 barcodes
OS052_time_0_barcodes_rc <- seq_complement(seq_reverse(OS052_time_0_barcodes_rc))

# Adding the t0 barcodes to the csv
OS052_time_0_barcodes_df <- as.data.frame(OS052_time_0_barcodes_rc)

# Removing 200 rows
sample_csv <- sample_csv[-(1:200),]


# merging the barcodes with the sample csv
sample_csv_1 <- merge(OS052_time_0_barcodes_df, sample_csv)


# Removing the 2nd column to add a new sequence column
sample_csv_1 <- sample_csv_1[,-2]


names(sample_csv_1)[1] <- 'sequence'



sample_csv_2 <- sample_csv_1[1:900,]

sample_csv_3 <- sample_csv_2[,-7]


sample_csv_3$target_gene_name <- sample_csv_3$sequence

new_ids <- paste0("LT_", 1:900)

# Replace the 'id' and 'name' columns with the new sequence
sample_csv_3$id <- new_ids
sample_csv_3$name <- new_ids

# writing the OS384 cellranger feature reference csv to Desktop
write.csv(sample_csv_3, "~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/scRNAseq_LT_analysis/LT_whitelist_barcodes_052_inVivo_rc.csv")



