
# reading in sample csv
sample_csv <- read.csv('~/Downloads/crispr_example_feature_ref.csv')


## change IDs to the actual barcode sequence


# Adding the sequence after the LT barcode to recognize the LT barcode
sample_csv$pattern <- '(BC)CCTAGAGGGCCCGTTTAAAC'


# Adding the t0 barcodes to the csv
time_0_barcodes_df <- as.data.frame(time_0_barcodes)

# merging the barcodes with the sample csv
sample_csv_1 <- merge(time_0_barcodes_df, sample_csv)



names(sample_csv_1)[1] <- 'sequence'



sample_csv_2 <- sample_csv_1[1:961,]

sample_csv_3 <- sample_csv_2[,-6]

# got an error when this was the sequence, hopefully the LT barcode sequence is still represented in the data somehow
sample_csv_3$target_gene_id <- 'Non-Targeting'


sample_csv_3$target_gene_name <- sample_csv_3$sequence


# reassigning all values in the id column to LT
sample_csv_3$id <- paste0("LT_", seq(1:961))


# changing the gene name to LT
sample_csv_3$name <- paste0("LT_", seq(1:961))


# writing the OS384 cellranger feature reference csv to Desktop
write.csv(sample_csv_3, "~/Desktop/LT_whitelist_barcodes_742_inVivo.csv")





