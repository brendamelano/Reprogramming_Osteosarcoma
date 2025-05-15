
##############


### Preparing Rds files for mave - need to move to another script

library(Matrix)

# Read the count matrix
count_matrix <- readMM('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS152_count_matrix.mtx')

# Read cell barcodes and gene names
cell_barcodes <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS152_barcodes.tsv')
gene_names <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS152genes.tsv')

# Assign row and column names
rownames(count_matrix) <- cell_barcodes
colnames(count_matrix) <- gene_names

# Save as .rds file
saveRDS(count_matrix, '/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/OS152_gene_barcode_matrix.rds')

# Check the dimensions
dim(count_matrix)

# View the first few rows and columns
count_matrix[1:5, 1:5]

# Verify row and column names
head(rownames(count_matrix))
head(colnames(count_matrix))


## 384 ##



# Read the count matrix
count_matrix <- readMM('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS384_count_matrix.mtx')

# Read cell barcodes and gene names
cell_barcodes <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS384_barcodes.tsv')
gene_names <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS384_genes.tsv')

# Assign row and column names
rownames(count_matrix) <- cell_barcodes
colnames(count_matrix) <- gene_names

# Save as .rds file
saveRDS(count_matrix, '/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/OS384_gene_barcode_matrix.rds')

# Check the dimensions
dim(count_matrix)

# View the first few rows and columns
count_matrix[1:5, 1:5]

# Verify row and column names
head(rownames(count_matrix))
head(colnames(count_matrix))



## 742 ##



# Read the count matrix
count_matrix <- readMM('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS742_count_matrix.mtx')

# Read cell barcodes and gene names
cell_barcodes <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS742_barcodes.tsv')
gene_names <- readLines('/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/files_for_rds_generation/OS742_genes.tsv')

# Assign row and column names
rownames(count_matrix) <- cell_barcodes
colnames(count_matrix) <- gene_names

# Save as .rds file
saveRDS(count_matrix, '/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/OS742_gene_barcode_matrix.rds')

# Check the dimensions
dim(count_matrix)

# View the first few rows and columns
count_matrix[1:5, 1:5]

# Verify row and column names
head(rownames(count_matrix))
head(colnames(count_matrix))
