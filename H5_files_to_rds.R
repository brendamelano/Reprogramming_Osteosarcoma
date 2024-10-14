# Load the package
library(hdf5r)
library(rhdf5)


# Alternatively, using rhdf5
data <- h5read("/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/h5_files/152_adata_PCA_subtypes.h5", "/X")
data <- h5read("/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/h5_files/384_adata_PCA_subtypes.h5", "/X")
data <- h5read("/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/h5_files/742_adata_PCA_subtypes.h5", "/X")

# Save the data as an .rds file
saveRDS(data, "/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/152_adata_PCA_subtypes.rds")
saveRDS(data, "/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/384_adata_PCA_subtypes.rds")
saveRDS(data, "/Users/brendamelano/Desktop/Reprogramming_Osteosarcoma/plain_scRNAseq_analysis/Brenda_data_for_Mehran/742_adata_PCA_subtypes.rds")
