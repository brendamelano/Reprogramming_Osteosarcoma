# reading in the libraries
library(tidyverse)
library(bioseq)
library(dplyr)


##### OS384   #########


#####   Depleted barcodes   ##########


# Reading in combined and trimmed fastq sequences
Fastq_sequences <- read.delim("~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/OS384_in_vivo_LT_barcodes.txt", header = F)


# reading in the depleted LT barcode sequences
depleted_barcodes <- read.csv("~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/depleted_barcodes_OS384_inVivo_LT.csv")
depleted_barcodes <- depleted_barcodes[,2]


# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(depleted_barcodes, collapse="|"), V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
cell_barcodes <- substr(filtered_fastqs$V1, 1, 16)


# Adding the "-1" to the barcode sequences
cell_barcodes <- paste(cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(cell_barcodes, "~/Desktop/depleted_cell_barcodesOS384_inVivo.csv")


## PD analysis of the depleted in Cluster 3

tenX_depleted_cluster_3 <- read.csv('~/Desktop/tenX_depleted_cluster_3.csv')
tenX_depleted_cluster_3_col <- tenX_depleted_cluster_3[,2]


# Filtering the fastq sequences based on the depleted barcode sequences
depleted_cluster_3_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(tenX_depleted_cluster_3_col, collapse="|"), V1))


# Function to find matching sequence
find_matching_barcode <- function(seq) {
  matching_barcodes <- sapply(depleted_barcodes, function(barcode) {
    if (grepl(barcode, seq)) {
      return(barcode)
    } else {
      return(NA)
    }
  })
  # Return the first matching barcode or NA if none found
  return(ifelse(any(!is.na(matching_barcodes)), matching_barcodes[!is.na(matching_barcodes)][1], NA))
}


# Apply the function to the dataframe
depleted_cluster_3_fastqs$LT_barcode <- sapply(depleted_cluster_3_fastqs$V1, find_matching_barcode)


# Save the LT_barcode column as a list without NAs
LT_depleted_cluster_3 <- as.vector(na.omit(depleted_cluster_3_fastqs$LT_barcode))


# Write the list out as a CSV
write.csv(LT_depleted_cluster_3, file="~/Desktop/LT_depleted_cluster_3.csv", row.names=FALSE)


#####   enriched barcodes   ##########


# reading in the depleted LT barcode sequences
enriched_barcodes <- read.csv("~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/enriched_barcodes_OS384_inVivo_LT.csv")
enriched_barcodes <- enriched_barcodes[,2]


# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(enriched_barcodes, collapse="|"), V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
enriched_cell_barcodes <- substr(filtered_fastqs$V1, 1, 16)


# adding the "-1" to the barcode sequences
enriched_cell_barcodes <- paste(enriched_cell_barcodes, "-1", sep = "")



# writing out the csv to upload into r on desktop
write.csv(enriched_cell_barcodes, "~/Desktop/enriched_cell_barcodesOS384_inVivo.csv")



#####   plotting barcodes for trajectory   ##########


# reading in the depleted LT barcode sequences
trajectory_barcodes <- read.csv("~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")
trajectory_barcodes <- trajectory_barcodes[,2]



# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs_trajectory <- dplyr::filter(Fastq_sequences, grepl(paste(trajectory_barcodes, collapse="|"), V1))
single_trajectory_fastq <- dplyr::filter(Fastq_sequences, grepl("GTTTTCATACATGCCATG", V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
trajectory_cell_barcodes <- substr(filtered_fastqs_trajectory$V1, 1, 16)
OS384_single_trajectory_cell_barcodes <- substr(single_trajectory_fastq$V1, 1, 16)


# adding the "-1" to the barcode sequences
trajectory_cell_barcodes <- paste(trajectory_cell_barcodes, "-1", sep = "")
OS384_single_trajectory_cell_barcodes <- paste(OS384_single_trajectory_cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(trajectory_cell_barcodes, "~/Desktop/scRNAseq_LT_analysis/trajectory_cell_barcodes.csv")
write.csv(OS384_single_trajectory_cell_barcodes, "~/Desktop/OS384_single_trajectory_cell_barcodes.csv")



##### Associating LT barcodes with 10X barcodes   #########


T0_384_barcodes <- read.csv("~/Desktop/OS384_time0_barcodes.csv")
T0_384_barcodes <- T0_384_barcodes[2]
names(T0_384_barcodes)[1] <- "LT_barcode"
LT_barcodes <- T0_384_barcodes$LT_barcode


df_sequences <- read.delim("~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/OS384_in_vivo_LT_barcodes.txt", header = F)
#LT_cell_barcodes <- data.frame(LT_cell_barcodes)


library(tidyverse)
tenX_cell_barcodes <- substr(LT_cell_barcodes$V1, 1, 16)
tenX_cell_barcodes <- data.frame(tenX_cell_barcodes)

# Function to find the LT barcode in a sequence
find_LT_barcode <- function(sequence, LT_barcodes) {
  LT_barcode_found <- NA  # Default is NA, changes if a barcode is found
  
  for (LT_barcode in LT_barcodes) {
    if (str_detect(sequence, LT_barcode)) {  # Check if the sequence contains the LT barcode
      LT_barcode_found <- LT_barcode
      break
    }
  }
  
  return(LT_barcode_found)
}


df_sequences <- df_sequences %>%
  rowwise() %>%
  mutate(LT_barcode = find_LT_barcode(V1, LT_barcodes))


final_df <- tenX_cell_barcodes %>%
  left_join(df_sequences, by = "V1") 



#######  OS742    #########




#####   Depleted barcodes   ##########


# reading in combined and trimmed fastq sequences
Fastq_sequences <- read.delim("~/Desktop/Osteo_Lineage_Tracing_Analysis/OS742_in_vivo_LT_barcodes.txt", header = F)


# reading in the depleted LT barcode sequences
depleted_barcodes <- read.csv("~/Desktop/Osteo_Lineage_Tracing_Analysis/15_day_treatments/depleted_barcodes_OS742_inVivo_LT.csv")
depleted_barcodes <- depleted_barcodes[,2]

# writing out the depleted barcodes into txt format
write(depleted_barcodes, file = "depleted_barcodes.txt", ncolumns = 1)


# The fastq sequences took forever to load, so I uploaded it to khayyam for OS742 and filtered on there
# grep -f OS742_depleted_barcodes.txt OS742_in_vivo_LT_barcodes.txt > OS742_depleted_FASTQ_10X_sequences.txt

filtered_fastqs <- read.delim("~/Desktop/Osteo_Lineage_Tracing_Analysis/scRNAseq_LT_analysis/OS742/OS742_depleted_FASTQ_10X_sequences.txt", header = F)

# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(depleted_barcodes, collapse="|"), V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
cell_barcodes <- substr(filtered_fastqs$V1, 1, 16)


# Adding the "-1" to the barcode sequences
cell_barcodes <- paste(cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(cell_barcodes, "~/Desktop/Osteo_Lineage_Tracing_Analysis/scRNAseq_LT_analysis/OS742/depleted_cell_barcodesOS742_inVivo.csv")


## PD analysis of the depleted in Cluster 3

tenX_depleted_cluster_3 <- read.csv('~/Desktop/tenX_depleted_cluster_3.csv')
tenX_depleted_cluster_3_col <- tenX_depleted_cluster_3[,2]


# Filtering the fastq sequences based on the depleted barcode sequences
depleted_cluster_3_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(tenX_depleted_cluster_3_col, collapse="|"), V1))


# Function to find matching sequence
find_matching_barcode <- function(seq) {
  matching_barcodes <- sapply(depleted_barcodes, function(barcode) {
    if (grepl(barcode, seq)) {
      return(barcode)
    } else {
      return(NA)
    }
  })
  # Return the first matching barcode or NA if none found
  return(ifelse(any(!is.na(matching_barcodes)), matching_barcodes[!is.na(matching_barcodes)][1], NA))
}


# Apply the function to the dataframe
depleted_cluster_3_fastqs$LT_barcode <- sapply(depleted_cluster_3_fastqs$V1, find_matching_barcode)


# Save the LT_barcode column as a list without NAs
LT_depleted_cluster_3 <- as.vector(na.omit(depleted_cluster_3_fastqs$LT_barcode))


# Write the list out as a CSV
write.csv(LT_depleted_cluster_3, file="~/Desktop/LT_depleted_cluster_3.csv", row.names=FALSE)


#####   enriched barcodes   ##########

# reading in combined and trimmed fastq sequences
Fastq_sequences <- read.delim("~/Desktop/Osteo_Lineage_Tracing_Analysis/scRNAseq_LT_analysis/OS742/OS742_in_vivo_LT_barcodes.txt", header = F)


# reading in the depleted LT barcode sequences
enriched_barcodes <- read.csv("~/Desktop/Osteo_Lineage_Tracing_Analysis/15_day_treatments/enriched_barcodes_OS742_inVivo_LT.csv")
enriched_barcodes <- enriched_barcodes[,2]

# writing out the depleted barcodes into txt format
write(enriched_barcodes, file = "~/Desktop/Osteo_Lineage_Tracing_Analysis/15_day_treatments/OS742_enriched_barcodes.txt", ncolumns = 1)


# The fastq sequences took forever to load, so I uploaded it to khayyam for OS742 and filtered on there
# grep -f OS742_enriched_barcodes.txt OS742_in_vivo_LT_barcodes.txt > OS742_enriched_FASTQ_10X_sequences.txt


filtered_fastqs <- read.delim("~/Desktop/Osteo_Lineage_Tracing_Analysis/scRNAseq_LT_analysis/OS742/OS742_enriched_FASTQ_10X_sequences.txt", header = F)


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
enriched_cell_barcodes <- substr(filtered_fastqs$V1, 1, 16)


# adding the "-1" to the barcode sequences
enriched_cell_barcodes <- paste(enriched_cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(enriched_cell_barcodes, "~/Desktop/Osteo_Lineage_Tracing_Analysis/scRNAseq_LT_analysis/OS742/enriched_cell_barcodesOS742_inVivo.csv")



#####   plotting barcodes for trajectory   ##########


# reading in the depleted LT barcode sequences
trajectory_barcodes <- read.csv("~/Desktop/scRNAseq_LT_analysis/OS384_trajectory_LT_barcodes.csv")
trajectory_barcodes <- trajectory_barcodes[,2]



# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs_trajectory <- dplyr::filter(Fastq_sequences, grepl(paste(trajectory_barcodes, collapse="|"), V1))
single_trajectory_fastq <- dplyr::filter(Fastq_sequences, grepl("GTTTTCATACATGCCATG", V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
trajectory_cell_barcodes <- substr(filtered_fastqs_trajectory$V1, 1, 16)
OS384_single_trajectory_cell_barcodes <- substr(single_trajectory_fastq$V1, 1, 16)


# adding the "-1" to the barcode sequences
trajectory_cell_barcodes <- paste(trajectory_cell_barcodes, "-1", sep = "")
OS384_single_trajectory_cell_barcodes <- paste(OS384_single_trajectory_cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(trajectory_cell_barcodes, "~/Desktop/scRNAseq_LT_analysis/trajectory_cell_barcodes.csv")
write.csv(OS384_single_trajectory_cell_barcodes, "~/Desktop/OS384_single_trajectory_cell_barcodes.csv")



##### Associating LT barcodes with 10X barcodes   #########


T0_384_barcodes <- read.csv("~/Desktop/OS384_time0_barcodes.csv")
T0_384_barcodes <- T0_384_barcodes[2]
names(T0_384_barcodes)[1] <- "LT_barcode"
LT_barcodes <- T0_384_barcodes$LT_barcode


df_sequences <- read.delim("~/Desktop/scRNAseq_LT_analysis/OS384_inVivo_scRNAseq_barcode_analysis/OS384_in_vivo_LT_barcodes.txt", header = F)
#LT_cell_barcodes <- data.frame(LT_cell_barcodes)


library(tidyverse)
tenX_cell_barcodes <- substr(LT_cell_barcodes$V1, 1, 16)
tenX_cell_barcodes <- data.frame(tenX_cell_barcodes)

# Function to find the LT barcode in a sequence
find_LT_barcode <- function(sequence, LT_barcodes) {
  LT_barcode_found <- NA  # Default is NA, changes if a barcode is found
  
  for (LT_barcode in LT_barcodes) {
    if (str_detect(sequence, LT_barcode)) {  # Check if the sequence contains the LT barcode
      LT_barcode_found <- LT_barcode
      break
    }
  }
  
  return(LT_barcode_found)
}


df_sequences <- df_sequences %>%
  rowwise() %>%
  mutate(LT_barcode = find_LT_barcode(V1, LT_barcodes))


final_df <- tenX_cell_barcodes %>%
  left_join(df_sequences, by = "V1") 



##### OS052   #########


#####   Depleted barcodes   ##########


# Reading in combined and trimmed fastq sequences
Fastq_sequences <- read.delim("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/OS052_LT_barcode_combined_sequences.txt", header = F)


# reading in the depleted LT barcode sequences
depleted_barcodes <- read.csv("~/Desktop/Reprogramming_Osteosarcoma/Lineage_Tracing/OS052/depleted_barcodes_OS052_LT.csv")
depleted_barcodes <- depleted_barcodes[,2]


# filtering the fastq sequences based on the depleted barcode sequences
filtered_fastqs <- dplyr::filter(Fastq_sequences, grepl(paste(depleted_barcodes, collapse="|"), V1))


# Getting the cell barcode sequences for the fastq sequences that remained (first 16 bases)
cell_barcodes <- substr(filtered_fastqs$V1, 1, 16)


# Adding the "-1" to the barcode sequences
cell_barcodes <- paste(cell_barcodes, "-1", sep = "")


# writing out the csv to upload into r on desktop
write.csv(cell_barcodes, "~/Desktop/depleted_cell_barcodesOS052_inVivo.csv")


