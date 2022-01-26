#!/bin/bash

# This script download the BAM files from TCGA dataset
# and extract the specific region of interest (ChrX:83000000-85000000)

# loading required packages
library(GenomicDataCommons)

# Token authentication 
  # GDC token txt file is required as input 
  # getwd()
setwd("/Users/chensisi/Documents/")
tokenfile <- "gdc-user-token.20220104T15_32_46.473Z.txt"
token <- function(tokenfile){
  tokentxt <- paste(readLines(tokenfile), collapse=" ")
  Sys.setenv(GDC_TOKEN = tokentxt)
  token = gdc_token()
}


### Download BAM files of all cancer types from TCGA 

# Create directory for each cancer type 
TCGA_list <- c("TCGA-BRCA","TCGA-GBM","TCGA-OV","TCGA-LUAD","TCGA-UCEC","TCGA-KIRC","TCGA-HNSC",
               "TCGA-LGG","TCGA-THCA","TCGA-LUSC","TCGA-PRAD","TCGA-SKCM","TCGA-COAD","TCGA-STAD",
               "TCGA-BLCA","TCGA-LIHC","TCGA-CESC","TCGA-KIRP","TCGA-SARC","TCGA-LAML","TCGA-ESCA",
               "TCGA-PAAD","TCGA-PCPG","TCGA-READ","TCGA-TGCT","TCGA-THYM","TCGA-KICH","TCGA-ACC",
               "TCGA-MESO","TCGA-UVM","TCGA-DLBC","TCGA-UCS","TCGA-CHOL")
TCGA_list <- as.data.frame(TCGA_list)

setwd("/Users/chensisi/Documents/test")
for (i in 1:nrow(TCGA_list)){
  CancerType <- TCGA_list[i,]
  dir.create(CancerType)
  print(paste0('Directory for ', CancerType, ' has been created!'))
}

# Download BAM files for each cancer type 
summary = data.frame()

for (i in 1:nrow(TCGA_list)){
  CancerType <- TCGA_list[i,]
  
  # Fetch TCGA file IDs for each cancer type 
  TCGA <- files() %>%
    filter(~ cases.project.project_id == CancerType &
             data_type == 'Aligned Reads' &
             experimental_strategy == 'RNA-Seq' &
             data_format == 'BAM')
  TCGA_sample_ids = TCGA %>% response_all() %>% ids()
  
  # Note down the file IDs and number of files in each cancer type 
  num = length(TCGA_sample_ids)
  print(paste0('Number of BAM files for ', CancerType , ' is ', num))
  samples <- data.frame(TCGA_sample_ids, CancerType)
  summary <- rbind(summary, samples)
  
  # File downloading 
  for (x in 1:num){
    dir <- file.path("/Users/chensisi/Documents/test", CancerType)
    bamfile = slicing(
      TCGA_sample_ids[x], regions="chrX:83000000-85000000", 
      destination = file.path(dir,paste0(TCGA_sample_ids[x],".bam")))
    print(paste0(x, " out of ", num, " BAM files for ", CancerType," has been downloaded!"))
  }
}

# Export a file with all sample IDs for each cancer type 
write.table(summary, "sampleIDs.txt",sep = ",",row.names = FALSE, col.names = TRUE)


