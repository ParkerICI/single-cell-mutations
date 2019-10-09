## TCRseq analysis
rm(list = ls())
library(data.table)
library(tidyverse)
library(dplyr)
library(magrittr)
setwd("~/vartrix-project")

## columns bam, bai, cell barcodes, vcf, sample name, output 
# for i in *.bam ; do echo $i ; samtools index $i $i".bai" ; done


## read gslinks in bucket folders and generate into files
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bams >> firecloud-tables/merged_bam_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bais >> firecloud-tables/merged_bai_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_cellbarcodes >> firecloud-tables/cellbarcodes.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged-vcfs >> firecloud-tables/vcf.txt")

## read-in files containing gs links
bam <- fread("firecloud-tables/merged_bam_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
dim(bam)
bai <- fread("firecloud-tables/merged_bai_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
dim(bai)
cell_barcodes <- fread("firecloud-tables/cellbarcodes.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
dim(cell_barcodes)
vcfFiles <- fread("firecloud-tables/vcf.txt", 
                       stringsAsFactors = FALSE,
                       check.names = FALSE,
                       sep = "\t")
dim(vcfFiles)
### prepare files
# BCC and SCC refer to basal cell carcinoma and squamos cell carcinoma
bam <- bam %>% rename(bam = colnames(.)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(bam)), 
         sample = gsub("\\..*" , "", basename(bam))) %>%
  select(subject, sample, bam)

bai <- bai %>% rename(bamIndex = colnames(.)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(bamIndex)), 
         sample = gsub("\\..*" , "", basename(bamIndex))) %>%
  select( sample, bamIndex)

cell_barcodes <- cell_barcodes %>% rename(cellBarcode = colnames(.)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(cellBarcode)), 
         sample = gsub("\\..*" , "", basename(cellBarcode)),
         sample = gsub("_RNA_barcodes","",sample)) %>%
  select( sample, cellBarcode)

vcf <- vcfFiles %>% rename(vcf = colnames(.)[1]) %>% 
  filter(!grepl(".gz.tbi", vcf)) %>%
  mutate(subject = paste("su00", gsub("_.*", "", basename(vcf)), sep = "")) %>%
  select(subject, vcf)

vcfIndex <- vcfFiles %>% rename(vcfIndex = colnames(.)[1]) %>%
  filter(grepl(".gz.tbi", vcfIndex)) %>%
  mutate(subject = paste("su00", gsub("_.*", "", basename(vcfIndex)), sep = "")) %>%
  select(subject, vcfIndex)

subjects <- unique(normal.bam$subject) ## filter based on existing contorl samples

##################
## samples file ##
##################
case_samples_file <- left_join(bam, bai, by = "sample") %>%
  left_join(., cell_barcodes, by = "sample") %>%
  left_join(., vcf, by = "subject") %>%
  left_join(., vcfIndex, by = "subject") %>%
  add_column(output = "",
             snv_loci = "")

samples_file <- case_samples_file %>%
  drop_na(bam, bamIndex, cellBarcode) %>% filter(subject %in% subjects) %>%
  select(`entity:sample_id` = sample,
         participant = subject,
         bam,
         bamIndex,
         cellBarcode,
         vcf,
         vcfIndex,
         output,
         snv_loci)


write.table(samples_file,
            "firecloud-tables/samples.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

#######################
## participants file ##
#######################
participants_file <- select(samples_file, `entity:participant_id` = participant)
participants_file <- unique(participants_file) %>% filter(`entity:participant_id` %in%  subjects)

write.table(participants_file,
            "firecloud-tables/participants.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")




