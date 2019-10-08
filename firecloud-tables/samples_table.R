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
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bams >> /merged_bam_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bais >> /merged_bai_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_cellbarcodes >> /cellbarcodes.txt")
bam <- fread("firecloud-tables/merged_bam_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
bai <- fread("firecloud-tables/merged_bai_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
cell_barcodes <- fread("firecloud-tables/cellbarcodes.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")


# BCC and SCC refer to basal cell carcinoma and squamos cell carcinoma
bam <- bam %>% rename(bam = colnames(bam)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(bam)), 
         sample = gsub("\\..*" , "", basename(bam))) %>%
  select(subject, sample, bam)

bai <- bai %>% rename(bamIndex = colnames(bai)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(bamIndex)), 
         sample = gsub("\\..*" , "", basename(bamIndex))) %>%
  select( sample, bamIndex)

cell_barcodes <- cell_barcodes %>% rename(cellBarcode = colnames(cell_barcodes)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(cellBarcode)), 
         sample = gsub("\\..*" , "", basename(cellBarcode)),
         sample = gsub("_RNA_barcodes","",sample)) %>%
  select( sample, cellBarcode)

subjects <- unique(normal.bam$subject) ## filter based on existing contorl samples

##################
## samples file ##
##################
case_samples_file <- left_join(bam, bai, by = "sample") %>%
  left_join(., cell_barcodes, by = "sample") %>%
  add_column(vcf = "",
             vcfIndex = "",
             fasta = "",
             fastaIndex = "",
             output = "",
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
         fasta,
         fastaIndex,
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




