## TCRseq analysis
rm(list = ls())
library(data.table)
library(tidyverse)
library(dplyr)
library(magrittr)
setwd("~/vartrix-project/data")

## columns bam, bai, cell barcodes, vcf, sample name, output 
# for i in *.bam ; do echo $i ; samtools index $i $i".bai" ; done


## read gslinks in bucket folders and generate into files
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bams >> /merged_bam_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_bais >> /merged_bai_gslinks.txt")
system("gsutil ls gs://stanford-sc-rnaseq-bcc/merged_cellbarcodes >> /cellbarcodes.txt")
bam <- fread("data/ansu-scRNA/processing/merged_bam_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
bai <- fread("data/ansu-scRNA/processing/merged_bai_gslinks.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")
cell_barcodes <- fread("data/ansu-scRNA/processing/cellbarcodes.txt", 
             stringsAsFactors = FALSE,
             check.names = FALSE,
             sep = "\t")


# BCC and SCC refer to basal cell carcinoma and squamos cell carcinoma
bam <- bam %>% rename(star_bam_file = colnames(bam)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(star_bam_file)), 
         sample = gsub("\\..*" , "", basename(star_bam_file))) %>%
  select(subject, sample, star_bam_file)

bai <- bai %>% rename(star_bai_file = colnames(bai)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(star_bai_file)), 
         sample = gsub("\\..*" , "", basename(star_bai_file))) %>%
  select( sample, star_bai_file)

cell_barcodes <- cell_barcodes %>% rename(cell.barcode = colnames(cell_barcodes)[1]) %>% 
  mutate(subject = gsub("_.*", "", basename(cell.barcode)), 
         sample = gsub("\\..*" , "", basename(cell.barcode)),
         sample = gsub("_RNA_barcodes","",sample)) %>%
  select( sample, cell.barcode)

subjects <- unique(normal.bam$subject) ## filter based on existing contorl samples

##################
## samples file ##
##################
case_samples_file <- left_join(bam, bai, by = "sample") %>%
  left_join(., cell_barcodes, by = "sample") %>%
  add_column(vcf = "",
             vcf.i = "",
             fasta = "",
             fasta.i = "",
             output = "",
             snv.loci = "")

samples_file <- case_samples_file %>%
  drop_na(star_bam_file, star_bai_file, cell.barcode) %>% filter(subject %in% subjects) %>%
  select(`entity:sample_id` = sample,
         participant = subject,
         star_bam_file,
         star_bai_file,
         cell.barcode,
         vcf,
         vcf.i,
         # fasta,
         # fasta.i,
         output,
         snv.loci)


write.table(samples_file,
            "data/ansu-scRNA/processing/samples.txt",
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
            "data/ansu-scRNA/processing/participants.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")




