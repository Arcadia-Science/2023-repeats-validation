library(tidyverse)
library(purrr)
library(readr)

counts_dir <- "transcriptome-workflow/results/test/counts/"
files <- dir(counts_dir, pattern = "*.htseq")

###########################################
# Read in all SRA info files per species, combine into a single table, filter for only Illumina RNA-Seq
# Basic filtering, tissue filtering in next step with species info
###########################################

# Custom function to read a file
read_file <- function(file_path) {
  data <- read_tsv(file_path, col_names = c("gene", "count"))
  return(data)
}

# Read and preprocess data
all_counts_table <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read_file(file.path(counts_dir, .)))) %>%
  unnest() %>% 
  mutate(filename = gsub(".htseq", "", filename))

all_counts_table <- all_counts_table %>% 
  separate(filename, into = c("genome_refseq_accession", "SRA_run_accession"), sep = "_vs_")

# metadata for SRA accessions, protein list with homologs of interest
sra_accessions <- read.csv("transcriptome-workflow/inputs/sra_refseq_download_table.csv")
protein_accessions <- read.csv("transcriptome-workflow/inputs/refseq_proteins_table.csv")

# join with the count information to add tissue metadata, filter for the proteins of interest
counts_table_info <- left_join(all_counts_table, sra_accessions) %>% 
  select(genome_refseq_accession, SRA_run_accession, species_name, gene, count, tissue)

# outlier removal prior to density plotting
counts_table_info %>% 
  group_by(SRA_run_accession) %>% 
  mutate(IQR = IQR(count),
         Q1 = quantile(count, 0.25),
         Q3 = quantile(count, 0.75)) %>% 
  filter(count >= (Q1 - 1.5 * IQR) & count <= (Q3 + 1.5 * IQR)) %>%
  filter(species_name == "Bengalese_finch") %>% 
  ggplot(aes(x=count)) +
  geom_density(adjust = 2)
