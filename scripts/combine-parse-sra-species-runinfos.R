library(dplyr)
library(purrr)
library(readr)

srainfo_dir <- "results/srainfos/"
srainfo_files <- dir(srainfo_dir, pattern = "*-srainfo.tsv")

###########################################
# Read in all SRA info files per species, combine into a single table, filter for only Illumina RNA-Seq
###########################################

# Custom function to read a file
read_file <- function(file_path) {
  data <- read_tsv(file_path, col_types = cols(.default = "c")) # Read all columns as character for conflicting character/date cells
  return(data)
}

# Read and preprocess data
srainfo_species <- data_frame(filename = srainfo_files) %>% 
  mutate(file_contents = map(filename, ~ read_file(file.path(srainfo_dir, .)))) %>% 
  unnest()

# Convert collection_date to date type
srainfo_species_clean <- srainfo_species %>%
  mutate(collection_date = as.Date(collection_date)) %>% 
  mutate(species_name = gsub("-srainfo.tsv", "", filename)) %>% 
  select(species_name, run_accession, study_accession, experiment_accession, experiment_desc, organism_taxid, organism_name, library_strategy, library_source, library_selection, library_layout, sample_accession, biosample, bioproject, instrument, instrument_model_desc, run_total_bases, dev_stage, sex, tissue, collection_date)

# filter for transcriptomic, illumina
srainfo_species_rna <- srainfo_species_clean %>% 
  filter(library_strategy == "RNA-Seq") %>% 
  filter(instrument_model_desc == 'ILLUMINA')

write.csv(srainfo_species_rna, "results/srainfo_spcies_rna_accesions.csv", quote=FALSE, row.names = FALSE)

###########################################
# Metadata for genomes with gene models
###########################################

genbank_genome_metadata <- read_tsv("metadata/20230227_genbank_genomes_cds.tsv") %>% 
  filter(genbank_group == "vertebrate_mammalian" | genbank_group == "vertebrate_other") %>% 
  filter(organism_name != "NA")

