library(tidyverse)
library(purrr)
library(readr)

srainfo_dir <- "results/srainfos/"
srainfo_files <- dir(srainfo_dir, pattern = "*-srainfo.tsv")

###########################################
# Read in all SRA info files per species, combine into a single table, filter for only Illumina RNA-Seq
# Basic filtering, tissue filtering in next step with species info
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
  select(species_name, run_accession, study_accession, experiment_accession, experiment_desc, organism_taxid, organism_name, library_strategy, library_source, library_selection, library_layout, sample_accession, biosample, bioproject, instrument, instrument_model_desc, total_size, run_total_bases, dev_stage, sex, tissue, collection_date, public_sratoolkit, public_url, ncbi_url)  %>%
  mutate(taxid = organism_taxid)  %>%
  select(-organism_taxid)  %>%
  mutate(species_clean = sub("-runinfo_[0-9]+", "", species_name))  %>%
  select(-species_name)  %>%
  filter(taxid != "9606") # remove human experiments that came down with some accessions like dogs

# filter for transcriptomic, illumina
srainfo_species_rna <- srainfo_species_clean %>%
  filter(library_strategy == "RNA-Seq") %>%
  filter(instrument_model_desc == 'ILLUMINA')  %>%
  filter(library_source == 'TRANSCRIPTOMIC')  %>%
  mutate(taxid = as.integer(taxid))  %>%
  filter(!is.na(run_total_bases))  %>%
  mutate(run_total_bases = as.numeric(run_total_bases))

write.csv(srainfo_species_rna, "results/srainfo_species_rna_accesions.csv", quote=FALSE, row.names = FALSE)

###########################################
# Metadata for genomes with Refseq accessions
# Combine with RNAseq experiments, filter further based on tissue
# Combine with data for FTP download links of Refseq accessions
# Download samplesheet is separate from the protein information samplesheet
###########################################

genome_refseq_ids <- read.csv("metadata/taxids_with_Refseq_accessions_for_expression_analysis.csv")  %>%
  mutate(taxid = Taxid)  %>%
  select(-Taxid)  %>%
  filter(taxid != "9606")

sra_genome_accessions <- left_join(srainfo_species_rna, genome_refseq_ids)  %>%
  filter(!is.na(refseq.accession))  %>%
  filter(!is.na(tissue))  %>%
  mutate(tissue_modf = tolower(tissue))  %>%
  filter(run_total_bases > 1000000)  %>%
  select(-tissue)

all_tissues <- sra_genome_accessions  %>% group_by(tissue_modf)  %>% count()  %>% arrange(desc(n))

# filter for brain/muscle/spinal cord experiments
sra_brain_muscle_accessions <- sra_genome_accessions  %>%
  filter(tissue_modf == "brain" | tissue_modf == "hippocampus" | tissue_modf == "spinal_cord" | tissue_modf == "hypothalamus" | tissue_modf == "skeletal muscle" | tissue_modf == "muscle" | tissue_modf == "cerebellum" | tissue_modf == "cortex" | tissue_modf == "dorsal root ganglion" | tissue_modf == "pituitary" | tissue_modf == "trigeminal ganglia" | tissue_modf == "prefrontal cortex" | tissue_modf == "dorsolateral prefrontal cortex" | tissue_modf == "brain tumor" | tissue_modf == "fetal brain" | tissue_modf == "midbrain" | tissue_modf == "brainstem" | tissue_modf == "brain (frontal cortex)" | tissue_modf == "cerebral cortex" | tissue_modf == "nucleus accumbens (brain)" | tissue_modf == "pituitary gland" | tissue_modf == "forebrain" | tissue_modf == "frontal lobe" | tissue_modf == "ganglia, spinal" | tissue_modf == "muscle, skeletal" | tissue_modf == "spinal cord" | tissue_modf == "brachial motor neurons" | tissue_modf == "brachial spinal cord")

sra_brain_muscle_accessions  %>%
  group_by(species_clean)  %>%
  count()  %>%
  arrange(desc(n))  %>%
  print(n=100)

brain_muscle_tissues <- sra_brain_muscle_accessions  %>%
  group_by(tissue_modf)  %>%
  count()  %>%
  arrange(desc(n))

# refseq assembly info from NCBI FTP site
refseq_assembly_info <- read_tsv("metadata/2023-10-05-assembly-summaries-refseq.tsv")  %>%
  select(refseq.accession, species_taxid, genome_size, scaffold_count, ftp_path)  %>%
  mutate(refseq_full_accession = str_extract(ftp_path, "(?<=\\/)[^\\/]+$"))

# join with the SRA table to get the FTP information
sra_brain_muscle_refseq_table <- left_join(sra_brain_muscle_accessions, refseq_assembly_info)

write.csv(sra_brain_muscle_refseq_table, "results/sra_brain_muscle_refseq_accessions.csv", quote=FALSE, row.names=FALSE)


###########################################
# Per RefSeq genome, RefSeq protein IDs identified as repeat-expansion homologs
###########################################

protein_results <- read.csv("metadata/2023-08-02-repeat-expansion-full-results.csv")  %>%
  mutate(taxid = Taxid)  %>%
  select(-Taxid) # remove human homologs so 1:1 mapping of homologs directly from the species

# pull the taxids that met the RNAseq criteria
taxid_list <- sra_brain_muscle_accessions  %>%
  group_by(taxid)  %>%
  filter(taxid != "9606")  %>%
  pull(taxid)  %>%
  unique()

# table of protein accessions for the corresponding species that have brain/muscle RNASeq experiments
# and also have a Refseq genome to pull the genome and GTF file from
proteins_brain_muscle_species_table <- protein_results  %>%
  filter(taxid %in% taxid_list)  %>%
  select(gene, Accession, Scientific.Name, taxid)  %>%
  left_join(genome_refseq_ids)

proteins_brain_muscle_species_table  %>%
  group_by(taxid)  %>%
  count()  %>%
  arrange(desc(n))  %>%
  print(n=50)

###########################################
# Output files for use in workflow
# First file is mapping file for each SRA accession to map to which Refseq genome, download locations for the SRA accession, and the Refseq FTP link for the genome and GTF file to create
# Second file is protein information for each genome, to extract from featurecounts
###########################################

sra_refseq_download_table <- sra_brain_muscle_refseq_table  %>%
  select(species_clean, run_accession, library_layout, tissue_modf, refseq_full_accession, ftp_path)  %>%
  mutate(species_name = species_clean)  %>%
  mutate(tissue = tissue_modf)  %>%
  mutate(genome_refseq_accession = refseq_full_accession)  %>%
  mutate(genome_ftp_path = ftp_path)  %>%
  mutate(SRA_run_accession = run_accession)  %>%
  select(species_name, SRA_run_accession, library_layout, tissue, genome_refseq_accession, genome_ftp_path)

refseq_proteins_table <- proteins_brain_muscle_species_table  %>%
  select(refseq.accession, gene, Accession)  %>%
  mutate(protein_refseq_accession = Accession)  %>%
  mutate(genome_refseq_accession = refseq.accession)  %>%
  mutate(protein_name = gene)  %>%
  select(genome_refseq_accession, protein_name, protein_refseq_accession)

# write these samplesheets to the workflow inputs folder to ingest as part of the workflow
write.csv(sra_refseq_download_table, "transcriptome-workflow/inputs/sra_refseq_download_table.csv", row.names = FALSE, quote = FALSE)

write.csv(refseq_proteins_table, "transcriptome-workflow/inputs/refseq_proteins_table.csv", row.names = FALSE, quote = FALSE)
