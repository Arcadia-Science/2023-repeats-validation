library(tidyverse)

###########################################
# Parse counts results and inspect homologs of interest
###########################################

###########################################
# Prep counts tables
###########################################

# paths
counts_dir <- "transcriptome-workflow/results/counts/"
files <- dir(counts_dir, pattern = "*.htseq")

# Custom function to read a file
read_file <- function(file_path) {
  data <- read_tsv(file_path, col_names = c("gene", "count"))
  data <- data %>% filter(!str_detect(gene, "^__"))
  return(data)
}

# Read and preprocess data
all_counts_table <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read_file(file.path(counts_dir, .)))) %>%
  unnest() %>% 
  mutate(filename = gsub(".htseq", "", filename))

all_counts_table <- all_counts_table %>% 
  separate(filename, into = c("genome_refseq_accession", "SRA_run_accession"), sep = "_vs_") %>% 
  filter(genome_refseq_accession != "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b") # remove chicken accessions that got mixed up with little skate

# metadata for SRA accessions
sra_accessions <- read.csv("transcriptome-workflow/inputs/sra_refseq_download_table.csv")

# join with the count information to add tissue metadata, filter for the proteins of interest
counts_table_info <- left_join(all_counts_table, sra_accessions) %>% 
  select(genome_refseq_accession, SRA_run_accession, species_name, gene, count, tissue) %>% 
  mutate(species_name = gsub("_", " ", species_name)) %>% 
  mutate(species_name = gsub("-", " ", species_name))
 
# outlier removal prior to density plotting
counts_table_info %>% 
  group_by(SRA_run_accession) %>% 
  mutate(IQR = IQR(count),
         Q1 = quantile(count, 0.25),
         Q3 = quantile(count, 0.75)) %>% 
  filter(count >= (Q1 - 1.5 * IQR) & count <= (Q3 + 1.5 * IQR)) %>%
  ggplot(aes(x=count)) +
  geom_density(adjust=1/5) + 
  facet_grid(~ species_name, scales = "free")

counts_table_info %>% 
  group_by(SRA_run_accession) %>% 
  mutate(IQR = IQR(count),
         Q1 = quantile(count, 0.25),
         Q3 = quantile(count, 0.75)) %>% 
  filter(count >= (Q1 - 1.5 * IQR) & count <= (Q3 + 1.5 * IQR)) %>%
  filter(species_name == "little_skate") %>% 
  ggplot(aes(x=count, fill=SRA_run_accession)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_density(adjust = 100, position = "stack", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "none")

# sum of counts for each SRA accession
counts_table_info %>% 
  group_by(SRA_run_accession, species_name) %>% 
  summarize(total_count = sum(count)) %>% 
  arrange(desc(total_count)) %>% 
  print(n=140)

# summary stats
run_stats <- counts_table_info %>% 
  group_by(SRA_run_accession) %>% 
  mutate(count_n = as.numeric(count)) %>% 
  summarise(
    Q1 = signif(1.5 * quantile(count, 0.25, na.rm = TRUE), digits=3),
    Median = signif(median(count, na.rm = TRUE), digits=3),
    Q3 = signif(1.5 * quantile(count, 0.75, na.rm = TRUE), digits=3)
  )

###########################################
# Protein accession information
###########################################

# filter for specific gene/protein accessions for homologs of interest
protein_accessions <- read.csv("metadata/2023-08-02-repeat-expansion-profiles.csv") %>% 
  select(Common.Name, Accession, gene) %>% 
  mutate(refseq_id = Accession) %>% 
  mutate(species_name = Common.Name) %>% 
  select(-Accession, -Common.Name)

# path for parsed gtf tables
gtf_dir <- "transcriptome-workflow/results/gtf_tables/"
gtf_files <- dir(gtf_dir, pattern = "*.tsv")

# Custom function to read a file
read_gtf <- function(file_path) {
  data <- read_tsv(file_path, col_names = c("gene", "refseq_id"))
  return(data)
}

# Read and preprocess data
all_gtf_tables <- data_frame(filename = gtf_files) %>%
  mutate(file_contents = map(filename, ~ read_gtf(file.path(gtf_dir, .)))) %>%
  unnest() %>% 
  mutate(filename = gsub("_gtf_table.tsv", "", filename))

# join gtf tables with proteins accessions of interest
hit_accessions <- left_join(all_gtf_tables, protein_accessions, by="refseq_id") %>% 
  drop_na() %>% 
  mutate(gene = gene.x) %>% 
  select(-gene.x, -gene.y)

###########################################
# Join protein accessions with counts
###########################################

# join accessions of homologs of interest with count table 
homolog_count_table <- left_join(counts_table_info, hit_accessions, by=c("species_name", "gene")) %>% 
  drop_na() %>% 
  select(species_name, SRA_run_accession, tissue, gene, count, refseq_id) %>% 
  group_by(SRA_run_accession, gene) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

# join homolog count table with run summary stats, add binary count compared to median
homolog_count_table_stats <- left_join(homolog_count_table, run_stats) %>% 
  mutate(binary_count = ifelse(count < Median, 0, 1))

###########################################
# Plotting
###########################################
