library(tidyverse)
library(viridisLite)
library(ggpubr)
library(stringr)
library(stringi)

###########################################
# Parse counts results and inspect homologs of interest
###########################################

###########################################
# Prep counts tables
###########################################

# paths
counts_dir <- "transcriptome-workflow/results/v2/counts/"
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
 
# outlier removal prior to density plotting to look at distribution of counts
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
  mutate(species_name = gsub("-", " ", Common.Name)) %>% 
  select(-Accession, -Common.Name)

# path for parsed gtf tables
gtf_dir <- "transcriptome-workflow/results/v2/gtf_tables/"
gtf_files <- dir(gtf_dir, pattern = "*.tsv")

# Custom function to read tsvs
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

# sample counts for each species and tissue type
sample_counts <- counts_table_info %>% 
  group_by(species_name, tissue) %>% 
  summarize(n_total_samples = n_distinct(SRA_run_accession))

###########################################
# Plotting
###########################################

# Arcadia color scheme
black <- "#09090A"
grape <- "#5A4596"
taffy <- "#E87485"
tangerine <- "#FFB984"
oat <- "#F5E4BE"

magma_colors <- c(black, grape, taffy, tangerine, oat)

magma_gradient <- colorRampPalette(magma_colors)

gradient_100 <- magma_gradient(100)

species_percent_samples_expression <- homolog_count_table_stats %>% 
  select(species_name, gene, binary_count, tissue) %>% 
  left_join(sample_counts, by = c('species_name', 'tissue')) %>% 
  mutate(species_upper = stri_trans_totitle(species_name, opts_brkiter = stri_opts_brkiter(type="sentence"))) %>% 
  mutate(tissue_upper = stri_trans_totitle(tissue, opts_brkiter = stri_opts_brkiter(type="sentence"))) %>% 
  group_by(species_upper, gene, tissue_upper) %>% 
  mutate(n_expressed_samples = sum(binary_count)) %>% 
  mutate(percent_expressed = n_expressed_samples / n_total_samples) %>% 
  ggplot(aes(x=tissue_upper, y=gene, fill=percent_expressed)) +
  geom_tile(color="white", linewidth=0.5) +
  facet_wrap(~species_upper, scales="free", nrow=2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust=1), legend.position = "bottom", axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_gradientn(colors=magma_gradient(100), name = "Percent Expressed")

species_percent_samples_expression

# write table to csv
species_expression_table <- species_percent_samples_expression <- homolog_count_table_stats %>% 
  left_join(sample_counts, by = c('species_name', 'tissue')) %>% 
  group_by(species_name, gene, tissue) %>% 
  mutate(n_expressed_samples = sum(binary_count)) %>% 
  mutate(percent_expressed = n_expressed_samples / n_total_samples)

write.csv(species_expression_table, "results/species-expression-counts-stats.csv", row.names = FALSE, quote = FALSE)

# save plot
ggsave("figs/all-species-tissue-expression-plots.png", species_percent_samples_expression, width=11, height=8, units=c("in"))
ggsave("figs/all-species-tissue-expression-plots.jpg", species_percent_samples_expression, width=11, height=8, units=c("in"))
ggsave("figs/all-species-tissue-expression-plots.pdf", species_percent_samples_expression, width=11, height=8, units=c("in"))
