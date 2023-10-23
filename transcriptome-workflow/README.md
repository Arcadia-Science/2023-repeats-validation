# RNASeq Mapping & Counting Workflow

This Nextflow workflow takes in a samplesheet of pairs of genome accessions and SRA run accessions that originated from the same species as the sequenced genome and:
1. Downloads each SRA run accession
2. Downloads the genome FASTA and GTF annotation file for each Refseq accession
3. Builds a STAR index of the Refseq genome accession
4. Maps the RNAseq reads using STAR to the corresponding species genome
5. Performs gene counts with `htseq`
6. Creates a table of gene_ids to Refseq protein_ids from the GTF annotation file for downstream use

To use the workflow, you need both conda and nextflow installed.

1. Install conda [according to the instructions for your operating system](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html)
2. The easiest way to install Nextflow without worrying about dependency issues on your machine is through a conda environment, and can [install according to the instructions for your operation system](https://docs.conda.io/en/latest/miniconda.html). This is included in the `environment.yml` file. You can access the `environment.yml` file and all files neccessary for running the workflow with:

```
nextflow run main.nf --samples inputs/sra_refseq_download_table.csv \\
    --outdir results
    --threads <CPUS>
```

The main results are the htseq count tables and the parsed GTF files for each genome. The FASTQ files for each SRA run accession are not output into the results directory to save space on the machine.
