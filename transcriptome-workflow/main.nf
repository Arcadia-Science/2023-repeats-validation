#! /usr/bin/env nextflow

// Description
// Download transcriptomic experiments and species genome/GTF files to map and quantify expression in select brain/muscle tissues from the SRA
// Assess the expression of repeat-expansion homologs identified in the species through protein sequence and structural comparison methods
// test with bengalese finch example and then process all ~3000 runs and ~145 species

nextflow.enable.dsl=2

params.threads=10
params.outdir=null

log.info """\
DOWNLOAD, MAP, AND PROCESS TRANSCRIPTOMES AGAINST SPECIES TO QUANTIFY EXPRESSION OF SELECT GENES
=========================================
samples            : $params.samples
proteins           : $params.proteins
outdir             : $params.outdir
"""

// channels for the CSV of samples for SRA accessions and genome Refseq FTP links to download, mapping queries

// CSV to provide to select proteins from to quantify expression in the final table

// if SRA accession is single or paired end, split up into different output channels to process the command differently
// download using SRA tools passing the SRA run accession

// download the Refseq genome assembly and the GTF file by grabbing the FTP link

// trim/QC reads

// map with HISAT2 either single-end with -U or paired-end with 1,2

// quantify with package

// select out proteins in the list for each genome

// binarize expression if 0=0, if greater than 0 then = 1 for just "expressed" in that tissue sample

// plot by species/tissue/gene binarized expression
