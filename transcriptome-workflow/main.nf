#! /usr/bin/env nextflow

// Description
// Download transcriptomic experiments and species genome/GTF files to map and quantify expression in select brain/muscle tissues from the SRA
// Assess the expression of repeat-expansion homologs identified in the species through protein sequence and structural comparison methods
// test with bengalese finch example and then process all ~3000 runs and ~145 species

// IMPORTANT: SRA-tools (fasterq-dump) must be installed locally in your path, the conda installation does not work

nextflow.enable.dsl=2

params.threads=10
params.outdir=null

log.info """\

DOWNLOAD, MAP, AND PROCESS TRANSCRIPTOMES AGAINST SPECIES GENOMES
TO QUANTIFY EXPRESSION OF SELECT GENES
=========================================
samples            : $params.samples
proteins           : $params.proteins
outdir             : $params.outdir
threads            : $params.threads
"""

// channels for the CSV of samples for SRA accessions and genome Refseq FTP links to download, mapping queries
// CSV to provide to select proteins from to quantify expression in the final table
// if SRA accession is single or paired end, split up into different output channels to process the command differently

// define channels and parameters
sample_csv = Channel.fromPath(params.samples)
    .splitCsv(header:true)

proteins_csv = Channel.fromPath(params.proteins)

// split SRA accessions based on if library_layout is SINGLE or PAIRED
paired_end_samples = sample_csv.filter { it.library_layout == "PAIRED" }.map { it.SRA_run_accession }
single_end_samples = sample_csv.filter { it.library_layout == "SINGLE" }.map { it.SRA_run_accession }

// channel linking each SRA run accesssion to the corresponding genome to map to
sra_genome_mapping_pairs = sample_csv.map { [ it.SRA_run_accession, it.genome_refseq_accession ] }

// channel linking genome refseq accession to FTP download path
refseq_genomes = sample_csv.map { [ it.genome_refseq_accession, it.genome_ftp_path ]}

workflow {
    downloaded_paired_end_reads = download_paired_SRA_runs(paired_end_samples)
    downloaded_single_end_reads = download_single_SRA_runs(single_end_samples)

    downloaded_refseq_genomes = download_refseq_files(refseq_genomes)

}
// download using SRA tools passing the SRA run accession
// paired-end reads process to split files
process download_paired_SRA_runs {
    // download each SRA run with SRAtools
    tag "${SRA_run_accession}_download"
    publishDir "${params.outdir}/sra_accessions", mode: 'copy', pattern:"*.fastq"

    input:
    val(SRA_run_accession)

    output:
    path("*.fastq"), emit: fastq

    script:
    """
    fastq-dump --gzip --split-3 ${SRA_run_accession}
    """
}

// download using SRA tools passing the SRA run accession
// single-end reads
process download_single_SRA_runs {
    // download each SRA run with SRAtools
    tag "${SRA_run_accession}_download"
    publishDir "${params.outdir}/sra_accessions", mode: 'copy', pattern:"*.fastq"

    input:
    val(SRA_run_accession)

    output:
    path("*.fastq"), emit: fastq

    script:
    """
    fastq-dump --gzip ${SRA_run_accession}
    """
}

// download the Refseq genome assembly and the GTF file by combining FTP link with Refseq accession
process download_refseq_files {
    tag "${genome_refseq_accession}_download"
    publishDir "${params.outdir}/refseq_assemblies", mode: 'copy', pattern: "*.gtf.gz"

    input:
    tuple val(genome_refseq_accession), val(genome_ftp_path)

    output:
    path("*.gtf.gz"), emit: gtf
    // path("*.fna.gz"), emit: fasta

    script:
    """
    wget ${genome_ftp_path}/${genome_refseq_accession}_genomic.gtf.gz
    """

}


// trim/QC reads trimmomatic


// map with HISAT2 either single-end with -U or paired-end with 1,2


// quantify with package


// select out proteins in the list for each genome
// python/R scripts

// binarize expression if 0=0, if greater than 0 then = 1 for just "expressed" in that tissue sample

// plot by species/tissue/gene binarized expression
//python/R scripts
