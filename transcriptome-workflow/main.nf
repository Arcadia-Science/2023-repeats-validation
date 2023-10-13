#! /usr/bin/env nextflow

// Description
// Download transcriptomic experiments and species genome/GTF files to map and quantify expression in select brain/muscle tissues from the SRA
// Assess the expression of repeat-expansion homologs identified in the species through protein sequence and structural comparison methods
// test with bengalese finch example and then process all ~3000 runs and ~145 species

// IMPORTANT: SRA-tools (fasterq-dump) must be installed locally in your path, the conda installation does not work

nextflow.enable.dsl=2

params.threads=6
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
// add the corresponding refseq_accession so carries through as a tuple
paired_end_samples = sample_csv.filter { it.library_layout == "PAIRED" }.map { it.SRA_run_accession }
single_end_samples = sample_csv.filter { it.library_layout == "SINGLE" }.map { it.SRA_run_accession }

// channel linking each SRA run accesssion to the corresponding genome to map to, this will be combined later
sra_genome_mapping_pairs = sample_csv.map { [ it.SRA_run_accession, it.genome_refseq_accession ] }

// channel linking genome refseq accession to FTP download path
refseq_genomes = sample_csv.map { [ it.genome_refseq_accession, it.genome_ftp_path ]}
    .unique() //deduplicate because samplesheet is set up for mapping pairs

workflow {
    // download SRA files
    downloaded_paired_end_reads = download_paired_SRA_runs(paired_end_samples)
    downloaded_single_end_reads = download_single_SRA_runs(single_end_samples)

    // download Refseq files
    download_refseq_files(refseq_genomes)
    downloaded_fasta = download_refseq_files.out.fasta
    downloaded_gtf = download_refseq_files.out.gtf

    // build STAR index, include the GTF and FASTA for exons
    // carry through the refseq genome accession as a tuple
    genome_index = build_star_index(downloaded_gtf, downloaded_fasta)
    genome_index.view()

    // for paired-end channel, combine the downloaded paired end with original links
    // then combine with the genome_index by the refseq_genome_accession joined
    // structure is genome_refseq_accession, SRA_run_accession, downloaded_paired_end_reads, genome_index ?

    // mapping jobs - link by the refseq genome accession tag assigned to both the SRA accession and the refseq files
    // mapping_SRA_to_genome = sra_genome_mapping_pairs
    //     .join(downloaded_refseq_genomes, by: 1)
    // mapping_SRA_to_genome.view()

}
// download using SRA tools passing the SRA run accession
// paired-end reads process to split files
process download_paired_SRA_runs {
    // download each SRA run with SRAtools
    tag "${SRA_run_accession}_download"
    publishDir "${params.outdir}/sra_accessions", mode: 'copy', pattern:"*.fastq.gz"

    input:
    val(SRA_run_accession)

    output:
    path("*.fastq.gz"), emit: fastq

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
    publishDir "${params.outdir}/sra_accessions", mode: 'copy', pattern:"*.fastq.gz"

    input:
    val(SRA_run_accession)

    output:
    path("*.fastq.gz"), emit: fastq

    script:
    """
    fastq-dump --gzip ${SRA_run_accession}
    """
}

// download the Refseq genome assembly and the GTF file by combining FTP link with Refseq accession
process download_refseq_files {
    tag "${genome_refseq_accession}_download"
    publishDir "${params.outdir}/refseq_assemblies", mode: 'copy', pattern: "*.gtf"

    input:
    tuple val(genome_refseq_accession), val(genome_ftp_path)

    output:
    tuple val (genome_refseq_accession), path("*.gtf"), emit: gtf
    path("*.fna"), emit: fasta

    script:
    """
    wget ${genome_ftp_path}/${genome_refseq_accession}_genomic.gtf.gz
    wget ${genome_ftp_path}/${genome_refseq_accession}_genomic.fna.gz

    gunzip ${genome_refseq_accession}_genomic.gtf.gz
    gunzip ${genome_refseq_accession}_genomic.fna
    """
}

// build STAR index with GTF file
process build_star_index {
    tag "${genome_refseq_accession}_build_index"
    publishDir "${params.outdir}/index", mode: 'copy', pattern:"*"

    conda "envs/star.yml"

    input:
    tuple val(genome_refseq_accession), path(genome_gtf)
    path(genome_fasta)

    output:
    tuple val(genome_refseq_accession), path("star"), emit: index

    script:
    """
    STAR --runThreadN ${params.threads} \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${genome_gtf} \\
        --sjdbGTFtagExonParentTranscript mRNA \\
        --sjdbOverhang 99
    """

}



// map with HISAT2 single-end with -U

// map paired with -1,2

// quantify with featurecounts (R subread), filter by counts, filter out by select proteins
