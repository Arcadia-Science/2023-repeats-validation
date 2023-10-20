#! /usr/bin/env nextflow

// Description
// Download RNAseq experiments from brain/muscle tissues and species genome/GTF files to map and quantify expression of genes
// Assess the expression of repeat-expansion homologs identified in the species through protein sequence and structural comparison methods


// outputs the count table for each mapping job and the GTF file for each genome
nextflow.enable.dsl=2

log.info """\

DOWNLOAD, MAP, AND PROCESS TRANSCRIPTOMES AGAINST SPECIES GENOMES
TO QUANTIFY EXPRESSION OF SELECT GENES
=========================================
samples            : $params.samples
outdir             : $params.outdir
threads            : $params.threads
"""

// channels for the CSV of samples for SRA accessions and genome Refseq FTP links to download, mapping queries
// CSV to provide to select proteins from to quantify expression in the final table
// if SRA accession is single or paired end, split up into different output channels to process the command differently

// define channels and parameters
sample_csv = Channel.fromPath(params.samples)
    .splitCsv(header:true)


// split SRA accessions based on if library_layout is SINGLE or PAIRED
// add the corresponding refseq_accession so carries through as a tuple
paired_end_samples = sample_csv.filter { it.library_layout == "PAIRED" }.map { [ it.genome_refseq_accession, it.SRA_run_accession ] }
single_end_samples = sample_csv.filter { it.library_layout == "SINGLE" }.map { [ it.genome_refseq_accession, it.SRA_run_accession ] }

// channel linking genome refseq accession to FTP download path
refseq_genomes = sample_csv.map { [ it.genome_refseq_accession, it.genome_ftp_path ]}
    .unique() //deduplicate because samplesheet is set up for mapping pairs

workflow {
    // download SRA files
    downloaded_paired_end_reads = download_paired_SRA_runs(paired_end_samples)
    downloaded_single_end_reads = download_single_SRA_runs(single_end_samples)

    // combined downloaded reads
    combined_reads = downloaded_paired_end_reads.concat(downloaded_single_end_reads)

    // download Refseq files
    download_refseq_files(refseq_genomes)
    downloaded_fasta = download_refseq_files.out.fasta
    downloaded_gtf = download_refseq_files.out.gtf

    // build STAR index, include the GTF and FASTA for exons
    // carry through the refseq genome accession as a tuple
    reference_files = downloaded_fasta
        .join(downloaded_gtf)
    genome_index = build_star_index(reference_files)

    // prepare mapping channels with correct refseq accession : SRA accession pairings
    mapping_samples = combined_reads
        .combine(genome_index, by: 0)

    mapped_BAMS = star_mapping(mapping_samples)

    // htseq count - combine mapped_BAMS with the corresponding genome GTF
    htseq_input = mapped_BAMS
        .combine(downloaded_gtf, by: 0)

    htseq_counts = htseq_count(htseq_input)

    // parse GTF files to simple tsv
    gtf_tables = parse_gtf(downloaded_gtf)

}
// download using SRA tools passing the SRA run accession
// paired-end reads process to split files
process download_paired_SRA_runs {
    // download each SRA run with SRAtools
    tag "${SRA_run_accession}_download"

    errorStrategy 'ignore' // ignore failed downloads

    conda "envs/sratoolkit.yml"

    input:
    tuple val(genome_refseq_accession), val(SRA_run_accession)

    output:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path("*.fastq.gz"), emit: fastq

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

    errorStrategy 'ignore' // ignore failed downloads

    conda "envs/sratoolkit.yml"

    input:
    tuple val(genome_refseq_accession), val(SRA_run_accession)

    output:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path("*.fastq.gz"), emit: fastq

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
    tuple val (genome_refseq_accession), path("*.fna"), emit: fasta

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
    label 'process_high'
    publishDir "${params.outdir}/index", mode: 'copy', pattern:"*"

    conda "envs/star.yml"

    input:
    tuple val(genome_refseq_accession), path(genome_fasta), path(genome_gtf)

    output:
    tuple val(genome_refseq_accession), path("*star"), emit: index

    script:
    """
    STAR --runThreadN ${params.threads} \\
        --runMode genomeGenerate \\
        --genomeDir ${genome_refseq_accession}_star/ \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${genome_gtf} \\
        --sjdbGTFtagExonParentTranscript mRNA \\
        --sjdbOverhang 99 \\
        --limitGenomeGenerateRAM=60000000000
    """

}

// map with STAR single end reads, index with samtools
process star_mapping {
    tag "${genome_refseq_accession}-vs-${SRA_run_accession}_mapping"

    conda "envs/star_samtools.yml"

    input:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path(reads), path(index)

    output:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path("*.bam"), path("*.bai"), emit: mapping_file

    script:
    """
    STAR --runThreadN ${params.threads} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat \\
        --outFilterType BySJout \\
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \\
         --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \\
         --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \\
         --alignIntronMax 1000000 --alignMatesGapMax 1000000  \\
         --outSAMattributes NH HI NM MD --outSAMtype BAM      \\
         SortedByCoordinate --outFileNamePrefix ${genome_refseq_accession}_vs_${SRA_run_accession}

    samtools index ${genome_refseq_accession}_vs_${SRA_run_accession}Aligned.sortedByCoord.out.bam
    """
}

// quantify with htseq count
process htseq_count {
    tag "${genome_refseq_accession}-vs-${SRA_run_accession}_count"
    publishDir "${params.outdir}/counts", mode: 'copy', pattern:"*"

    conda "envs/htseq.yml"

    input:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path(bam_file), path(bai_file), path(genome_gtf)

    output:
    tuple val(genome_refseq_accession), val(SRA_run_accession), path("*.htseq"), emit: htseq_counts

    script:
    """
    htseq-count -f bam ${bam_file} ${genome_gtf} -r pos > ${genome_refseq_accession}_vs_${SRA_run_accession}.htseq
    """

}

// parse GTF files to TSV with gene_id : protein_id
process parse_gtf {
    tag "${genome_refseq_accession}_gtf_parse"
    publishDir "${params.outdir}/gtf_tables", mode: 'copy', pattern:"*.tsv"

    input:
    tuple val(genome_refseq_accession), path(gtf)

    output:
    path("*.tsv"), emit: gtf_table

    script:
    """
    python3 ${baseDir}/bin/parse-gtf-files.py ${gtf} ${genome_refseq_accession}_gtf_table.tsv
    """
}
