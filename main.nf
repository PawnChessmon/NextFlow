#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline for A549 and H129 Cell Lines
 * Compares control vs overexpression conditions
 */

nextflow.enable.dsl = 2

/*
 * Pipeline parameters
 */
params.samplesheet = 'nfcore_rnaseq_samplesheet.csv'
params.genome_fasta = 'reference_files/GRCh38.primary_assembly.genome.fa'
params.gtf = 'reference_files/gencode.v46.primary_assembly.annotation.gtf'
params.outdir = 'results'
params.help = false

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

def helpMessage() {
    log.info"""
    ========================================
    RNA-seq Analysis Pipeline
    ========================================
    
    Usage:
    nextflow run main.nf --samplesheet <path_to_samplesheet.csv>
    
    Required Arguments:
    --samplesheet    Path to samplesheet CSV file
    
    Optional Arguments:
    --genome_fasta   Path to genome FASTA file (default: reference_files/GRCh38.primary_assembly.genome.fa)
    --gtf            Path to GTF annotation file (default: reference_files/gencode.v46.primary_assembly.annotation.gtf)
    --outdir         Output directory (default: results)
    --help           Show this help message
    
    Example:
    nextflow run main.nf --samplesheet nfcore_rnaseq_samplesheet.csv
    """.stripIndent()
}

/*
 * Create channels from samplesheet
 */
process PARSE_SAMPLESHEET {
    input:
    path samplesheet
    
    output:
    path "samples.txt"
    
    script:
    """
    # Skip header and extract sample info
    tail -n +2 ${samplesheet} | cut -d',' -f1,2 > samples.txt
    """
}

/*
 * Download FASTQ files from SRA
 */
process SRA_DOWNLOAD {
    tag "$sample_id"
    publishDir "${params.outdir}/01_raw_data", mode: 'copy'
    
    input:
    tuple val(sample_id), val(sra_id)
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")
    
    script:
    """
    # Download from SRA
    prefetch ${sra_id}
    fasterq-dump ${sra_id} --outfile ${sample_id}.fastq
    gzip ${sample_id}.fastq
    
    # Clean up
    rm -rf ${sra_id}
    """
}

/*
 * Quality control with FastQC
 */
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/02_fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path(fastq), path("*_fastqc.{zip,html}")
    
    script:
    """
    fastqc ${fastq} --threads ${task.cpus}
    """
}

/*
 * Read trimming with fastp
 */
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/03_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq), path(fastqc)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), path("${sample_id}_fastp.json")
    
    script:
    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${sample_id}_trimmed.fastq.gz \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_tail \\
        --cut_window_size 4 \\
        --cut_mean_quality 15 \\
        --qualified_quality_phred 15 \\
        --unqualified_percent_limit 40 \\
        --length_required 36
    """
}

/*
 * Create STAR genome index
 */
process STAR_INDEX {
    publishDir "${params.outdir}/04_genome_index", mode: 'copy'
    
    input:
    path genome_fasta
    path gtf
    
    output:
    path "star_index"
    
    script:
    """
    mkdir star_index
    
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 149 \\
        --runThreadN ${task.cpus}
    """
}

/*
 * Align reads with STAR
 */
process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/05_alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed_fastq), path(fastp_json)
    path star_index
    
    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), path("${sample_id}_Log.final.out")
    
    script:
    """
    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${trimmed_fastq} \\
        --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${sample_id}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --outSAMattributes Standard \\
        --quantMode GeneCounts
    """
}

/*
 * Index BAM files
 */
process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/05_alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(star_log)
    
    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), path(star_log)
    
    script:
    """
    samtools index ${bam}
    """
}

/*
 * Count reads with featureCounts
 */
process FEATURE_COUNTS {
    publishDir "${params.outdir}/06_counts", mode: 'copy'
    
    input:
    path bam_files
    path gtf
    
    output:
    path "gene_counts.txt"
    path "gene_counts.txt.summary"
    
    script:
    def bam_list = bam_files.collect { it.toString() }.join(' ')
    """
    featureCounts \\
        -a ${gtf} \\
        -o gene_counts.txt \\
        -T ${task.cpus} \\
        -s 2 \\
        -t exon \\
        -g gene_id \\
        ${bam_list}
    """
}

/*
 * MultiQC report
 */
process MULTIQC {
    publishDir "${params.outdir}/07_multiqc", mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

/*
 * Create R script for basic differential expression analysis
 */
process CREATE_DE_SCRIPT {
    publishDir "${params.outdir}/08_differential_expression", mode: 'copy'
    
    output:
    path "run_deseq2.R"
    
    script:
    """
    cat > run_deseq2.R << 'EOF'
# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Read count data
counts <- read.table("../06_counts/gene_counts.txt", header = TRUE, row.names = 1, sep = "\\t")
# Remove first 5 columns (Chr, Start, End, Strand, Length) to keep only count columns
counts <- counts[, 6:ncol(counts)]

# Create sample information
sample_names <- colnames(counts)
# Extract cell line and condition from sample names
cell_line <- ifelse(grepl("A549", sample_names), "A549", "H129")
condition <- ifelse(grepl("control", sample_names), "control", "overexp")

colData <- data.frame(
  sample = sample_names,
  cell_line = factor(cell_line),
  condition = factor(condition, levels = c("control", "overexp"))
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = colData,
                             design = ~ cell_line + condition)

# Pre-filtering: keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Results for condition effect (overexp vs control)
res <- results(dds, contrast = c("condition", "overexp", "control"))
res_ordered <- res[order(res\$pvalue),]

# Write results
write.csv(as.data.frame(res_ordered), "DESeq2_results_condition.csv")

# MA plot
pdf("MA_plot.pdf")
plotMA(res, main = "MA Plot: Overexpression vs Control")
dev.off()

# Volcano plot
pdf("volcano_plot.pdf")
with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano Plot"))
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
dev.off()

# PCA plot
vsd <- vst(dds, blind = FALSE)
pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = c("condition", "cell_line"))
dev.off()

# Heatmap of top variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
pdf("heatmap_top_genes.pdf")
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[, c("cell_line", "condition")]))
dev.off()

print("Analysis complete!")
print(paste("Found", sum(res\$padj < 0.05, na.rm = TRUE), "significantly differential genes"))
EOF
    """
}

/*
 * Main workflow
 */
workflow {
    // Check if reference files exist
    if (!file(params.genome_fasta).exists()) {
        error "Genome FASTA file not found: ${params.genome_fasta}. Please run ./setup_references.sh first."
    }
    
    if (!file(params.gtf).exists()) {
        error "GTF file not found: ${params.gtf}. Please run ./setup_references.sh first."
    }
    
    // Parse samplesheet
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
    PARSE_SAMPLESHEET(samplesheet_ch)
    
    // Create sample channel
    samples_ch = PARSE_SAMPLESHEET.out
        .splitCsv()
        .map { row -> tuple(row[0], row[1]) }
    
    // Download SRA data
    SRA_DOWNLOAD(samples_ch)
    
    // Quality control
    FASTQC(SRA_DOWNLOAD.out)
    
    // Trim reads
    FASTP(FASTQC.out)
    
    // Create genome index
    genome_fasta_ch = Channel.fromPath(params.genome_fasta, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    STAR_INDEX(genome_fasta_ch, gtf_ch)
    
    // Align reads
    STAR_ALIGN(FASTP.out, STAR_INDEX.out)
    
    // Index BAM files
    SAMTOOLS_INDEX(STAR_ALIGN.out)
    
    // Extract BAM files for counting
    bam_files_ch = SAMTOOLS_INDEX.out.map { sample_id, bam, bai, log -> bam }.collect()
    
    // Count features
    FEATURE_COUNTS(bam_files_ch, gtf_ch)
    
    // Collect all QC files for MultiQC
    qc_files_ch = FASTQC.out.map { sample_id, fastq, qc -> qc }
        .mix(FASTP.out.map { sample_id, fastq, json -> json })
        .mix(SAMTOOLS_INDEX.out.map { sample_id, bam, bai, log -> log })
        .mix(FEATURE_COUNTS.out[1])
        .collect()
    
    // Generate MultiQC report
    MULTIQC(qc_files_ch)
    
    // Create DE analysis script
    CREATE_DE_SCRIPT()
}