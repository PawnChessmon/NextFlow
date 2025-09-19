# RNA-seq Analysis Pipeline

This repository contains a complete Nextflow pipeline for RNA-seq analysis of A549 and H129 cell lines comparing control vs overexpression conditions.

## Dataset Information

- **Cell Lines**: A549 and H129
- **Conditions**: Control vs Overexpression
- **Samples**: 12 total (6 per cell line, 3 replicates per condition)
- **Library Type**: Strand-specific (reverse stranded)
- **Data Source**: SRA public datasets

## Pipeline Features

### Analysis Steps
1. **SRA Download**: Automatically downloads FASTQ files from SRA
2. **Quality Control**: FastQC analysis of raw reads
3. **Read Trimming**: fastp for adapter removal and quality trimming
4. **Genome Indexing**: STAR genome index creation
5. **Read Alignment**: STAR alignment to GRCh38
6. **Read Counting**: featureCounts for gene expression quantification
7. **Quality Report**: MultiQC comprehensive report
8. **Differential Expression**: DESeq2 R script generation

### Key Files
- `main.nf` - Main Nextflow pipeline script
- `nextflow.config` - Configuration file with resource settings
- `nfcore_rnaseq_samplesheet.csv` - Sample information for the pipeline
- `setup_references.sh` - Script to download reference genome and annotations

## Quick Start

### 1. Prerequisites
Install the following software:
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

### 2. Setup References
```bash
./setup_references.sh
```

### 3. Run the Pipeline

#### Using Conda (Recommended)
```bash
nextflow run main.nf -profile standard
```

#### Using Docker
```bash
nextflow run main.nf -profile docker
```

#### Using Singularity
```bash
nextflow run main.nf -profile singularity
```

#### Custom Parameters
```bash
nextflow run main.nf \
    --samplesheet nfcore_rnaseq_samplesheet.csv \
    --outdir my_results \
    -profile docker
```

### 4. Results
The pipeline will create a `results/` directory with:
- `01_raw_data/` - Downloaded FASTQ files
- `02_fastqc/` - Quality control reports
- `03_trimmed/` - Trimmed FASTQ files
- `04_genome_index/` - STAR genome index
- `05_alignment/` - BAM alignment files
- `06_counts/` - Gene count matrix
- `07_multiqc/` - Comprehensive quality report
- `08_differential_expression/` - R script for DE analysis

### 5. Differential Expression Analysis
After the pipeline completes, run the generated R script:
```bash
cd results/08_differential_expression
Rscript run_deseq2.R
```

## Resource Requirements

### Minimum System Requirements
- **CPU**: 8+ cores recommended
- **RAM**: 32+ GB recommended
- **Storage**: ~100 GB free space
- **Time**: ~6-12 hours for 12 samples

### Compute Profiles
- `standard`: Uses Conda for software management
- `docker`: Uses Docker containers
- `singularity`: Uses Singularity containers
- `test`: Reduced resources for testing

## Alternative: nf-core/rnaseq

You can also use the established nf-core/rnaseq pipeline:
```bash
nextflow run nf-core/rnaseq \
    --input nfcore_rnaseq_samplesheet.csv \
    --genome GRCh38 \
    --stranded reverse \
    --aligner star_salmon \
    --save_reference \
    --deseq2_vst \
    -profile docker
```

## Sample Design

| Sample ID | Cell Line | Condition | SRA Accession |
|-----------|-----------|-----------|---------------|
| A549_control_1 | A549 | Control | SRR32033552 |
| A549_control_2 | A549 | Control | SRR32033551 |
| A549_control_3 | A549 | Control | SRR32033550 |
| A549_overexp_1 | A549 | Overexpression | SRR32033549 |
| A549_overexp_2 | A549 | Overexpression | SRR32033548 |
| A549_overexp_3 | A549 | Overexpression | SRR32033547 |
| H129_control_1 | H129 | Control | SRR32033546 |
| H129_control_2 | H129 | Control | SRR32033545 |
| H129_control_3 | H129 | Control | SRR32033544 |
| H129_overexp_1 | H129 | Overexpression | SRR32033543 |
| H129_overexp_2 | H129 | Overexpression | SRR32033542 |
| H129_overexp_3 | H129 | Overexpression | SRR32033541 |

## Getting Started

1. **Setup references** (optional if using --genome GRCh38):
   ```bash
   ./setup_references.sh
   ```

2. **Run with nf-core/rnaseq**:
   ```bash
   nextflow run nf-core/rnaseq \
     --input nfcore_rnaseq_samplesheet.csv \
     --genome GRCh38 \
     --stranded reverse \
     --aligner star_salmon \
     --save_reference \
     --deseq2_vst
   ```

3. **Or use on Seqera Platform** with the provided samplesheet and parameters.