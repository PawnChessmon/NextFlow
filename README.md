# RNA-seq Analysis Pipeline

This repository contains setup files for RNA-seq analysis of A549 and H129 cell lines comparing control vs overexpression conditions.

## Dataset Information

- **Cell Lines**: A549 and H129
- **Conditions**: Control vs Overexpression
- **Samples**: 12 total (6 per cell line, 3 replicates per condition)
- **Library Type**: Strand-specific (reverse stranded)
- **Data Source**: SRA public datasets

## Files

### Sample Information
- `ss.csv` - Original sample sheet format
- `nfcore_rnaseq_samplesheet.csv` - Formatted for nf-core/rnaseq pipeline
- `SRR_Acc_List.txt` - SRA accession numbers for data download

### Reference Setup
- `setup_references.sh` - Script to download GRCh38 genome and GENCODE v46 annotations

## Usage with nf-core/rnaseq

### Recommended Parameters
```
--input nfcore_rnaseq_samplesheet.csv
--genome GRCh38
--stranded reverse
--aligner star_salmon
--save_reference true
--deseq2_vst true
```

### Expected Analysis
- Quality control (FastQC, MultiQC)
- Alignment with STAR
- Quantification with Salmon
- Differential expression with DESeq2

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