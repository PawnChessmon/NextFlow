#!/bin/bash

# Script to download and set up reference files for the RNA-seq pipeline
# Run this script from your NextFlow directory

echo "Setting up reference files for RNA-seq pipeline"
echo "==============================================="

# Create reference_files directory if it doesn't exist
mkdir -p reference_files
cd reference_files

echo "1. Downloading reference genome (GRCh38)..."
if [ ! -f "GRCh38.primary_assembly.genome.fa" ]; then
    echo "   Downloading genome file..."
    curl -L -o GRCh38.primary_assembly.genome.fa.gz \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    
    echo "   Extracting genome file..."
    gunzip GRCh38.primary_assembly.genome.fa.gz
    echo "   ✅ Genome file ready: $(pwd)/GRCh38.primary_assembly.genome.fa"
else
    echo "   ✅ Genome file already exists: $(pwd)/GRCh38.primary_assembly.genome.fa"
fi

echo ""
echo "2. Downloading GTF annotation..."
if [ ! -f "gencode.v46.primary_assembly.annotation.gtf" ]; then
    echo "   Downloading GTF file..."
    curl -L -o gencode.v46.primary_assembly.annotation.gtf.gz \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
    
    echo "   Extracting GTF file..."
    gunzip gencode.v46.primary_assembly.annotation.gtf.gz
    echo "   ✅ GTF file ready: $(pwd)/gencode.v46.primary_assembly.annotation.gtf"
else
    echo "   ✅ GTF file already exists: $(pwd)/gencode.v46.primary_assembly.annotation.gtf"
fi

cd ..

echo ""
echo "Reference files are ready!"
echo "=========================="
echo "Genome: $(pwd)/reference_files/GRCh38.primary_assembly.genome.fa"
echo "GTF:    $(pwd)/reference_files/gencode.v46.primary_assembly.annotation.gtf"
echo ""
echo "You can now run the pipeline with:"
echo "./nextflow run main.nf \\"
echo "  --fasta reference_files/GRCh38.primary_assembly.genome.fa \\"
echo "  --gtf reference_files/gencode.v46.primary_assembly.annotation.gtf"
echo ""
echo "File sizes:"
ls -lh reference_files/