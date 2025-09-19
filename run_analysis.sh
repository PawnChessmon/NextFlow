#!/bin/bash

# Quick setup and run script for RNA-seq analysis
# This script will setup everything and run the pipeline

set -e

echo "========================================="
echo "RNA-seq Analysis Pipeline Setup"
echo "========================================="

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "❌ Nextflow is not installed"
    echo "Please install Nextflow first:"
    echo "curl -s https://get.nextflow.io | bash"
    echo "sudo mv nextflow /usr/local/bin/"
    exit 1
else
    echo "✅ Nextflow is installed: $(nextflow -version | head -n 1)"
fi

# Check for container runtime
PROFILE=""
if command -v docker &> /dev/null; then
    echo "✅ Docker found - will use Docker profile"
    PROFILE="docker"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity found - will use Singularity profile"
    PROFILE="singularity"
elif command -v conda &> /dev/null; then
    echo "✅ Conda found - will use Conda profile"
    PROFILE="standard"
else
    echo "❌ No container runtime found (Docker/Singularity) or Conda"
    echo "Please install one of: Docker, Singularity, or Conda"
    exit 1
fi

echo ""
echo "========================================="
echo "Setting up reference files"
echo "========================================="

# Run reference setup
if [ ! -f "reference_files/GRCh38.primary_assembly.genome.fa" ]; then
    echo "📥 Downloading reference genome and annotations..."
    ./setup_references.sh
else
    echo "✅ Reference files already exist"
fi

echo ""
echo "========================================="
echo "Pipeline Information"
echo "========================================="
echo "📊 Dataset: 12 RNA-seq samples (A549 and H129 cell lines)"
echo "🔬 Analysis: Control vs Overexpression comparison"  
echo "⚙️  Profile: $PROFILE"
echo "📁 Output: results/"
echo "⏱️  Estimated time: 6-12 hours"

echo ""
read -p "Do you want to start the analysis now? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "========================================="
    echo "Starting RNA-seq Analysis"
    echo "========================================="
    
    # Run the pipeline
    nextflow run main.nf -profile $PROFILE -resume
    
    echo ""
    echo "========================================="
    echo "Analysis Complete!"
    echo "========================================="
    echo "📁 Results are in: results/"
    echo "📊 MultiQC report: results/07_multiqc/multiqc_report.html"
    echo "🧬 Gene counts: results/06_counts/gene_counts.txt"
    echo ""
    echo "To run differential expression analysis:"
    echo "cd results/08_differential_expression"
    echo "Rscript run_deseq2.R"
    
else
    echo ""
    echo "Analysis not started. To run manually:"
    echo "nextflow run main.nf -profile $PROFILE"
fi