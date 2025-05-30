#!/usr/bin/bash

# Setup script to install dependencies and configure the RNA-Seq pipeline environment

set -e  # Exit on any error

# Define directories and variables
BASE_DIR="/home/program/RNASeq/Vrunda_Kallisto_Sleuth"
CONDA_ENV="rna_quanti"
CONDA_SH="/home/program/miniconda3/etc/profile.d/conda.sh"
LOG_FILE="$BASE_DIR/setup.log"

# Create base directory and log file
mkdir -p "$BASE_DIR"
touch "$LOG_FILE"

echo "Starting setup process..." | tee -a "$LOG_FILE"

# Checklist function
print_checklist() {
    local step=$1
    local status=$2
    if [ "$status" = "success" ]; then
        echo "[✔] $step completed successfully" | tee -a "$LOG_FILE"
    else
        echo "[✘] $step failed. Check $LOG_FILE for details." | tee -a "$LOG_FILE"
        exit 1
    fi
}

# Step 1: Install Miniconda
echo "Checking for Miniconda..." | tee -a "$LOG_FILE"
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..." | tee -a "$LOG_FILE"
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh >> "$LOG_FILE" 2>&1
    bash miniconda.sh -b -p /home/program/miniconda3 >> "$LOG_FILE" 2>&1
    rm miniconda.sh
    [ $? -eq 0 ] && print_checklist "Miniconda installation" "success" || print_checklist "Miniconda installation" "failed"
else
    print_checklist "Miniconda already installed" "success"
fi

# Initialize Conda
source "$CONDA_SH"

# Step 2: Create Conda environment
echo "Creating Conda environment: $CONDA_ENV..." | tee -a "$LOG_FILE"
if ! conda env list | grep -q "$CONDA_ENV"; then
    conda create -n "$CONDA_ENV" python=3.8 -y >> "$LOG_FILE" 2>&1
    [ $? -eq 0 ] && print_checklist "Conda environment creation" "success" || print_checklist "Conda environment creation" "failed"
else
    print_checklist "Conda environment already exists" "success"
fi

# Activate Conda environment
conda activate "$CONDA_ENV"
[ $? -eq 0 ] && print_checklist "Conda environment activation" "success" || print_checklist "Conda environment activation" "failed"

# Step 3: Install dependencies
echo "Installing dependencies..." | tee -a "$LOG_FILE"
conda install -c bioconda sra-tools kallisto fastqc trim-galore multiqc -y >> "$LOG_FILE" 2>&1
conda install -c r r r-base r-essentials -y >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "Bioconda and R dependencies installation" "success" || print_checklist "Bioconda and R dependencies installation" "failed"

# Step 4: Install R packages
echo "Installing R packages (sleuth, dplyr)..." | tee -a "$LOG_FILE"
Rscript -e "install.packages(c('dplyr', 'sleuth'), repos='http://cran.rstudio.com/')" >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "R packages installation" "success" || print_checklist "R packages installation" "failed"

# Step 5: Install Python pandas
echo "Installing Python pandas..." | tee -a "$LOG_FILE"
pip install pandas >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "Python pandas installation" "success" || print_checklist "Python pandas installation" "failed"

# Step 6: Create directory structure
echo "Setting up directory structure..." | tee -a "$LOG_FILE"
mkdir -p "$BASE_DIR/input" "$BASE_DIR/sra_data" "$BASE_DIR/fastqc" "$BASE_DIR/trimmed" \
         "$BASE_DIR/kallisto_out" "$BASE_DIR/targetf/For_sleuth" "$BASE_DIR/sleuth_out" \
         "$BASE_DIR/log_files" "$BASE_DIR/multiqc" "$BASE_DIR/reference"
[ $? -eq 0 ] && print_checklist "Directory structure creation" "success" || print_checklist "Directory structure creation" "failed"

# Step 7: Download reference transcriptome
echo "Downloading reference transcriptome..." | tee -a "$LOG_FILE"
if [ ! -f "$BASE_DIR/reference/transcriptome.fa" ]; then
    wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O "$BASE_DIR/reference/transcriptome.fa.gz" >> "$LOG_FILE" 2>&1
    gunzip "$BASE_DIR/reference/transcriptome.fa.gz" >> "$LOG_FILE" 2>&1
    [ $? -eq 0 ] && print_checklist "Reference transcriptome download" "success" || print_checklist "Reference transcriptome download" "failed"
else
    print_checklist "Reference transcriptome already exists" "success"
fi

# Step 8: Build Kallisto index
echo "Building Kallisto index..." | tee -a "$LOG_FILE"
if [ ! -f "$BASE_DIR/reference/transcriptome.idx" ]; then
    kallisto index -i "$BASE_DIR/reference/transcriptome.idx" "$BASE_DIR/reference/transcriptome.fa" >> "$LOG_FILE" 2>&1
    [ $? -eq 0 ] && print_checklist "Kallisto index creation" "success" || print_checklist "Kallisto index creation" "failed"
else
    print_checklist "Kallisto index already exists" "success"
fi

echo "Setup completed successfully!" | tee -a "$LOG_FILE"
conda deactivate
