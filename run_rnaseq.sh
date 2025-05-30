#!/usr/bin/bash

# RNA-Seq pipeline script
# Usage: bash run_rnaseq.sh --input <SRA_IDs_file or FASTQ_dir>

set -e  # Exit on any error

# Parse command-line arguments
INPUT=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift 2 ;;
        *) echo "Unknown option: $1" | tee -a "$LOG_FILE"; exit 1 ;;
    esac
done

if [ -z "$INPUT" ]; then
    echo "Error: --input is required (SRA IDs file or FASTQ directory)" | tee -a "$LOG_FILE"
    exit 1
fi

# Define directories and variables
BASE_DIR="/home/program/RNASeq/Vrunda_Kallisto_Sleuth"
CONDA_ENV="rna_quanti"
CONDA_SH="/home/program/miniconda3/etc/profile.d/conda.sh"
LOG_FILE="$BASE_DIR/log_files/pipeline.log"
INDEX="$BASE_DIR/reference/transcriptome.idx"

# Create log directory and file
mkdir -p "$BASE_DIR/log_files"
touch "$LOG_FILE"

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

# Initialize Conda
source "$CONDA_SH"
conda activate "$CONDA_ENV"
[ $? -eq 0 ] && print_checklist "Conda environment activation" "success" || print_checklist "Conda environment activation" "failed"

# Step 1: Input processing
echo "Processing input..." | tee -a "$LOG_FILE"
if [ -f "$INPUT" ]; then
    # Input is a file with SRA IDs
    echo "Downloading SRA data..." | tee -a "$LOG_FILE"
    while read -r sra_id; do
        fasterq-dump "$sra_id" -O "$BASE_DIR/sra_data" --split-files >> "$LOG_FILE" 2>&1
        [ $? -eq 0 ] || { print_checklist "SRA download for $sra_id" "failed"; exit 1; }
        gzip "$BASE_DIR/sra_data/$sra_id"*.fastq >> "$LOG_FILE" 2>&1
    done < "$INPUT"
    print_checklist "SRA data download" "success"
elif [ -d "$INPUT" ]; then
    # Input is a directory with FASTQ files
    echo "Copying FASTQ files from $INPUT..." | tee -a "$LOG_FILE"
    cp "$INPUT"/*.fastq.gz "$BASE_DIR/sra_data/" >> "$LOG_FILE" 2>&1
    [ $? -eq 0 ] && print_checklist "FASTQ files copy" "success" || print_checklist "FASTQ files copy" "failed"
else
    echo "Error: Input must be a file with SRA IDs or a directory with FASTQ files" | tee -a "$LOG_FILE"
    exit 1
fi

# Detect paired-end or single-end
PAIRED_END=false
for fastq in "$BASE_DIR/sra_data"/*.fastq.gz; do
    if [[ "$fastq" == *"_1.fastq.gz" ]] && [[ -f "${fastq/_1.fastq.gz/_2.fastq.gz}" ]]; then
        PAIRED_END=true
        break
    fi
done
echo "Detected data type: $( [ "$PAIRED_END" = true ] && echo "Paired-end" || echo "Single-end" )" | tee -a "$LOG_FILE"

# Step 2: Quality control with FastQC
echo "Running FastQC..." | tee -a "$LOG_FILE"
mkdir -p "$BASE_DIR/fastqc"
fastqc "$BASE_DIR/sra_data"/*.fastq.gz -o "$BASE_DIR/fastqc" >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "FastQC analysis" "success" || print_checklist "FastQC analysis" "failed"

# Step 3: Adapter trimming with Trim Galore
echo "Running Trim Galore..." | tee -a "$LOG_FILE"
mkdir -p "$BASE_DIR/trimmed"
for fastq in "$BASE_DIR/sra_data"/*.fastq.gz; do
    if [ "$PAIRED_END" = true ] && [[ "$fastq" == *"_1.fastq.gz" ]]; then
        base=$(basename "$fastq" _1.fastq.gz)
        if [[ -f "$BASE_DIR/sra_data/${base}_2.fastq.gz" ]]; then
            trim_galore --paired "$BASE_DIR/sra_data/${base}_1.fastq.gz" "$BASE_DIR/sra_data/${base}_2.fastq.gz" \
                        -o "$BASE_DIR/trimmed" >> "$LOG_FILE" 2>&1
            [ $? -eq 0 ] || { print_checklist "Trim Galore for $base" "failed"; exit 1; }
        fi
    elif [ "$PAIRED_END" = false ]; then
        trim_galore "$fastq" -o "$BASE_DIR/trimmed" >> "$LOG_FILE" 2>&1
        [ $? -eq 0 ] || { print_checklist "Trim Galore for $(basename "$fastq")" "failed"; exit 1; }
    fi
done
print_checklist "Trim Galore processing" "success"

# Step 4: Kallisto quantification
echo "Running Kallisto quantification..." | tee -a "$LOG_FILE"
mkdir -p "$BASE_DIR/kallisto_out"
for fastq in "$BASE_DIR/trimmed"/*.fastq.gz; do
    base=$(basename "$fastq" .fastq.gz)
    out_dir="$BASE_DIR/kallisto_out/$base"
    if [ "$PAIRED_END" = true ] && [[ "$base" == *"_1" ]]; then
        base_pair=$(echo "$base" | sed 's/_1$//')
        if [[ -f "$BASE_DIR/trimmed/${base_pair}_2.fastq.gz" ]]; then
            kallisto quant -i "$INDEX" -o "$out_dir" --bias \
                          "$BASE_DIR/trimmed/${base_pair}_1.fastq.gz" "$BASE_DIR/trimmed/${base_pair}_2.fastq.gz" \
                          >> "$LOG_FILE" 2>&1
            [ $? -eq 0 ] || { print_checklist "Kallisto for $base_pair" "failed"; exit 1; }
        fi
    elif [ "$PAIRED_END" = false ]; then
        kallisto quant -i "$INDEX" -o "$out_dir" --single -l 200 -s 20 --bias "$fastq" \
                      >> "$LOG_FILE" 2>&1
        [ $? -eq 0 ] || { print_checklist "Kallisto for $base" "failed"; exit 1; }
    fi
done
print_checklist "Kallisto quantification" "success"

# Step 5: Generate metadata for Sleuth
echo "Generating metadata for Sleuth..." | tee -a "$LOG_FILE"
cat > "$BASE_DIR/generate_metadata.py" << 'EOF'
import pandas as pd
import os
import glob

base_dir = "/home/program/RNASeq/Vrunda_Kallisto_Sleuth"
kallisto_dir = f"{base_dir}/kallisto_out"
output_dir = f"{base_dir}/targetf/For_sleuth"
os.makedirs(output_dir, exist_ok=True)

samples = [os.path.basename(d) for d in glob.glob(f"{kallisto_dir}/*") if os.path.isdir(d)]
data = []
for sample in samples:
    # Adjust condition logic based on your naming convention
    condition = "control" if "control" in sample.lower() else "treatment"
    data.append({"sample": sample, "condition": condition, "path": f"{kallisto_dir}/{sample}"})

df = pd.DataFrame(data)
df.to_csv(f"{output_dir}/metadata.csv", index=False)
print("Metadata generated at:", f"{output_dir}/metadata.csv")
EOF
python3 "$BASE_DIR/generate_metadata.py" >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "Metadata generation" "success" || print_checklist "Metadata generation" "failed"

# Step 6: Run Sleuth
echo "Running Sleuth differential expression analysis..." | tee -a "$LOG_FILE"
cat > "$BASE_DIR/sleuth.R" << 'EOF'
library(sleuth)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript sleuth.R input_target.csv sleuthO_output.csv")
}

targetf <- args[1]
sleuth_out <- args[2]

s2c <- read.csv(targetf, header = TRUE, stringsAsFactors = FALSE)
cond <- unique(s2c$condition)
if (length(cond) != 2) {
  stop("Exactly two conditions are required for differential expression analysis")
}

s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, ref = cond[1])
group_var <- paste0("condition", cond[2])

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, num_cores = 4)
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")
so <- sleuth_wt(so, which_beta = group_var, which_model = "full")

sleuth_table <- sleuth_results(so, group_var, "wt", show_all = FALSE)
write.csv(sleuth_table, file = sleuth_out, row.names = FALSE)
EOF
mkdir -p "$BASE_DIR/sleuth_out"
Rscript "$BASE_DIR/sleuth.R" "$BASE_DIR/targetf/For_sleuth/metadata.csv" "$BASE_DIR/sleuth_out/sleuthO_results.csv" >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "Sleuth analysis" "success" || print_checklist "Sleuth analysis" "failed"

# Step 7: Generate MultiQC report
echo "Generating MultiQC report..." | tee -a "$LOG_FILE"
multiqc "$BASE_DIR" -o "$BASE_DIR/multiqc" >> "$LOG_FILE" 2>&1
[ $? -eq 0 ] && print_checklist "MultiQC report generation" "success" || print_checklist "MultiQC report generation" "failed"

echo "Pipeline completed successfully!" | tee -a "$LOG_FILE"
conda deactivate
