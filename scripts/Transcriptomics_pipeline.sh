#!/bin/bash
# Usage: bash qc_trim_pipeline.sh reads output_dir

# Initialize conda for this shell
eval "$(conda shell.bash hook)"


# activate envirnment once 
conda activate Transcriptomics_pipeline

############------------------------############

## To run this script using Arguements 
# Arguement1 is reads folder
# Arguement2 is Outputs 
# Arguement3 is reference folder 


#------------------------------------#

INPUT_DIR=$1   # Arguement1
OUT_DIR=$2

if [ $# -lt 2 ]; then
       echo "Usage: bash qc_trim_pipeline.sh <reads_dir> <output_dir>"
     exit 1
     fi
     mkdir -p "$OUT_DIR/fastqc_raw" "$OUT_DIR/trimmed"
     echo "=== QC & Trimming Pipeline Started ==="
     echo "Input folder : $INPUT_DIR"
     echo "Output folder: $OUT_DIR"
     echo

# -----------------------------

# Step 1: FastQC on raw reads

# -----------------------------
echo "Running FastQC on raw reads..."

for fq in "$INPUT_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)

    # Check if FastQC HTML report already exists

    if [ -f "$OUT_DIR/fastqc_raw/${base}_fastqc.html" ]; then

        echo "FastQC already done for $base. Skipping..."
    else
        echo "Running FastQC for $base..."
        fastqc -o "$OUT_DIR/fastqc_raw" -t 8 "$fq"
    fi
done

     echo "FastQC completed."
     echo

     # ----------------------------
     # Step 2: Trim Galore
     # -----------------------------
     echo "Running Trim Galore..."

     # Activate Trim Galore env

     

########## Trimming for the paired ends #############  

    for f in "$INPUT_DIR"/*_1.fastq.gz; do

    base=$(basename "$f" _1.fastq.gz)

    r="$INPUT_DIR/${base}_2.fastq.gz"

    # Check if both R1 and R2 exist

    if [ ! -f "$r" ]; then
        echo "Paired file missing for $base, skipping..."
        continue
    fi

    ## check files already exist 

     if [ -f "$OUT_DIR/trimmed/${base}_1_val_1.fq.gz" ] && [ -f "$OUT_DIR/trimmed/${base}_2_val_2.fq.gz" ]; then
     echo "Trimmed paired-end files already exist for $base. Skipping..."
     continue
     fi


    echo "Processing paired-end sample: $base"
    trim_galore --paired -o "$OUT_DIR/trimmed" "$f" "$r" -j 8
    done

     # === Process single-end reads ===

     for file in "$INPUT_DIR"/*.fastq.gz; do

     # Skip files already handled as paired-end

     if [[ "$file" == *_1*.fastq.gz ]] || [[ "$file" == *_2*.fastq.gz ]]; then
     continue
     fi



    ##=== if single files already exist ==

      base=$(basename "$file" .fastq.gz)

    # Check if single-end trimmed file already exists

    if [ -f "$OUT_DIR/trimmed/${base}_trimmed.fq.gz" ]; then
        echo "Trimmed single-end file exists for $base. Skipping..."
        continue
    fi

    ## == trimming of single end reads == 

     base=$(basename "$file" .fastq.gz)
     echo "Processing single-end sample: $base"
     trim_galore -o "$OUT_DIR/trimmed" "$file" -j 8
     done

     echo "upto trimming done"
     
     
     # -----------------------------
# Step 3: HISAT2 indexing and alignment
# -----------------------------


REF_DIR=$3  # Third argument: folder containing reference genome (.fa) and GTF (.gtf)
mkdir -p "$OUT_DIR/hisat2_index" "$OUT_DIR/alignment"

# Check if reference fasta exists

GENOME_FA=$(echo "$REF_DIR"/*.fa "$REF_DIR"/*.fasta | awk '{print $1}')
if [ -z "$GENOME_FA" ]; then
    echo "Reference genome FASTA not found in $REF_DIR"
    exit 1
fi

# Build HISAT2 index if not already present
INDEX_PREFIX="$OUT_DIR/hisat2_index/genome_index"
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then      # ! is used here if file missing
    echo "Building HISAT2 index..."
    hisat2-build "$GENOME_FA" "$INDEX_PREFIX"
else
    echo "HISAT2 index already exists. Skipping..."
fi


# Alignment of trimmed reads
echo "Starting alignment with HISAT2..."

# Paired-end alignment

for f in "$OUT_DIR/trimmed/"*_1_val_1.fq.gz; do
    base=$(basename "$f" _1_val_1.fq.gz)
    r="$OUT_DIR/trimmed/${base}_2_val_2.fq.gz"

    # Skip if BAM already exists

    if [ -f "$OUT_DIR/alignment/${base}.bam" ]; then
        echo "Alignment already done for $base. Skipping..."
        continue
    fi

    echo "Aligning paired-end sample: $base"

    hisat2 -x "$INDEX_PREFIX" -1 "$f" -2 "$r" -S "$OUT_DIR/alignment/${base}.sam" -p 8
    samtools view -bS "$OUT_DIR/alignment/${base}.sam" | samtools sort -o "$OUT_DIR/alignment/${base}.bam"
    samtools index "$OUT_DIR/alignment/${base}.bam"
    rm "$OUT_DIR/alignment/${base}.sam"
done

# Single-end alignment
for f in "$OUT_DIR/trimmed/"*_trimmed.fq.gz; do
    base=$(basename "$f" _trimmed.fq.gz)

    # Skip if BAM already exists
    if [ -f "$OUT_DIR/alignment/${base}.bam" ]; then
        echo "Alignment already done for $base. Skipping..."
        continue
    fi

    echo "Aligning single-end sample: $base"
    hisat2 -x "$INDEX_PREFIX" -U "$f" -S "$OUT_DIR/alignment/${base}.sam" -p 8
    samtools view -bS "$OUT_DIR/alignment/${base}.sam" | samtools sort -o "$OUT_DIR/alignment/${base}.bam"
    samtools index "$OUT_DIR/alignment/${base}.bam"
    rm "$OUT_DIR/alignment/${base}.sam"
done


echo "HISAT2 alignment completed."





###### starting of Transcriptome assembly 




# -----------------------------
# Step 4: StringTie assembly (GTF per sample only) + count matrix generation
# -----------------------------


# Detect reference annotation (.gtf or .gff)

REF_ANNOT=$(ls "$REF_DIR"/*.{gtf,gff,gff3} 2>/dev/null | head -n 1)


# Check if a file was found
if [ -z "$REF_ANNOT" ]; then
echo "Reference annotation (.gtf or .gff) not found in $REF_DIR"
exit 1
fi

echo "Using annotation file: $REF_ANNOT"

#######################

# Output directories
ASSEMBLY_DIR="$OUT_DIR/stringtie_assembly"
COUNT_DIR="$OUT_DIR/stringtie_counts"
mkdir -p "$ASSEMBLY_DIR" "$COUNT_DIR"

echo "Starting StringTie transcript assembly (GTFs) using: $(basename "$REF_ANNOT")"

for bam in "$OUT_DIR/alignment/"*.bam; do
base=$(basename "$bam" .bam)

# Skip if GTF already exists
if [ -f "$ASSEMBLY_DIR/${base}.gtf" ]; then
echo "GTF already exists for $base. Skipping..."
continue
fi

echo "Processing sample: $base"
stringtie "$bam" \
-p 8 \
-G "$REF_ANNOT" \
-o "$ASSEMBLY_DIR/${base}.gtf" \
-e
done
echo "✅ Per-sample GTF files generated."


# -----------------------------
# Step 5: Generate combined count matrix using deprep.py
# -----------------------------
# -----------------------------

#-------- Create samples.txt automatically --------#

SAMPLES_FILE="$OUT_DIR/stringtie_assembly/samples.txt"

echo "Generating samples.txt for prepDE.py3..."
> "$SAMPLES_FILE"  # clear file if exists

for gtf in "$OUT_DIR/stringtie_assembly/"*.gtf; do
    base=$(basename "$gtf" .gtf)
    # Write: sample_name <tab> full_path_to_gtf
    echo -e "${base}\t${gtf}" >> "$SAMPLES_FILE"
done

echo "samples.txt created at: $SAMPLES_FILE"


#-------- run prepDE.py3 using this file:  --------#

python3  prepDE.py3 \
  -i $SAMPLES_FILE 

echo "gene and transcript count matrix were generated "
echo "gene and transcript count matrix were generated "

conda deactivate 

