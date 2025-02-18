Stickleback whole genome poolseq data from 27 Canada and Alaska populations. I would like to use mtDNA consensus sequences (extracted from the Poolseq data) for my final project. I’ve trimmed, mapped, checked coverage, converted to bam.files, deduplicates, marked Indel targets and realignment with GATK. 

1.Trim raw_data with bbduk 

#!/bin/bash

# Set input and output directories and parameters
INPUT_DIR="/work/cyu/poolseq/raw_data"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
ADAPTER_REF="/home/cyu/adapters.fa"
THREADS=48

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all paired R1 and R2 files
for R1_FILE in "$INPUT_DIR"/*_R1_001.fastq; do
    # Get the corresponding R2 file
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    
    # Extract the file prefix (e.g., 16_AMO_S16)
    PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001.fastq//')
    
    # Set output file paths
    OUT_R1="$OUTPUT_DIR/trimmed_R1_${PREFIX}.fastq"
    OUT_R2="$OUTPUT_DIR/trimmed_R2_${PREFIX}.fastq"
    LOG_FILE="$OUTPUT_DIR/${PREFIX}_bbduk_log.txt"
    
    # Run bbduk.sh
    bbduk.sh \
        in1="$R1_FILE" \
        in2="$R2_FILE" \
        out1="$OUT_R1" \
        out2="$OUT_R2" \
        ref="$ADAPTER_REF" \
        ktrim=rl \
        trimq=20 \
        minlength=25 \
        ftl=10 \
        tossbrokenreads=t \
        threads="$THREADS" > "$LOG_FILE" 2>&1
done

2. Map trimmed data to reference mt genome with bowtie2

#!/bin/bash
# Set input and output paths
TRIMMED_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
MAPPED_DIR="/work/cyu/PPalign_output/mapped"
REFERENCE="/work/cyu/chrM_reference"

# Create output directory (if it does not exist)
mkdir -p "$MAPPED_DIR"

# Loop through all R1 files and find the corresponding R2 files
for R1_FILE in "$TRIMMED_DIR"/trimmed_R1_*.fastq; do
    # Get the base name of the file (remove the path and prefix)
    BASENAME=$(basename "$R1_FILE" | sed 's/trimmed_R1_//; s/.fastq//')
    R2_FILE="$TRIMMED_DIR/trimmed_R2_$BASENAME.fastq"
    
    # Check if the R2 file exists
    if [[ -f "$R2_FILE" ]]; then
        # Output file paths
        SAM_FILE="$MAPPED_DIR/${BASENAME}.sam"
        LOG_FILE="$MAPPED_DIR/${BASENAME}_bowtie2.log"
        
        # Perform mapping and save logs
        echo "Mapping $BASENAME..."
        bowtie2 -x "$REFERENCE" \
            -1 "$R1_FILE" \
            -2 "$R2_FILE" \
            -p 48 \
            --very-sensitive-local \
            --no-mixed \
            --no-discordant \
            -X 2000 \
            -S "$SAM_FILE" > "$LOG_FILE" 2>&1
    else
        echo "Warning: R2 file for $BASENAME not found, skipping." | tee -a "$MAPPED_DIR/mapping_warnings.log"
    fi
done

echo "Mapping completed."

3. Convert mapped sam files to bam and index

#!/bin/bash

# Input and output paths
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"

# Loop through all SAM files
for SAM_FILE in "${INPUT_DIR}"/*.sam; do
    # Get the base name of the file (remove the path and extension)
    BASENAME=$(basename "$SAM_FILE" .sam)
    
    # Set the output file path
    BAM_FILE="${OUTPUT_DIR}/${BASENAME}_sorted.bam"
    
    echo "Processing $SAM_FILE..."
    
    # Convert to BAM format, filter unmapped reads, then sort
    samtools view -Sbu -F 4 "$SAM_FILE" | samtools sort -o "$BAM_FILE" -
    
    # Generate BAM index
    samtools index "$BAM_FILE"
    
    echo "Finished processing $SAM_FILE: Sorted BAM and index created."
done

echo "All SAM files have been processed."



4. Check the coverage/depth

#!/bin/bash

# Input and output directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"

# Coverage output directory
COVERAGE_DIR="/work/cyu/poolseq/PPalign_output/coverage"
mkdir -p "$COVERAGE_DIR"

# Process all sorted BAM files
for BAM_FILE in "${OUTPUT_DIR}"/*_sorted.bam; do
    # Get the basename (removing directory and extension)
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Coverage output file
    COVERAGE_FILE="${COVERAGE_DIR}/${BASENAME}_coverage.txt"
    
    echo "Calculating depth for $BAM_FILE..."
    
    # Calculate depth
    samtools depth "$BAM_FILE" > "$COVERAGE_FILE"
    
    echo "Coverage file created: $COVERAGE_FILE"
done

echo "All sorted BAM files have been processed for depth calculation."

5. Piscard mark duplicates
#!/bin/bash

# Directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/marked_duplicates"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all sorted BAM files
for BAM_FILE in "${INPUT_DIR}"/*_sorted.bam; do
    # Get the base name without extension
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    
    # Output BAM file path for the marked duplicates
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_marked_duplicates.bam"
    
    # Run Picard MarkDuplicates
    echo "Processing $BAM_FILE..."
    picard MarkDuplicates \
        I="$BAM_FILE" \
        O="$OUTPUT_BAM" \
        M="${OUTPUT_DIR}/${BASENAME}_metrics.txt" \
        REMOVE_DUPLICATES=true
    
    echo "Finished processing $BAM_FILE: Duplicates marked."
done

echo "All BAM files have been processed


6. Piscard add group name
#!/bin/bash

# Directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/marked_duplicates"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all marked BAM files (after dedup BAM 文件)
for BAM_FILE in "${INPUT_DIR}"/*_marked_duplicates.bam; do
    # Get the base name without extension
    BASENAME=$(basename "$BAM_FILE" _marked_duplicates.bam)
    
    # Output BAM file path
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_rg.bam"
    
    # Run Picard AddOrReplaceReadGroups to add read groups
    echo "Processing $BAM_FILE..."
    java -Xmx16g -jar /home/cyu/picard/picard.jar AddOrReplaceReadGroups \
        -I "$BAM_FILE" \
        -O "$OUTPUT_BAM" \
        -RGID "$BASENAME" \
        -RGLB "$BASENAME" \
        -RGPL ILLUMINA \
        -RGPU "${BASENAME}_Unit1" \
        -RGSM "$BASENAME"
    
    echo "Finished processing $BAM_FILE: Read groups added."
done

echo "All BAM files have been processed."


7. Index group name added files

#!/bin/bash

# Input directory
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"  

# Loop through all _rg.bam files
for BAM_FILE in "${INPUT_DIR}"/*_rg.bam; do
    # Get the sample name (remove the extension)
    SAMPLE_NAME=$(basename "$BAM_FILE" _rg.bam)
    
    # Index the BAM file
    echo "Indexing BAM file for ${SAMPLE_NAME}..."
    samtools index "$BAM_FILE"  # Create the BAM file index
    
    echo "Finished indexing ${SAMPLE_NAME}"
done

echo "All BAM files have been indexed."



8. GATK3 generate Indel targets and run realignment
#!/bin/bash

# Input and output directories
INPUT_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"  
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/realigned"
REFERENCE="/work/cyu/chrM_reference.fa"  # Reference genome
GATK_JAR="/home/cyu/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"  # Path to GATK JAR

# Create output directory (if it does not exist)
mkdir -p "$OUTPUT_DIR"

# Loop through all _rg.bam files
for BAM_FILE in "${INPUT_DIR}"/*_rg.bam; do
    # Get the sample name (remove the extension)
    SAMPLE_NAME=$(basename "$BAM_FILE" _rg.bam)
    
    # Define output file paths
    INTERVAL_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.intervals"
    REALIGNED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_realigned.bam"
    LOG_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_realignment_log.txt"
    
    # Step 1: Run RealignerTargetCreator to generate the InDel target list
    echo "Running RealignerTargetCreator for ${SAMPLE_NAME}..."
    java -Xmx8g -jar "$GATK_JAR" \
        -T RealignerTargetCreator \
        -R "$REFERENCE" \
        -I "$BAM_FILE" \
        -o "$INTERVAL_FILE" > "$LOG_FILE" 2>&1
    
    # Step 2: Run IndelRealigner for realignment
    echo "Running IndelRealigner for ${SAMPLE_NAME}..."
    java -Xmx8g -jar "$GATK_JAR" \
        -T IndelRealigner \
        -R "$REFERENCE" \
        -I "$BAM_FILE" \
        -targetIntervals "$INTERVAL_FILE" \
        --maxReadsForRealignment 1000000 \
        -o "$REALIGNED_BAM" >> "$LOG_FILE" 2>&1
    
    echo "Finished processing ${SAMPLE_NAME}. Logs saved to $LOG_FILE"
done

echo "All BAM files have been processed. Logs are saved in $OUTPUT_DIR."

9. Create bam.list
ls -1 /work/cyu/poolseq/PPalign_output/realigned/*.bam > /work/cyu/poolseq/PPalign_output/realigned/bamlist.txt

# merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20
samtools mpileup -B \
    -f /work/cyu/chrM_reference.fa \
    -b bamlist.txt \
    -q 20 \
    -Q 20 \
    | gzip > fish.mpileup.gz

#call SNPs with PoolSNP



10. Run PoolSNP for SNP calling

bash /home/cyu/.conda/envs/poolseq_env/bin/PoolSNP.sh \
mpileup=/work/cyu/poolseq/PPalign_output/realigned/fish.mpileup.gz \
reference=/work/cyu/chrM_reference.fa \
names=10_THE,11_JOE,12_BEA,13_MUC,14_PYE,16_AMO,17_SAY,18_GOS,19_ROB,1_FG,20_BOOT,21_ECHO,22_FRED,23_LAW,24_PACH,25_RS,26_SC,27_LB,28_CH,2_LG,3_SR,4_SL,5_TL,6_WB,7_WT,8_WK,9_SWA \
max-cov=0.999 \
min-cov=10 \
min-count=3 \
min-freq=0.001 \
miss-frac=0.1 \
jobs=1 \
BS=1 \
output=/work/cyu/poolseq/PPalign_output/realigned/snp \



11. Detect indels using https://github.com/capoony/DrosEU_pipeline.git scripts
python DetectIndels.py \
--mpileup /work/cyu/poolseq/PPalign_output/realigned/fish.mpileup.gz \
--minimum-count 5 \
--mask 5 \
| gzip > /work/cyu/poolseq/PPalign_output/realigned/InDel-positions_20.txt.gz



12. Filter indel using https://github.com/capoony/DrosEU_pipeline.git scripts
python FilterPosFromVCF.py \
--indel /work/cyu/poolseq/PPalign_output/realigned/InDel-positions_20.txt.gz \
--vcf /work/cyu/poolseq/PPalign_output/realigned/snp.vcf.gz \
| gzip > /work/cyu/poolseq/PPalign_output/realigned/SNPs_clean.vcf.gz


13. Navigate to the vcf.file directory
cd /work/cyu/poolseq/PPalign_output/realigned/

14. Decompress the VCF file and output it as an uncompressed VCF
gunzip -c /work/cyu/poolseq/PPalign_output/realigned/SNPs_clean.vcf.gz > SNPs_clean.vcf

15. Recompress the VCF file using bgzip (required for indexing with tabix)
bgzip -c SNPs_clean.vcf > /work/cyu/poolseq/PPalign_output/realigned/SNPs_clean.vcf.gz

16. Index the compressed VCF file using tabix
tabix -p vcf /work/cyu/poolseq/PPalign_output/realigned/SNPs_clean.vcf.gz

17. Bcftools generated consensus sequences

for sample in 10_THE 11_JOE 12_BEA 13_MUC 14_PYE 16_AMO 17_SAY 18_GOS 19_ROB 1_FG 20_BOOT 21_ECHO 22_FRED 23_LAW 24_PACH 25_RS 26_SC 27_LB 28_CH 2_LG 3_SR 4_SL 5_TL 6_WB 7_WT 8_WK 9_SWA
do
  echo "Processing sample: ${sample}"
  bcftools consensus -f /work/cyu/chrM_reference.fa -s ${sample} /work/cyu/poolseq/PPalign_output/realigned/SNPs_clean.vcf.gz > /work/cyu/poolseq/PPalign_output/realigned/consensus_sequences/${sample}.fa
done


18. Revise header names to population names since the headers should all be chrM
#!/bin/bash

# Define input directory
FASTA_DIR="/work/cyu/poolseq/PPalign_output/realigned/consensus_sequences"

# Navigate to the directory
cd "$FASTA_DIR" || exit

# Process all .fa files
for file in *.fa; do
    # Extract population name (remove numeric prefix and underscore)
    pop_name=$(echo "$file" | sed -E 's/^[0-9]+_//; s/.fa$//')

    # Modify the FASTA header
    sed -i "1s/>.*/>${pop_name}/" "$file"

    echo "✔ Modified FASTA header for $file to: >${pop_name}"
done

echo "All FASTA headers have been updated!"



