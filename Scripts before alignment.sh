# mtDNA Variant Calling and Consensus Construction from Pool-seq Data
# Description: Reproducible pipeline from raw Pool-seq FASTQ to consensus mtDNA FASTA sequences

#############################################
# STEP 0: INSTALL REQUIRED SOFTWARE & TOOLS
#############################################
# Requires: Conda (https://docs.conda.io/en/latest/miniconda.html)
# Install all software tools with conda:

# Create and activate environment
conda create -n poolseq_env python=3.10 -y
conda activate poolseq_env

# Install all necessary tools
conda install -c bioconda \
    bbmap=38.90 \
    fastqc \
    bowtie2=2.5.4 \
    bwa \
    samtools=1.17 \
    bcftools=1.17 \
    vcftools \
    gatk4=4.6.1.0 \
    picard -y

# Optional: Java 17 may be needed for Picard if not included
# e.g., conda install -c conda-forge openjdk=17

#############################################
# STEP 1: TRIM RAW DATA WITH BBDUK (bbmap 38.9)
#############################################
# Tool: bbduk.sh (from BBMap suite)
# Purpose: Adapter trimming, quality filtering

INPUT_DIR="/work/cyu/poolseq/raw_data"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
ADAPTER_REF="/home/cyu/adapters.fa"
THREADS=48
mkdir -p "$OUTPUT_DIR"

for R1_FILE in "$INPUT_DIR"/*_R1_001.fastq; do
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001.fastq//')
    OUT_R1="$OUTPUT_DIR/trimmed_R1_${PREFIX}.fastq"
    OUT_R2="$OUTPUT_DIR/trimmed_R2_${PREFIX}.fastq"
    LOG_FILE="$OUTPUT_DIR/${PREFIX}_bbduk_log.txt"

    bbduk.sh \
        in1="$R1_FILE" \
        in2="$R2_FILE" \
        out1="$OUT_R1" \
        out2="$OUT_R2" \
        ref="$ADAPTER_REF" \
        ktrim=rl trimq=20 minlength=25 ftl=10 \
        tossbrokenreads=t threads="$THREADS" > "$LOG_FILE" 2>&1
done

mkdir -p /work/cyu/poolseq/PPalign_output/quality_after_trim/
for r1 in $OUTPUT_DIR/trimmed_R1_*.fastq; do
    r2="${r1/R1/R2}"
    fastqc "$r1" "$r2" -o /work/cyu/poolseq/PPalign_output/quality_after_trim/
done

#############################################
# STEP 2: MAP TO MITOCHONDRIAL GENOME (bowtie2 2.5.4)
#############################################
# Tool: bowtie2
# Purpose: Align trimmed reads to mitochondrial reference genome

TRIMMED_DIR="$OUTPUT_DIR"
MAPPED_DIR="/work/cyu/poolseq/PPalign_output/mapped"
REFERENCE="/work/cyu/chrM_index"
mkdir -p "$MAPPED_DIR"

for R1_FILE in "$TRIMMED_DIR"/trimmed_R1_*.fastq; do
    BASENAME=$(basename "$R1_FILE" | sed 's/trimmed_R1_//; s/.fastq//')
    R2_FILE="$TRIMMED_DIR/trimmed_R2_$BASENAME.fastq"
    [ -f "$R2_FILE" ] || continue

    SAM_FILE="$MAPPED_DIR/${BASENAME}.sam"
    LOG_FILE="$MAPPED_DIR/${BASENAME}_bowtie2.log"

    bowtie2 -x "$REFERENCE" -1 "$R1_FILE" -2 "$R2_FILE" \
        -p 48 --very-sensitive-local --no-mixed --no-discordant -X 2000 \
        -S "$SAM_FILE" > "$LOG_FILE" 2>&1
done

#############################################
# STEP 3: CONVERT SAM TO SORTED BAM (samtools 1.17)
#############################################

for SAM_FILE in "$MAPPED_DIR"/*.sam; do
    BASENAME=$(basename "$SAM_FILE" .sam)
    BAM_FILE="$MAPPED_DIR/${BASENAME}.bam"
    SORTED_BAM_FILE="$MAPPED_DIR/${BASENAME}_sorted.bam"
    samtools view -b -q 20 "$SAM_FILE" > "$BAM_FILE"
    samtools sort -o "$SORTED_BAM_FILE" "$BAM_FILE"
done

#############################################
# STEP 4: CALCULATE DEPTH (samtools depth)
#############################################

DEPTH_OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/depth"
mkdir -p "$DEPTH_OUTPUT_DIR"

for BAM_FILE in "$MAPPED_DIR"/*_sorted.bam; do
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    samtools depth "$BAM_FILE" > "$DEPTH_OUTPUT_DIR/${BASENAME}_depth.txt"
done

#############################################
# STEP 5: OPTIONAL DEDUPLICATION (picard MarkDuplicates)
#############################################
# Skip for mtDNA due to high copy number and unique mapping patterns

#############################################
# STEP 6: ADD READ GROUPS (picard AddOrReplaceReadGroups)
#############################################

RG_DIR="/work/cyu/poolseq/PPalign_output/mapped_with_rg"
mkdir -p "$RG_DIR"

for BAM_FILE in "$MAPPED_DIR"/*_sorted.bam; do
    BASENAME=$(basename "$BAM_FILE" _sorted.bam)
    OUTPUT_BAM="$RG_DIR/${BASENAME}_rg.bam"
    java -Xmx16g -jar "/home/cyu/picard/picard.jar" AddOrReplaceReadGroups \
        -I "$BAM_FILE" -O "$OUTPUT_BAM" \
        -RGID "$BASENAME" -RGLB "$BASENAME" -RGPL ILLUMINA \
        -RGPU "${BASENAME}_Unit1" -RGSM "$BASENAME"
done

#############################################
# STEP 7: INDEX BAM FILES (samtools index)
#############################################

for BAM_FILE in "$RG_DIR"/*_rg.bam; do
    samtools index "$BAM_FILE"
done

#############################################
# STEP 8: VARIANT CALLING WITH GATK (HaplotypeCaller)
#############################################

POOL_INFO="/work/cyu/poolseq/pool_info.txt"
REFERENCE="/work/cyu/sequence.fasta"
GATK_OUT="/work/cyu/poolseq/PPalign_output/gatk4_vcf"
mkdir -p "$GATK_OUT"

tail -n +2 "$POOL_INFO" | while read -r SAMPLE SIZE; do
    CLEAN_SIZE=$(echo "$SIZE" | tr -cd '0-9')
    BAM_FILE="$RG_DIR/${SAMPLE}_rg.bam"
    gatk HaplotypeCaller \
        -R "$REFERENCE" -I "$BAM_FILE" \
        -O "$GATK_OUT/${SAMPLE}.g.vcf.gz" -ERC GVCF \
        --ploidy "$CLEAN_SIZE" --minimum-mapping-quality 20 -mbq 13
    done

#############################################
# STEP 9: VARIANT CALLING WITH BCFTOOLS (bcftools mpileup/call)
#############################################

MTDNA_ID="MH205729.1"
VCF_DIR="/work/cyu/poolseq/PPalign_output/vcf"
BAM_DIR="/work/cyu/poolseq/PPalign_output/mtDNA_bam"
mkdir -p "$BAM_DIR" "$VCF_DIR"

for BAM in "$MAPPED_DIR"/*_sorted.bam; do
    SAMPLE=$(basename "$BAM" _sorted.bam)
    samtools view -b -q 20 -o "$BAM_DIR/${SAMPLE}_mtDNA.bam" "$BAM" "$MTDNA_ID"
    samtools index "$BAM_DIR/${SAMPLE}_mtDNA.bam"
    bcftools mpileup -Ou -f "$REFERENCE" "$BAM_DIR/${SAMPLE}_mtDNA.bam" | \
    bcftools call --ploidy 1 -mv -Oz -o "$VCF_DIR/${SAMPLE}_mtDNA_raw.vcf.gz"
    bcftools norm -m -any -f "$REFERENCE" -Oz \
        -o "$VCF_DIR/${SAMPLE}_mtDNA.vcf.gz" "$VCF_DIR/${SAMPLE}_mtDNA_raw.vcf.gz"
    bcftools index "$VCF_DIR/${SAMPLE}_mtDNA.vcf.gz"
done

#############################################
# STEP 10: INTERSECT GATK AND BCFTOOLS VARIANTS (bcftools isec)
#############################################

COMPARE_OUT="/work/cyu/poolseq/PPalign_output/compare_results"
mkdir -p "$COMPARE_OUT"

for BCF_VCF in "$VCF_DIR"/*_mtDNA.vcf.gz; do
    SAMPLE=$(basename "$BCF_VCF" _mtDNA.vcf.gz)
    GATK_VCF="$GATK_OUT/single_sample_vcf/${SAMPLE}_final.vcf.gz"
    OUTDIR="$COMPARE_OUT/${SAMPLE}_isec"
    mkdir -p "$OUTDIR"
    bcftools isec "$BCF_VCF" "$GATK_VCF" -p "$OUTDIR" -Oz
done

#############################################
# STEP 11: FILTER OVERLAPPING SNPS (vcftools)
#############################################

OVERLAP_VCF="/work/cyu/poolseq/PPalign_output/overlap.vcf"
mkdir -p "$OVERLAP_VCF"

for DIR in "$COMPARE_OUT"/*_isec; do
    SAMPLE=$(basename "$DIR" _isec)
    INFILE="$DIR/0002.vcf.gz"
    vcftools --gzvcf "$INFILE" --minDP 10 --maxDP 5000 --max-missing 1 --minQ 30 \
      --recode --recode-INFO-all --out "$OVERLAP_VCF/${SAMPLE}_filtered"
    bgzip -c "$OVERLAP_VCF/${SAMPLE}_filtered.recode.vcf" > "$OVERLAP_VCF/${SAMPLE}_overlap.vcf.gz"
    bcftools index "$OVERLAP_VCF/${SAMPLE}_overlap.vcf.gz"
done

#############################################
# STEP 12: CALL CONSENSUS FASTA SEQUENCES (bcftools consensus)
#############################################

CONSENSUS_DIR="$OVERLAP_VCF/consensus"
mkdir -p "$CONSENSUS_DIR"

for VCF in "$OVERLAP_VCF"/*_overlap.vcf.gz; do
    SAMPLE=$(basename "$VCF" _overlap.vcf.gz)
    bcftools consensus -f "$REFERENCE" "$VCF" > "$CONSENSUS_DIR/${SAMPLE}_consensus.fasta"
    POP_NAME=$(echo "$SAMPLE" | cut -d'_' -f2)
    awk -v name="$POP_NAME" '/^>/{print ">" name; next} {print}' "$CONSENSUS_DIR/${SAMPLE}_consensus.fasta" > tmp && mv tmp "$CONSENSUS_DIR/${SAMPLE}_consensus.fasta"
done

cat "$CONSENSUS_DIR"/*.fasta > "$CONSENSUS_DIR/all_mt_overlap.fasta"
echo "Pipeline complete. Final consensus FASTA: $CONSENSUS_DIR/all_mt_overlap.fasta"
