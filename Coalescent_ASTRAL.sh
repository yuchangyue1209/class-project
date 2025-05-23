# STEP 16: COALESCENT SPECIES TREE INFERENCE (ASTRAL)
# Software: ASTRAL v5.7.8 (https://github.com/smirarab/ASTRAL)
# Installation:
#   git clone https://github.com/smirarab/ASTRAL.git
#   Requires: Java 8+ (e.g., conda install -c conda-forge openjdk)

# Description:
#   ASTRAL is a statistically consistent method for inferring species trees from unrooted gene trees using the multi-species coalescent (MSC) model.
#   It infers the species tree that agrees with the largest number of quartet trees induced from the input gene trees.

# Strengths:
#   - Handles incomplete lineage sorting (ILS) well
#   - Scalable and fast, no need for sequence alignment inputs
#   - Ideal for summarizing across many gene trees

# Weaknesses:
#   - Sensitive to gene tree errors
#   - Ignores gene duplication/loss and horizontal transfer
#   - Only considers topologies (not branch lengths)

# Assumptions:
#   - Input gene trees are accurate and independent
#   - Discordance is due to ILS
#   - Genes are unlinked and non-recombining

# User Choices:
#   - Used IQ-TREE to infer gene trees from aligned CDS regions (13 PCGs from mitochondrial genome)
#   - Ran ASTRAL on all .treefile outputs from these gene trees



# 1. Extract CDS Coordinates from GFF
# Extract all CDS regions from a GFF file and save as BED format
grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
awk '{print $1, $4, $5}' OFS='\t' > /work/cyu/poolseq/PPalign_output/ann/pcg.bed

#2.Extract All 13 PCG Sequences from Full mtDNA Alignment
#!/bin/bash

# Input paths
fasta="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.fasta"
bed="/work/cyu/poolseq/PPalign_output/ann/pcg.bed"
output_dir="/work/cyu/poolseq/PPalign_output/consensus/per_gene_aligned_from_full_fasta"
mkdir -p "$output_dir"

# Number and format each gene region line as gene001, gene002, ...
awk '{print $2"-"$3}' "$bed" | nl -n rz -w 3 > tmp_gene_regions.txt

# Loop through each region and extract sequences from full alignment
while read geneid region; do
    start=$(echo "$region" | cut -d- -f1)
    end=$(echo "$region" | cut -d- -f2)
    gene="gene${geneid}"
    outfile="$output_dir/${gene}_all_pops.fasta"
    > "$outfile"

    awk -v gene="$gene" -v start="$start" -v end="$end" '
    /^>/ {
        if (seq) {
            subseq = substr(seq, start, end - start + 1)
            print ">" gene "_" name >> out
            print subseq >> out
        }
        name = substr($0, 2)
        seq = ""
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq) {
            subseq = substr(seq, start, end - start + 1)
            print ">" gene "_" name >> out
            print subseq >> out
        }
    }
    ' out="$outfile" "$fasta"

    echo "✅ Extracted: $gene ($start-$end)"
done < tmp_gene_regions.txt

rm tmp_gene_regions.txt

#3. Align and Build Gene Trees with MAFFT + IQ-TREE
#!/bin/bash

# Input and output
fasta="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_overlap.fasta"
bed="/work/cyu/poolseq/PPalign_output/ann/pcg.bed"
outdir="/work/cyu/poolseq/PPalign_output/consensus/per_gene_aligned_from_full_fasta"
mkdir -p "$outdir"

awk '{print $2"-"$3}' "$bed" | nl -n rz -w 3 > tmp_gene_regions.txt

while read geneid region; do
    start=$(echo "$region" | cut -d- -f1)
    end=$(echo "$region" | cut -d- -f2)
    gene="gene${geneid}"
    fasta_out="$outdir/${gene}_all_pops.fasta"
    aln_out="$outdir/${gene}_aligned.fasta"

    > "$fasta_out"

    # Extract
    awk -v gene="$gene" -v start="$start" -v end="$end" '
    /^>/ {
        if (seq) {
            subseq = substr(seq, start, end - start + 1)
            print ">" gene "_" name >> out
            print subseq >> out
        }
        name = substr($0, 2)
        seq = ""
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq) {
            subseq = substr(seq, start, end - start + 1)
            print ">" gene "_" name >> out
            print subseq >> out
        }
    }
    ' out="$fasta_out" "$fasta"

    echo "🔄 MAFFT aligning: $gene"
    mafft --auto "$fasta_out" > "$aln_out"

    echo "IQ-TREE inferring: $gene"
    iqtree2 -s "$aln_out" -m GTR+G -bb 1000 -nt AUTO &> "$outdir/${gene}_iqtree.log"

    echo "✅ Done: $gene"
done < tmp_gene_regions.txt

rm tmp_gene_regions.txt


#4. install and run ASTRAL to infer species tree
git clone https://github.com/smirarab/ASTRAL.git

#!/bin/bash

# Set input/output paths
treedir="/work/cyu/poolseq/PPalign_output/consensus/per_gene_aligned_from_full_fasta"
astral_jar="$HOME/ASTRAL/Astral/astral.5.7.8.jar"
outdir="$treedir/ASTRAL_output"
mkdir -p "$outdir"

# Merge all gene trees into a single file
cat "$treedir"/*.treefile > "$outdir/all_gene_trees.tre"
echo "Merged all trees: $outdir/all_gene_trees.tre"

# Run ASTRAL
echo "Running ASTRAL..."
java -Xmx8G -jar "$astral_jar" -i "$outdir/all_gene_trees.tre" -o "$outdir/astral_species_tree.tre"

echo "ASTRAL species tree completed: $outdir/astral_species_tree.tre"
