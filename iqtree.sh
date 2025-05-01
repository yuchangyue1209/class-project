# mtDNA Variant Calling and Consensus Construction from Pool-seq Data
# Description: Reproducible pipeline from raw Pool-seq FASTQ to consensus mtDNA FASTA sequences

#############################################
# STEP 13: PHYLOGENETIC TREE CONSTRUCTION (IQ-TREE)
#############################################
# Software: IQ-TREE (http://www.iqtree.org/)
# Installation:
#   conda install -c bioconda iqtree -y
#   (or download from https://github.com/iqtree/iqtree2)
#
# Description:
#   IQ-TREE is a fast and efficient maximum likelihood (ML) method for phylogenetic tree inference
#   based on aligned sequence data. Here, the GTR+G model is used.
#
# Model:
#   GTR (General Time Reversible) allows for different substitution rates between all pairs of nucleotides.
#   +G adds a gamma distribution to account for rate variation across sites.
#
# Assumptions:
#   - The input alignment file is correctly formatted (e.g., FASTA) and sequences are aligned.
#   - Substitution model adequately reflects biological reality.
#   - Sites evolve independently and under homogeneous conditions.
#
# Strengths:
#   - Fast and accurate, model selection included (ModelFinder)
#   - Robust bootstrapping options (e.g., ultrafast bootstraps)
#
# Limitations:
#   - Sensitive to alignment quality and model choice.
#   - Assumes no recombination and site independence.
#
# Example command:
14. IQtree 
iqtree -s aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
iqtree -s mc_aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
# Output: IQ-TREE generates .treefile (ML tree), .log (log file), and .iqtree (summary of the run)
