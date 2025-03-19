#This script uses IQ-TREE for maximum likelihood (ML) phylogenetic inference.
#   IQ-TREE is a fast and efficient ML method that infers phylogenetic trees 
#   based on aligned sequence data. In this example, we use the GTR+G model.
#
#   The GTR (General Time Reversible) model assumes that substitution rates 
#   between nucleotides are reversible, and the +G indicates that rate 
#   variation among sites is modeled using a gamma distribution.
#
# Assumptions:
#   - The input alignment file (aligned_mt.fasta) is correctly formatted and aligned.
#   - The GTR+G model adequately represents the evolutionary process of the dataset.
#   - Sites are assumed to evolve independently under homogeneous conditions.
#
# Limitations:
#   - Errors in the alignment or model misspecification can affect the results.
#   - The assumption of site-independence might not hold for all datasets.
#   - Larger or more complex datasets will require more computational resources.

iqtree -s aligned_mt.fasta -m GTR+G -nt AUTO