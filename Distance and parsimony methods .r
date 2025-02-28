# Phylogenetic Tree Estimation using NJ and MP methods
library(ape)
library(phangorn)
library(adegenet)
# 1. Set Working Directory & Load Data
setwd("/users/changyueyu/desktop/")  # Update the path if needed

# Read the aligned mitochondrial sequences (FASTA format)
aligned_mt <- read.dna("aligned_mt.fasta", format = "fasta")

# 2. Distance-Based Phylogenetic Tree (NJ)


# Compute genetic distances using the Tamura-Nei 1993 (TN93) model
D <- dist.dna(aligned_mt, model = "TN93")

# Construct the Neighbor-Joining (NJ) tree
nj_tree <- nj(D)

# Ladderize for better visualization
nj_tree <- ladderize(nj_tree)

# Plot the NJ tree
plot(nj_tree, cex = 0.6)
title("Neighbor-Joining Tree")

# Save NJ tree in Newick format
write.tree(nj_tree, file = "nj_tree.newick")

# -------- Description of NJ Algorithm -------- #
# A distance-based method that constructs a tree by clustering taxa based on genetic distances. It is widely used for large datasets due to its speed and efficiency.
# Strengths: Fast and computationally efficient, especially for large datasets.
# Weaknesses: 1. Does not account for multiple substitutions at the same site (homoplasy) 2. Can be affected by long-branch attraction. 3.  Accuracy depends heavily on the distance model used.
# The genetic distance between taxa is additive (i.e., total distance along a path is the sum of distances between nodes).
# Evolutionary rates are constant across lineages (not always true).	- Choice of distance model (e.g., Tamura-Nei 1993 (TN93), Jukes-Cantor, Kimura 2-parameter).

# --------------------------------------
# 3. Parsimony-Based Phylogenetic Tree (MP)

# Convert DNA alignment into a phyDat object for parsimony analysis
phy_data <- as.phyDat(aligned_mt)  # FIXED: Removed incorrect 'type = "DNA"' argument

# Generate a starting tree using NJ
start_tree <- nj(D)

# Compute parsimony score of the starting tree
initial_parsimony_score <- parsimony(start_tree, phy_data)
cat("Initial parsimony score:", initial_parsimony_score, "\n")

# Optimize the tree using Maximum Parsimony (MP)
mp_tree <- optim.parsimony(start_tree, phy_data)

# Plot the MP tree
plot(mp_tree, cex = 0.6)
title("Maximum Parsimony Tree")

# Save MP tree in Newick format
write.tree(mp_tree, file = "mp_tree.newick")

# -------- Description of MP Algorithm -------- #
# - Maximum Parsimony (MP) is a character-based method that searches for the tree requiring the fewest evolutionary changes. It assumes that mutations are rare and prefers the simplest explanation 
# Conceptually simple and does not require a substitution model. Can provide accurate trees when mutation rates are low. Suitable for small datasets with well-aligned sequences.	
# Weakness: Can be misleading due to long-branch attraction (groups taxa with many shared changes, even if they are not closely related).  Computationally expensive for large datasets. Does not account for different mutation rates across sites. The simplest explanation (fewest mutations) is the most likely (Occam’s Razor).
# All nucleotide substitutions are equally probable.
# Does not model rate variation across sites.	- Choice of starting tree (often NJ).
# Whether to perform tree search optimizations like NNI (Nearest Neighbor Interchange).
# --------------------------------------------- #

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")
