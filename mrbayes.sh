#This script is for NEXUS format converting and run MrBayes
#Software-MrBayes Bayesian Phylogenetic Inference. MrBayes is a software tool that performs Bayesian inference of phylogeny using Markov Chain Monte Carlo (MCMC) methods to estimate the posterior probability distribution of trees. It allows the use of flexible models of sequence evolution, such as GTR+G (General Time Reversible with gamma-distributed rate variation among sites), and provides posterior probabilities for each clade, offering a measure of confidence in the inferred topology.
#Assumptions: Input sequences are properly aligned. The chosen substitution model (e.g., GTR+G) accurately reflects the evolutionary process. Sites evolve independently and identically under the same model. The MCMC chain converges given enough generations. Prior distributions are appropriate and do not overly influence the posterior.
#Strengths: Provides posterior probabilities to quantify uncertainty in tree topology. Supports complex models, including partitioning and relaxed clock models. Can incorporate prior biological knowledge via priors.
#Limitations: Computationally intensive and slower than maximum likelihood methods. Requires careful convergence diagnostics (e.g., PSRF, ESS). Results may be sensitive to prior choice and model specification.

# 1. Convert FASTA format to NEXUS format
# This script converts a multiple sequence alignment in FASTA format
# into NEXUS format, which is compatible with PopART, BEAST, MrBayes, and PAUP.

from Bio import SeqIO

def fasta_to_nexus(input_fasta, output_nexus):
    """
    Convert a FASTA file to NEXUS format for PopART, BEAST, MrBayes, and PAUP.
    Ensures all sequences have the same length and adds a MrBayes block at the end.
    """

    # Read FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Ensure all sequences have the same length
    seq_lengths = {len(seq.seq) for seq in sequences}
    if len(seq_lengths) > 1:
        raise ValueError("Error: Sequences have different lengths. Alignment is required.")

    n_tax = len(sequences)              # Number of taxa
    n_char = len(sequences[0].seq)     # Sequence length

    # Write to NEXUS file
    with open(output_nexus, "w") as nexus:
        # NEXUS header
        nexus.write("#NEXUS\n\n")

        # TAXA block
        nexus.write("BEGIN TAXA;\n")
        nexus.write(f"  DIMENSIONS NTAX={n_tax};\n")
        nexus.write("  TAXLABELS\n")
        for seq in sequences:
            nexus.write(f"    {seq.id}\n")
        nexus.write("  ;\nEND;\n\n")

        # CHARACTERS block
        nexus.write("BEGIN CHARACTERS;\n")
        nexus.write(f"  DIMENSIONS NCHAR={n_char};\n")
        nexus.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            nexus.write(f"    {seq.id} {str(seq.seq)}\n")
        nexus.write("  ;\nEND;\n\n")

        # MrBayes block
        nexus.write("BEGIN MRBAYES;\n")
        nexus.write("  set autoclose=yes nowarn=yes;\n")
        nexus.write("  lset nst=6 rates=gamma;\n")
        nexus.write("  mcmc ngen=1000000 samplefreq=100 printfreq=100 nchains=4 savebrlens=yes diagnfreq=1000;\n")
        nexus.write("  sumt burnin=2500;\n")
        nexus.write("  sump burnin=2500;\n")
        nexus.write("END;\n")

# Set input & output file paths
input_fasta = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt.fasta"
output_nexus = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_mb.nex"

# Run the conversion
fasta_to_nexus(input_fasta, output_nexus)

print(f"✅ Nexus file with MrBayes block created: {output_nexus}")

#2. run mrbayes
mb /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_mb.nex


print(f"✅ Nexus file created: {output_nexus}")
