# STEP 15: CONVERT FASTA TO NEXUS & RUN MRBAYES
#############################################
# Software: MrBayes (http://nbisweden.github.io/MrBayes/)
# Installation:
#   conda install -c bioconda mrbayes -y
#   pip install biopython  # Required for conversion script

# Description:
#   MrBayes performs Bayesian phylogenetic inference using MCMC.
#   This step converts aligned FASTA to NEXUS and runs MrBayes with a standard GTR+G model.

# Python conversion script (uses BioPython)
# Converts FASTA alignment into NEXUS format with embedded MrBayes block

# Input: /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt.fasta
# Output: /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_mb.nex

# MrBayes block includes:
# - GTR+G model (nst=6, rates=gamma)
# - 1,000,000 generations, 4 chains, sampled every 100 generations
# - Burn-in of 2,500 samples (25%)

# Assumptions:
#   - Sequences are properly aligned
#   - Substitution model matches evolutionary process
#   - Priors and MCMC settings are appropriate

# Strengths:
#   - Supports model complexity, uncertainty estimation
#   - Posterior probabilities give branch support

# Limitations:
#   - Slower than ML methods
#   - Requires convergence checks (ESS, PSRF)

# Example Python code for conversion and execution:

from Bio import SeqIO

def fasta_to_nexus(input_fasta, output_nexus):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    seq_lengths = {len(seq.seq) for seq in sequences}
    if len(seq_lengths) > 1:
        raise ValueError("Sequences have different lengths. Alignment is required.")
    n_tax = len(sequences)
    n_char = len(sequences[0].seq)

    with open(output_nexus, "w") as nexus:
        nexus.write("#NEXUS\n\n")
        nexus.write("BEGIN TAXA;\n")
        nexus.write(f"  DIMENSIONS NTAX={n_tax};\n")
        nexus.write("  TAXLABELS\n")
        for seq in sequences:
            nexus.write(f"    {seq.id}\n")
        nexus.write("  ;\nEND;\n\n")

        nexus.write("BEGIN CHARACTERS;\n")
        nexus.write(f"  DIMENSIONS NCHAR={n_char};\n")
        nexus.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus.write("  MATRIX\n")
        for seq in sequences:
            nexus.write(f"    {seq.id} {str(seq.seq)}\n")
        nexus.write("  ;\nEND;\n\n")

        nexus.write("BEGIN MRBAYES;\n")
        nexus.write("  set autoclose=yes nowarn=yes;\n")
        nexus.write("  lset nst=6 rates=gamma;\n")
        nexus.write("  mcmc ngen=1000000 samplefreq=100 printfreq=100 nchains=4 savebrlens=yes diagnfreq=1000;\n")
        nexus.write("  sumt burnin=2500;\n")
        nexus.write("  sump burnin=2500;\n")
        nexus.write("END;\n")

# Call the function
input_fasta = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt.fasta"
output_nexus = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_mb.nex"
fasta_to_nexus(input_fasta, output_nexus)


# run mrbayes
mb /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_mb.nex


print(f"✅ Nexus file created: {output_nexus}")
