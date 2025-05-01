#The alignment tools are muscle and mafft, before align the sequences, we need merge the sequences into one file.
cat /work/cyu/poolseq/PPalign_output/direct_consensus/*.fasta > all_mt.fasta

#Software	MUSCLE (Multiple Sequence Comparison by Log-Expectation)
#Description	Iterative alignment algorithm that improves alignment by refining progressive alignment trees based on pairwise k-mer distances and log-expectation scoring.
#Strengths	High accuracy for moderate-sized datasets; good speed; widely used; requires no external libraries.
#Weaknesses	Slower and less scalable than MAFFT for large datasets (>1000 seqs); fewer optimization strategies for long gapped regions.
#Assumptions	Assumes global homology across the whole alignment; substitution patterns follow a roughly consistent model.
#User choices	- -maxiters: controls refinement rounds (default is 16) - -diags: faster for closely related sequences- -quiet: suppresses standard output - Input should be unaligned FASTA
#Muscle install and rub

conda install muscle
muscle -in all_mt_overlap.fasta -out mc_aligned_mt_overlap.fasta

#Software	MAFFT 
#Description	Fast and flexible multiple sequence alignment tool using FFT-accelerated alignment. Supports multiple algorithms optimized for dataset size and divergence.
#Strengths	Excellent speed–accuracy tradeoff; highly scalable; many alignment strategies (--auto, L-INS-i, FFT-NS-2, etc.); robust for both short and long sequences.
#Weaknesses	May produce slightly less accurate alignments than structure-aware aligners on small, divergent datasets; colinearity assumption can break in genome rearrangements.
#Assumptions	Assumes sequences are homologous and colinear; model-free distance estimation.
#User choices	- --auto: automatically chooses the best algorithm - --maxiterate: number of refinement iterations - --thread: number of CPUs - Input should be unaligned FASTA
#Mafft install and run
conda install mafft
mafft --auto all_mt.fasta > aligned_mt.fasta