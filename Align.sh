#The alignment tool is MAFFT since the data is huge. MAFFT is fast but usually less accurate 
cat /work/cyu/poolseq/PPalign_output/direct_consensus/*.fasta > all_mt.fasta
mafft --auto all_mt.fasta > aligned_mt.fasta