

## Workflow

There are already sequenced genomes of _A. mellifera_ and _B. terrestris_ and _L. niger_ in NCBI. We'll have to make something for _L. flavus_ from scratch.

We also need to determine possible isoforms for the genes in the data set, which we will do using tophat for the genomes that have references. 

### _L. niger_ workflow

This specis has no isoform data in NCBI, so we create our own. First we create bowtie2 indexed for the reference genomes using `bowtie2-build`, and then tophat and cufflinks are executed using `tophat.sh` using reference-based mapping to identify alternative splicing.


### _L flavus workflow_
`trinity.sh` assemble lf from raw reads. We are not using a genome guided assembly, since the ln genome is fragmented, as per [recommendation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) of the trinity authors
`transdecoder.sh` predict proteins from the transcripts. 


### Gene expression analysis

We use kallisto for gene expression analysis of transcripts from the NCBI data bases for _A. mellifera_ and _B. terrestis_, for the TopHat assembly of _L. niger_ and for the Trinity assembly of _L. flavus_. 


## Orthology

We handle this by reciprocal blastp hit on the proteins vs the honey bee genome (the best annotated of the bunch) _e.g.,_ `blastp -num_threads 12 -query lf.fa -db amel -outfmt 6 -evalue 1e-4 -max_target_seqs 1`