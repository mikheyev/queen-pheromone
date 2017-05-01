

## Workflow

There are already sequenced genomes of _A. mellifera_ and _B. terrestris_ and _L. niger_ in NCBI. We'll have to make something for _L. flavus_ from scratch.

We also need to determine possible isoforms for the genes in the data set, which we will do using tophat for the genomes that have references. 

### _A. mellifera_ and _B. terrestris_ and _L. niger_ workflow

first we create bowtie2 indexed for the reference genomes using `bowtie2-build`


### _L flavus workflow_
`trinity.sh` assemble lf from raw reads. We are not using a genome guided assembly, since the ln genome is fragmented, as per [recommendation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) of the trinity authors
`transdecoder.sh` predict proteins from 
