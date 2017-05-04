# Bioinformatic analysis starting from raw reads.

## Workflow

There are already sequenced genomes of _A. mellifera_ and _B. terrestris_ and _L. niger_ in NCBI. We'll have to make something for _L. flavus_ from scratch.

We also need to determine possible isoforms for the genes in the data set, which we will do using tophat for the genomes that have references. 

### _L. niger_ workflow

This specis has no isoform data in NCBI, so we create our own. First we create bowtie2 indexed for the reference genomes using `bowtie2-build`, and then tophat and cufflinks are executed using `tophat.sh` using reference-based mapping to identify alternative splicing.


### _L flavus_ workflow
`trinity.sh` assemble lf from raw reads. We are not using a genome guided assembly, since the ln genome is fragmented, as per [recommendation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) of the trinity authors
`transdecoder.sh` predict proteins from the transcripts. 


### Gene expression analysis

**kallisto** We use kallisto for gene expression analysis of transcripts from the NCBI data bases for _A. mellifera_ and _B. terrestis_, for the TopHat assembly of _L. niger_ and for the Trinity assembly of _L. flavus_. 

**rsem and ebseq** This is the more traditional approach using the same data sources. We prepare references as ngvector files from the predicted transcripts, as per instructions.

## Orthology

We handle this by reciprocal blastp hit on the proteins vs the honey bee genome (the best annotated of the bunch) _e.g.,_ `blastp -num_threads 12 -query lf.fa -db amel -outfmt 6 -evalue 1e-4 -max_target_seqs 1`

First, we have to generate proteins with the same names as the genes

	gffread GCF_000002195.4_Amel_4.5_genomic.gff -g GCF_000002195.4_Amel_4.5_genomic.fna -F -y - | perl -ne '(/GeneID:(\w+)/ && print ">$1\n") || print' | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//' | grep -v "\." | tr "\t" "\n" > am_prot.fa
	gffread GCF_000214255.1_Bter_1.0_genomic.gff -g GCF_000214255.1_Bter_1.0_genomic.fna -F -y - | perl -ne '(/GeneID:(\w+)/ && print ">$1\n") || print' | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//' | grep -v "\." | tr "\t" "\n" > bt_prot.fa
	gffread GCA_001045655.1_ASM104565v1_genomic.gff -g GCA_001045655.1_ASM104565v1_genomic.fna -F -y - | perl -ne '(/gene=(\w+)/ && print ">$1\n") || print' | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//' | grep -v "\." | tr "\t" "\n" > ln_prot.fa
	cat ../data/assembly/trinity_lf/Trinity.fasta.transdecoder.pep |sed 's/_i.*//' | tr -d '*' > lf_prot.fa


	for i in  *.txt; do cat $i | awk -v OFS=, '{print $1"|"$2,$(NF-1)}' | sort -t , -k1,2 -g | tr "|" "," |  awk -F, '!x[$1$2]++' > `basename $i txt`csv; done

### Generating tables for analysis

#### Gene, transcript and protein IDs

For the same gene we have different IDs for proteins, transcripts and the gene itself. We need to make tables interrelating all these data sources from the gff files.

Pull out honey bee transcripts, genes and bee base id

	cat ref/GCF_000002195.4_Amel_4.5_genomic.gff | awk '$3=="mRNA"' | perl  -ne 'print "$1,$2,$3\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+),BEEBASE:(GB[0-9]+)/' | sort -u > data/tables/amel_mrna.csv
	cat ref/GCF_000002195.4_Amel_4.5_genomic.gff | awk '$3=="CDS"' | perl  -ne 'print "$1,$2,$3\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+),BEEBASE:(GB[0-9]+)/' | sort -u > data/tables/amel_prot.csv
	
We can do the same for bumblebees

	cat ref/GCF_000214255.1_Bter_1.0_genomic.gff | awk '$3=="mRNA"' | perl  -ne 'print "$1,$2\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+)/' > data/tables/bter_mrna.csv
	cat ref/GCF_000214255.1_Bter_1.0_genomic.gff | awk '$3=="CDS"' | perl  -ne 'print "$1,$2\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+)/' > data/tables/bter_prot.csv

For _L. niger_, we obtain all these data from the genbank gff, and the tophat gtf. This is made a bit more difficult by the fact that TopHat used a different gene identifier, which needs to be merged.

	awk '$3=="gene"' ref/GCA_001045655.1_ASM104565v1_genomic.gff | perl  -ne 'print "$1,$2\n" if /ID=(\w+);Name=(\w+)/'  |sed 's/gene/rna/' > data/tables/ln_genes.csv
	grep "NCBI_GP" ref/GCA_001045655.1_ASM104565v1_genomic.gff | perl  -ne 'print "$1,$2\n" if /Parent=(rna\d+);Dbxref=NCBI_GP:(\w+\.\w+)/' |sort -u > data/tables/ln_proteins.csv
	cat data/assembly/tophat/ln/merged.gtf | perl  -ne 'print "$2,$1\n" if /transcript_id "(\w+)".* gene_name "(\w+)"/' | sort -u > data/tables/ln_transcripts.csv
	join -1 1 -2 2 -t , <(sort -k1,1 -t , data/tables/ln_transcripts.csv ) <(sort -k2,2 -t , data/tables/ln_genes.csv )  |sort -u > data/tables/ln_gene_transcript_rna.csv
	join -1 1 -2 2 -t , <(sort -k1,1 -t , data/tables/ln_transcripts.csv ) <(sort -k2,2 -t , data/tables/ln_genes.csv )  |sort -u > data/tables/ln_gene_transcript_rna_protein.csv

For _L. flavus_ we extract all the relevant details from the TransDecoder output

	grep ">" data/assembly/trinity_lf/Trinity.fasta.transdecoder.pep |tr -d '>' | awk '{split($1, a, "|"); gene=a[1]; sub(/_i[0-9]+/,"",gene); print gene","a[1]","$1}' > data/tables/lf_mrna_prot.csv


join -j 1 -t , <(sort -k1,1 -t , data/tables/ln_transcripts.csv ) <(sort -k1,1 -t , data/tables/ln_proteins.csv) 