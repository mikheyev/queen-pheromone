# R-based analysis of queen pheromone data

`data/kallisto`  gene expression results
`data/orthology` results from reciprocal best blast


### Generating tables for analysis

#### Gene, transcript and protein IDs

For the same gene we have different IDs for proteins, transcripts and the gene itself. We need to make tables interrelating all these data sources from the gff files

Pull out honey bee transcripts, genes and bee base id

	cat ref/GCF_000002195.4_Amel_4.5_genomic.gff | awk '$3=="mRNA"' | perl  -ne 'print "$1,$2,$3\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+),BEEBASE:(GB[0-9]+)/' > data/tables/amel_mrna.csv
	cat ref/GCF_000002195.4_Amel_4.5_genomic.gff | awk '$3=="CDS"' | perl  -ne 'print "$1,$2,$3\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+),BEEBASE:(GB[0-9]+)/' > data/tables/amel_prot.csv

We can do the same for bumblebees

	cat ref/GCF_000214255.1_Bter_1.0_genomic.gff | awk '$3=="mRNA"' | perl  -ne 'print "$1,$2\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+)/' > data/tables/bter_mrna.csv
	cat ref/GCF_000214255.1_Bter_1.0_genomic.gff | awk '$3=="CDS"' | perl  -ne 'print "$1,$2\n" if /GeneID:([0-9]+),Genbank:(.._[0-9]+\.[0-9]+)/' > data/tables/bter_prot.csv

For _L. niger_

For _L. flavus_ we extract all the relevant details from the TransDecoder output

	grep ">" ~/pj/qp/data/assembly/trinity_lf/Trinity.fasta.transdecoder.pep |tr -d '>' | awk '{split($1, a, "|"); gene=a[1]; sub(/_i[0-9]+/,"",gene); print gene","a[1]","$1}' > data/tables/lf_mrna_prot.csv

