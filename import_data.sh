#!/bin/bash

# import data for downstream statistica analysis

## treatments
echo "DROP TABLE IF EXISTS treatments; CREATE TABLE treatments (id varchar(4) not null, species varchar(2) not null, colony integer not null, treatment varchar, PRIMARY KEY (id));" > temp.sql
echo ".separator ," >> temp.sql
echo '.import "data/treatments.csv" treatments' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 

## bee go terms derived from uniprot
echo "DROP TABLE IF EXISTS bee_go; CREATE TABLE bee_go (gene varchar not null, GO varchar not null);" > temp.sql
echo ".separator ," >> temp.sql
echo '.import "data/am go.csv" bee_go' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 

## bee gene names
echo "DROP TABLE IF EXISTS bee_names; CREATE TABLE bee_names (gene varchar not null, name varchar not null);" > temp.sql
echo ".separator ," >> temp.sql
echo '.import "data/am names.csv" bee_names' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 

## reciprocal best blast results

for i in data/orthologs/*csv; do
	name1=`basename $i | cut -c-2`
	name2=`basename $i | cut -c4,5`
	base=`basename $i .csv`  
	echo "DROP TABLE IF EXISTS $base; CREATE TABLE $base ($name1 varchar(30) not null, $name2 varchar(30) not null, evalue real not null, PRIMARY KEY ($name1, $name2));" > temp.sql
	echo ".separator ," >> temp.sql
	echo ".import \"$i\" $base" >> temp.sql
	cat temp.sql | sqlite3 data/queen\ pheromone.db 
done

## isoforms 
for i in data/isoforms/*csv; do
	base="isoforms_"`basename $i .csv`
	echo "DROP TABLE IF EXISTS $base; CREATE TABLE $base (gene varchar not null, isoform varchar not null, PRIMARY KEY (isoform));" > temp.sql
	echo ".separator ," >> temp.sql
	echo ".import \"$i\" $base" >> temp.sql
	cat temp.sql | sqlite3 data/queen\ pheromone.db
done

# ## rsem results (TPM)
echo "DROP TABLE IF EXISTS rsem_am; CREATE TABLE rsem_am (gene varchar, am2 real, am4 real, am6 real, am1 real, am3 real, am5 real, PRIMARY KEY (gene));" > temp.sql
echo '.separator "\t"' >> temp.sql
echo '.import "data/rsem/am.genes.matrix" rsem_am' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db

echo "DROP TABLE IF EXISTS rsem_bt; CREATE TABLE rsem_bt (gene varchar, bt2 real, bt4 real, bt5 real, bt8 real, bt10 real, bt1 real, bt3 real, bt6 real, bt7 real, bt9 real, PRIMARY KEY (gene));" > temp.sql
echo '.separator "\t"' >> temp.sql
echo '.import "data/rsem/bt.genes.matrix" rsem_bt' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 	

echo "DROP TABLE IF EXISTS rsem_ln; CREATE TABLE rsem_ln (gene varchar, ln1 real, ln3 real, ln5 real, ln7 real, ln9 real, ln11 real, ln2 real, ln4 real, ln6 real, ln8 real, ln10 real, ln12 real, PRIMARY KEY (gene));" > temp.sql
echo '.separator "\t"' >> temp.sql
echo '.import "data/rsem/ln.genes.matrix" rsem_ln' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 	

echo "DROP TABLE IF EXISTS rsem_lf; CREATE TABLE rsem_lf (gene varchar, lf10 real, lf11 real, lf12 real, lf13 real, lf14 real, lf15 real, lf16 real, lf1 real, lf2 real, lf3 real, lf4 real, lf5 real, lf6 real, lf7 real, lf8 real, PRIMARY KEY (gene));" > temp.sql
echo '.separator "\t"' >> temp.sql
echo '.import "data/rsem/lf.genes.matrix" rsem_lf' >> temp.sql
cat temp.sql | sqlite3 data/queen\ pheromone.db 	

## ebseq results by by gene, but without correction

for i in  data/rsem/??.genes.ebseq; do
	base="ebseq_gene_"`basename $i .genes.ebseq`
	echo "DROP TABLE IF EXISTS $base; CREATE TABLE $base (gene varchar, PPEE real, PPDE real, PostFC real, RealFC real, PRIMARY KEY (gene));" > temp.sql
	echo '.separator "\t"' >> temp.sql
	echo ".import \"$i\" $base" >> temp.sql
	cat temp.sql | sqlite3 data/queen\ pheromone.db 	
done

## ebseq results by isoform
for i in  data/rsem/??.padj.ebseq; do
	base="ebseq_padj_isoform_"`basename $i .padj.ebseq`
	echo "DROP TABLE IF EXISTS $base; CREATE TABLE $base (isoform varchar, PPEE real, PPDE real, PostFC real, RealFC real, PRIMARY KEY (isoform));" > temp.sql
	echo '.separator "\t"' >> temp.sql
	echo ".import \"$i\" $base" >> temp.sql
	cat temp.sql | sqlite3 data/queen\ pheromone.db 	
done

## ebseq results by gene
for i in  data/rsem/??.genes.padj.ebseq; do
	base="ebseq_padj_gene_"`basename $i .genes.padj.ebseq`
	echo "DROP TABLE IF EXISTS $base; CREATE TABLE $base (gene varchar, PPEE real, PPDE real, PostFC real, RealFC real, PRIMARY KEY (gene));" > temp.sql
	echo '.separator "\t"' >> temp.sql
	echo ".import \"$i\" $base" >> temp.sql
	cat temp.sql | sqlite3 data/queen\ pheromone.db 	
done
rm temp.sql
