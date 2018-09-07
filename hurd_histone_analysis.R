library(ape)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(EBSeq)
library(reshape2)

########### Analysis of the ChIPseq data (refernce number is GSE110640, https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE110640)

# Open the ChIPseq data from Hurd et al. 2018 Genome Biology
files <- list.files("data/hurd_histone_data/GSE110640_RAW/")

# Make a list of all the ChIPseq samples
chipseq_samples <- data.frame(sample = 1:length(files), strsplit(files, split = "_") %>% 
  do.call("rbind", .) %>% as.data.frame() %>% 
    mutate(V3 = gsub("[.]txt[.]gz", "", V3),
           V2 = substr(V2,1,4))) %>%
  rename(sample=V1, caste=V2, histone=V3)

# Open all the ChIP files into a list
x <- lapply(list.files("data/hurd_histone_data/GSE110640_RAW/", full.names = TRUE), 
            function(x) read_tsv(x, col_names = c("seq_id", "start", "end", "value")))
# Make sure there is one observation per sample per site, by averaging across any duplicates
x <- lapply(x, function(x) group_by(x, seq_id, start, end) %>% summarise(value = mean(value)))

# Concatenate all the ChIP-seq data into one data frame using full_join
histone_data <- full_join(x[[1]], x[[2]], by = c("seq_id", "start", "end")) %>%
  dplyr::rename(sample1 = value.x, sample2 = value.y) 
for(i in 3:length(files)){
  print(i)
  histone_data <- full_join(histone_data, x[[i]], by = c("seq_id", "start", "end")) 
  names(histone_data)[names(histone_data) == "value"] <- paste("sample", i, sep="")
}
histone_data$annotation <- vector(mode = "list", length = nrow(histone_data))


# Open the A. mellifera 4.5 genome annotation (same one as used by Hurd et al)
apis_gff <- read.gff("data/hurd_histone_data/GCF_000002195.4_Amel_4.5_genomic.gff.gz") %>% 
  filter(type == "gene") %>% select(seqid, start, end, attributes) %>%
  mutate(attributes = lapply(
    lapply(str_extract_all(attributes, "GB[:digit:]+"), unique), 
    function(x) ifelse(length(x) > 0, x, NA))) %>%
  rename(seq_id = seqid, gene = attributes) %>%
  filter(!is.na(gene), end - start < 50000) # only keep measurements where we know a gene nearby

# Remove any contigs that are in the ChIP-seq data but not the annotation
histone_data <- histone_data[histone_data$seq_id %in% apis_gff$seq_id, ] %>% ungroup()

# Annotate the ChIPseq data with gene IDs
histone_data <- parallel::mclapply(1:nrow(histone_data), function(i, apis_gff, histone_data){
  foc_histone <- histone_data[i,]
  foc_gff <- apis_gff %>% filter(seq_id == foc_histone$seq_id)
  hits <- which(
    (foc_histone$start > foc_gff$start & foc_histone$start < foc_gff$end) |
      (foc_histone$end > foc_gff$start & foc_histone$end < foc_gff$end) 
  )
  num_hits <- length(hits)
  if(num_hits == 1) return(data.frame(gene = as.character(apis_gff$gene[hits]), histone_data[i, ] %>% select(starts_with("sample"))))
  if(num_hits > 1)  return(data.frame(gene = as.character(apis_gff$gene[hits]), histone_data[rep(i, num_hits), ] %>% select(starts_with("sample"))))
  else return(NULL)
}, mc.cores = 8, apis_gff = apis_gff, histone_data = histone_data) %>% 
  bind_rows() %>%
  filter(!is.na(gene)) # get rid of epimarks with no known gene nearby

histone_data_by_gene <- histone_data %>% ungroup() %>%
   group_by(gene) %>%
  summarise_all(sum, na.rm = T) %>% arrange(gene)
genes <- histone_data_by_gene$gene
histone_data_by_gene <- as.matrix(histone_data_by_gene[,-1])
rownames(histone_data_by_gene) <- genes

run_ebseq_histone <- function(histone, chipseq_samples, histone_data_by_gene){
  histone_data_by_gene <- histone_data_by_gene[, chipseq_samples$sample[chipseq_samples$histone == histone]]
  
  mean_queens <- data.frame(histone_data_by_gene[,1:2]) %>% tibble::rownames_to_column()
  names(mean_queens) <- c("gene", "Q1", "Q2")
  mean_workers <- data.frame(histone_data_by_gene[,3:4]) %>% tibble::rownames_to_column()
  names(mean_workers) <- c("gene", "W1", "W2")
  
  mean_queens <- mean_queens %>%
    group_by(gene) %>%
    summarise(mean_queens = round(Q1 + Q2, 2))
  
  mean_workers <- mean_workers %>%
    group_by(gene) %>%
    summarise(mean_workers = round(W1 + W2, 2))
  
  Sizes <- MedianNorm(histone_data_by_gene)
  EBOut <- EBTest(Data = histone_data_by_gene, 
                  Conditions = as.factor(rep(c("96hQ","96hW"), each = 2)),
                  sizeFactors = Sizes, maxround = 50)
  
  histone_ebseq <- GetDEResults(EBOut, FDR = 0.05)
  GeneFC <- PostFC(EBOut)
  left_join(melt(GeneFC$PostFC) %>% tibble::rownames_to_column() %>% rename(PostFC = value),
                               melt(GeneFC$RealFC) %>% tibble::rownames_to_column() %>% rename(RealFC = value),  by = "rowname") %>%
    
    left_join(as.data.frame(histone_ebseq[[2]]) %>% tibble::rownames_to_column(),  by = "rowname") %>%
    rename(gene = rowname) %>% arrange(-abs(log2(PostFC))) %>%
    left_join(mean_queens, by = "gene") %>% left_join(mean_workers, by = "gene")
}

# EB-seq on the CHipSeq results
H3K4me3_ebseq <- run_ebseq_histone("H3K4me3", chipseq_samples, histone_data_by_gene)
H3K27ac_ebseq <- run_ebseq_histone("H3K27ac", chipseq_samples, histone_data_by_gene)
H3K36me3_ebseq <- run_ebseq_histone("H3K36me3", chipseq_samples, histone_data_by_gene)




########### Analysis of the RNA-seq data (refernce number is GSE110641, https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE110641)

my_db <- src_sqlite("data/queen_pheromone.db")

files <- list.files("data/hurd_histone_data/GSE110641_RAW/")
sampleID <- data.frame(sample = 1:length(files), strsplit(files, split = "_") %>% 
                         do.call("rbind", .) %>% as.data.frame()) %>% 
  mutate(V2 = substr(V2,1,7)) %>% rename(id = V1, caste = V2)

expression_data_isoforms <- lapply(list.files("data/hurd_histone_data/GSE110641_RAW/", full.names = TRUE), 
                          function(x) read_tsv(x) %>% select(target_id, tpm)) %>% 
  purrr::reduce(left_join, by = "target_id") 

expression_data_isoforms <- expression_data_isoforms[rowSums(expression_data_isoforms[,-1]) > 10^-5, ]

expression_data_genes <- left_join(expression_data_isoforms,
                                   tbl(my_db, "isoforms_am") %>% 
                                     collect(), 
                                   by = c("target_id" = "isoform")) %>%
  select(-target_id) %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise_all(sum)

gene <- expression_data_genes$gene
expression_data_genes <- as.matrix(expression_data_genes[,-1])
rownames(expression_data_genes) <- gene
colnames(expression_data_genes) <- sampleID$id

isoforms <- expression_data_isoforms$target_id
expression_data_isoforms <- as.matrix(expression_data_isoforms[,-1])
rownames(expression_data_isoforms) <- isoforms
colnames(expression_data_isoforms) <- sampleID$id


# Gene-level
Sizes <- MedianNorm(expression_data_genes)
EBOut <- EBTest(Data=expression_data_genes, 
                Conditions = as.factor(rep(c("QLH","WLH"), times = c(5,4))),
                sizeFactors=Sizes, maxround=50)
hurd_ebseq_gene <- GetDEResults(EBOut, FDR = 0.05)
GeneFC <- PostFC(EBOut)
hurd_ebseq_gene <- left_join(melt(GeneFC$PostFC) %>% tibble::rownames_to_column() %>% rename(PostFC = value),
                             melt(GeneFC$RealFC) %>% tibble::rownames_to_column() %>% rename(RealFC = value),  by = "rowname") %>%
  left_join(as.data.frame(hurd_ebseq_gene[[2]]) %>% tibble::rownames_to_column(),  by = "rowname") %>%
  rename(gene = rowname) %>% 
  filter(!is.na(PPEE)) %>%
  arrange(-PPDE)

# Isoform level
isoforms_am <- tbl(my_db, "isoforms_am") %>% collect()
NgList <- GetNg(isoforms_am$isoform, isoforms_am$gene)
IsoNgTrun <- NgList$IsoformNgTrun
IsoSizes <- MedianNorm(expression_data_isoforms)

IsoEBOut <- EBTest(Data = expression_data_isoforms, 
                   NgVector = IsoNgTrun,
                   Conditions = as.factor(rep(c("QLH","WLH"), times = c(5,4))),
                   sizeFactors = IsoSizes, maxround = 50)
hurd_ebseq_isoform <- GetDEResults(IsoEBOut, FDR = 0.05)
IsoformFC <- PostFC(IsoEBOut)
hurd_ebseq_isoform  <- left_join(melt(IsoformFC$PostFC) %>% tibble::rownames_to_column() %>% rename(PostFC = value),
                                 melt(IsoformFC$RealFC) %>% tibble::rownames_to_column() %>% rename(RealFC = value),  by = "rowname") %>%
  left_join(as.data.frame(hurd_ebseq_isoform[[2]]) %>% tibble::rownames_to_column(),  by = "rowname") %>%
  rename(gene = rowname) %>% arrange(-PPDE)


H3K4me3_ebseq
H3K27ac_ebseq 
H3K36me3_ebseq

hurd_ebseq_gene
hurd_ebseq_isoform

my_db$con %>% db_drop_table(table = "H3K4me3_ebseq")
copy_to(my_db, H3K4me3_ebseq, "H3K4me3_ebseq", temporary = FALSE)

my_db$con %>% db_drop_table(table = "H3K27ac_ebseq")
copy_to(my_db, H3K27ac_ebseq, "H3K27ac_ebseq", temporary = FALSE)

my_db$con %>% db_drop_table(table = "H3K36me3_ebseq")
copy_to(my_db, H3K36me3_ebseq, "H3K36me3_ebseq", temporary = FALSE)

my_db$con %>% db_drop_table(table = "hurd_ebseq_gene")
copy_to(my_db, hurd_ebseq_gene, "hurd_ebseq_gene", temporary = FALSE)


tbl(my_db, "bee_names") %>% filter(gene %in% hurd_ebseq_gene$gene[hurd_ebseq_gene$PPDE>=0.95]) %>% collect %>% as.data.frame() %>% .$name
