# Simple script to make the table 'go_meanings', which has ontology and terms for all the GOs in the bee_go table
# Luke wrote this because loading GO.db causes many annoying conflicts with other packages, e.g. dplyr, that we use a lot.
# Restart R after running this, so that GO.db is definitely gone!

library(GO.db)
library(dplyr)
my_db <- src_sqlite("data/queen pheromone.db")

go.meanings <- suppressMessages(AnnotationDbi::select(GO.db, 
                                                      unique((tbl(my_db, "bee_go") %>% 
                                                                dplyr::select(GO) %>% 
                                                                collect(n=Inf) %>% 
                                                                as.data.frame())[,1]), c("GOID", "ONTOLOGY", "TERM")))
names(go.meanings) <- c("GO", "ontology", "term")

# Add/overwrite the new table to the database
my_db$con %>% db_drop_table(table = "go_meanings")
copy_to(my_db, go.meanings, "go_meanings", temporary = FALSE)



# For the same reason (conflicts), let's use this script to make the 'GeneSetCollection' object 
# that we need to do GO term enrichment on non-model organisms. This follows the instructions from:
# http://bioconductor.org/packages/2.11/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf

library(GSEABase)
library(AnnotationDbi)

goframeData <- tbl(my_db, "bee_go") %>% 
  mutate(evidence = "ND") %>% 
  dplyr::select(GO, evidence, gene) %>%
  as.data.frame()
goFrame <- GOFrame(goframeData, organism="Homo sapiens")
goAllFrame <- GOAllFrame(goFrame)
gene_set_collection <- GeneSetCollection(goAllFrame, setType = GOCollection())
save(gene_set_collection, file="data/gene_set_collection.RData")



####### Now to get KEGG pathway annotations for each Apis gene

# First, we need to get the Entrez names for each gene (many are given by their Beebase names)
# I found the beebase-entrez mappings using help from this thread, giving me the am.gene_info.txt file: http://seqanswers.com/forums/archive/index.php/t-35472.html

# get entrez to beebase mappings
entrez.tbl <- read.delim("./data/am.gene_info.txt", stringsAsFactors = FALSE)[,c(2,5,6)]
names(entrez.tbl) <- c("entrez.id", "beebase1", "beebase2")
entrez.tbl$beebase1 <- str_extract(entrez.tbl$beebase1, "GB[:digit:]*")
entrez.tbl$beebase2 <- gsub("BEEBASE:", "", entrez.tbl$beebase2) 
entrez.tbl$beebase2[entrez.tbl$beebase2 == "-"] <- NA
entrez.tbl <- entrez.tbl %>% 
  filter(!(is.na(beebase1) & is.na(beebase2))) %>% 
  filter(!is.na(entrez.id))
entrez.tbl <- tbl(my_db, "rsem_am") %>% 
  select(gene) %>% collect(n=Inf) %>% 
  left_join(entrez.tbl, by = c("gene" = "beebase2")) %>% 
  select(gene, entrez.id) 
entrez.tbl$entrez.id[is.na(entrez.tbl$entrez.id)] <- entrez.tbl$gene[is.na(entrez.tbl$entrez.id)]

# Get the kegg to entrez mapping
kegg.tbl <- melt(kegg.gsets(species = "ame", id.type = "kegg")$kg.sets)
kegg.tbl <- data.frame(entrez.id = kegg.tbl$value,
                       kegg = substr(kegg.tbl$L1, 1, 8),
                       meaning = substr(kegg.tbl$L1, 10, 
                                        nchar(kegg.tbl$L1)), stringsAsFactors = FALSE)

# Make a table of all the kegg terms we know for the genes in our set
kegg.tbl <- entrez.tbl %>% 
  left_join(kegg.tbl) %>% 
  filter(!is.na(kegg)) %>% 
  select(-entrez.id) %>% 
  mutate(kegg = gsub("ame", "", kegg))
copy_to(my_db, kegg.tbl, "bee_kegg", temporary = FALSE)

# Make the GeneSetCollection object needed to do KEGG enrichment with the GOstats package
library(KEGG.db)
keggframeData <- tbl(my_db, "bee_kegg") %>% 
  dplyr::select(kegg, gene) %>% 
  collect(n=Inf) %>% as.data.frame()
keggFrame <- KEGGFrame(keggframeData, organism="Homo sapiens")
gene_set_collection_kegg <- GeneSetCollection(keggFrame, setType = KEGGCollection())
save(gene_set_collection_kegg, file="data/gene_set_collection_kegg.RData")


