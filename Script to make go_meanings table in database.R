# Simple script to make the table 'go_meanings', which has ontology and terms for all the GOs in the bee_go table
# Luke wrote this because loading GO.db causes many annoying conflicts with other packages, e.g. dplyr, that we use a lot
# Restart R after running this, so that GO.db is definitely gone!

library(GO.db)
library(dplyr)
my_db <- src_sqlite("data/queen pheromone.db")

go.meanings <- suppressMessages(AnnotationDbi::select(GO.db, unique((tbl(my_db, "bee_go") %>% dplyr::select(GO) %>% collect(n=Inf) %>% as.data.frame())[,1]), c("GOID", "ONTOLOGY", "TERM")))
names(go.meanings) <- c("GO", "ontology", "term")

# Add/overwrite the new table to the database
my_db$con %>% db_drop_table(table = "go_meanings")
copy_to(my_db, go.meanings, "go_meanings", temporary = FALSE)



# For the same reason (conflicts), let's use this script to make the 'GeneSetCollection' object 
# that we need to do GO term enrichment on non-model organisms. This follows the instructions from:
# http://bioconductor.org/packages/2.11/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf

library(GSEABase)
library(AnnotationDbi)

goframeData <- tbl(my_db, "bee_go") %>% mutate(evidence = "ND") %>% dplyr::select(GO, evidence, gene) %>% as.data.frame()
goFrame <- GOFrame(goframeData, organism="Homo sapiens")
goAllFrame <- GOAllFrame(goFrame)
gene_set_collection <- GeneSetCollection(goAllFrame, setType = GOCollection())
save(gene_set_collection, file="data/gene_set_collection.RData")

