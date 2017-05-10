# Simple script to make the table 'go_meanings', which has ontology and terms for all the GOs in the bee_go table
# Luke wrote this because loading GO.db causes many annoying conflicts with other packages
# Restart R after running this, so that GO.db is definitely gone!

library(GO.db)
library(dplyr)
go.meanings <- suppressMessages(AnnotationDbi::select(GO.db, unique((tbl(my_db, "bee_go") %>% dplyr::select(GO) %>% collect(n=Inf) %>% as.data.frame())[,1]), c("GOID", "ONTOLOGY", "TERM")))
names(go.meanings) <- c("GO", "ontology", "term")

# Add/overwrite the new table to the database
my_db$con %>% db_drop_table(table = "go_meanings")
copy_to(my_db, go.meanings, "go_meanings", temporary = FALSE)
