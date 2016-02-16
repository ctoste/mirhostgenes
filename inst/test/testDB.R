## some basic tests for methods of the database.
detachem <- function(x){
    NS <- loadedNamespaces()
    if(any(NS==x)){
        pkgn <- paste0("package:", x)
        detach(pkgn, unload=TRUE, character.only=TRUE)
    }
}
Pkgs <- c("MirhostDb.Hsapiens.v75.v20", "mirhostgenes")
tmp <- sapply(Pkgs, detachem)
tmp <- sapply(Pkgs, library, character.only=TRUE)


## just plain show on the database
MirhostDb.Hsapiens.v75.v20

DB <- MirhostDb.Hsapiens.v75.v20


## what versions?
version(DB, "ensembl")
version(DB, "mirbase")


## list all tables:
listTables(DB)

## columns
listColumns(DB, table="pre_mirna")

listArrays(DB)

organism(DB)

listDatabases(DB)

listGenebiotypes(DB)

listTxbiotypes(DB)

