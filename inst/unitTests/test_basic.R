## test basic stuff...
library(MirhostDb.Hsapiens.v75.v20)
DB <- MirhostDb.Hsapiens.v75.v20

test_basics <- function(){
    ## this are just basic info methods; we don't expect any special stuff there.
    checkEquals(ensemblVersion(DB), "75")
    checkEquals(mirbaseVersion(DB), "v20")
    listTables(DB)
    listColumns(DB, table="host_tx")
    checkEquals(listColumns(DB, table="adfdkfdf"), NULL)
    listDatabases(DB)
    listGenebiotypes(DB)
    listTxbiotypes(DB)
    listArrays(DB)
    metadata(DB)
    genome(DB)
    seqinfo(DB)
    return(TRUE)
}

