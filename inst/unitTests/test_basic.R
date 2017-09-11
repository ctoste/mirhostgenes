## test basic stuff...
## library(MirhostDb.Hsapiens.v75.v20)
## DB <- MirhostDb.Hsapiens.v75.v20
library(MirhostDb.Hsapiens.v75.v20)
DB <- MirhostDb.Hsapiens.v75.v20

test_basics <- function(){
    ## this are just basic info methods; we don't expect any special stuff there.
    checkEquals(ensemblVersion(DB), "75")
    checkEquals(mirbaseVersion(DB), "v20")
    listTables(DB)
    listColumns(DB, table="host_tx")
    checkEquals(listColumns(DB, table="adfdkfdf"), NULL)
    checkEquals(length(listDatabases(DB)), 3)
    checkTrue(length(listGenebiotypes(DB)) > 1)
    checkTrue(length(listTxbiotypes(DB)) > 1)
    checkTrue(length(listArrays(DB)) > 1)
    metadata(DB)
    genome(DB)
    seqinfo(DB)
}

