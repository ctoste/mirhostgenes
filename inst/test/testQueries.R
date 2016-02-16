## check namespace.
detachem <- function(x){
    NS <- loadedNamespaces()
    if(any(NS==x)){
        pkgn <- paste0("package:", x)
        detach(pkgn, unload=TRUE, character.only=TRUE)
    }
}
Pkgs <- c("MirhostDb.Hsapiens.v75.v20", "mirhostgenes", "ensembldb")
tmp <- sapply(Pkgs, detachem)
tmp <- sapply(Pkgs, library, character.only=TRUE)

DB <- MirhostDb.Hsapiens.v75.v20

######********************************************
## testing the queries. that's really an internal function that should never be called
## directly.
#####
## 1)
tojoin <- c("pre_mirna", "host_gene")
## adds host_tx
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q
library(RSQLite)
R1 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## OK.
## the same, but starting from the host_gene table!
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join",
                                      start.table="host_gene")
Q
R2 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## compare the results
nrow(R1)
nrow(R2)
## have more for R1 where we started with pre_mirna
sum(is.na(R1$gene_id))
sum(is.na(R2$gene_id))

#####
## 2)
tojoin <- c("mat_mirna", "host_gene")
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q
R1 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## the same but starting from host_gene
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join",
                                      start.table="host_gene")
Q
R2 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## what if we started from mat_mirna?
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join",
                                      start.table="mat_mirna")
Q
## compare the results
nrow(R1)
nrow(R2)
## have more for R1 where we started with pre_mirna
sum(is.na(R1$gene_id))
sum(is.na(R2$gene_id))

#####
## 3)
tojoin <- c("mirfam", "mat_mirna")
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q
R1 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## start from mirfam
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join",
                                      start.table="mirfam")
Q
R2 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## compare the results
nrow(R1)
nrow(R2)
## have more for R1 where we started with pre_mirna
sum(is.na(R1$mirfam_id))
sum(is.na(R2$mirfam_id))

#####
## 4)
## adds pre_mirna, host_tx
tojoin <- c("mat_mirna", "mirfam", "host_gene")
## adds pre_mirna, host_tx
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q
R1 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## same but start from mirfam...
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join",
                                      start.table="mirfam")
Q
R2 <- dbGetQuery(dbconn(DB), paste("select * from", Q))
## we're supposed to get less in R2
## compare the results
nrow(R1)
nrow(R2)
## have more for R1 where we started with pre_mirna
sum(is.na(R1$mirfam_id))
sum(is.na(R2$mirfam_id))

#####
## 5)
##
tojoin <- c("mat_mirna", "array_feature")
## adds pre_mirna, host_tx
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q
tojoin <- c( "host_gene", "array_feature" )
Q <- mirhostgenes:::joinQueryOnTables(DB, tojoin, join="left join")
Q

##
######********************************************




## testing the queries to see whether we get what we want...
## Note: these functions are internal.
Query <- mirhostgenes:::.buildQuery(DB, columns=c("gene_name", "sequence",
                                             "pre_mirna_name", "pre_mirna_id",
                                             "mat_mirna_name"),
                                    filter=list(SeqendFilter(123, condition=">",
                                        feature="mat_mirna"),
                                        GenenameFilter("SMC4")),
                                    order.by="mat_mirna_seq_end")
Query

Res <- dbGetQuery(dbconn(DB), Query)
Res

## now we want to get all mature_mirnas for host gene SMC4, along with its host transcripts, ordered by mat_mirna_seq_start.
Query <- mirhostgenes:::.buildQuery(DB, columns=c("tx_id", "gene_name",
                                             "pre_mirna_id", "pre_mirna_name",
                                             "mat_mirna_name", "mat_mirna_seq_start"),
                                    filter=list(GenenameFilter("SMC4")),
                                    order.by="mat_mirna_seq_start")
Query
Res <- dbGetQuery(dbconn(DB), Query)
Res

## get all pre_mirnas
Query <- mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "pre_mirna_name"))
Allpres <- dbGetQuery(dbconn(DB), Query)
nrow(Allpres)
length(unique(Allpres$pre_mirna_id))

## now build a query joining pre_mirnas with host_tx
Query <- mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "pre_mirna_name", "tx_id"))
Res.join <- dbGetQuery(dbconn(DB), Query)
nrow(Res.join)
length(unique(Res.join$pre_mirna_id))

## what if we used a full join??? doesn't work, only left joins!!!
Query <- mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "pre_mirna_name", "tx_id"),
                                    join="left join")
Res.left <- dbGetQuery(dbconn(DB), Query)
nrow(Res.left)
length(unique(Res.left$pre_mirna_id))
sum(is.na(Res.left$tx_id))
## OK, so we get also some without tx_id

## if we say now that we want to start with host_tx we will end up with the same as above.
Query <- mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "pre_mirna_name", "tx_id"),
                                    join="left join", start.table="host_tx")
Res.join <- dbGetQuery(dbconn(DB), Query)
nrow(Res.join)
length(unique(Res.join$pre_mirna_id))

Query <- mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "pre_mirna_name",
                                             "tx_id", "mirfam_name"))
Query
Res <- dbGetQuery(dbconn(DB), Query)
nrow(Res)
length(unique(Res$pre_mirna_id))
length(unique(Res$mirfam_name))

## what if we used a left join here?
Query <- gsub(Query, pattern="join", replacement="left join")
Res.left <- dbGetQuery(dbconn(DB), Query)
nrow(Res.left)
length(unique(Res.left$pre_mirna_id))
length(unique(Res.left$mirfam_name))

dbGetQuery(dbconn(DB), "select count(distinct mirfam_name) from mirfam")
dbGetQuery(dbconn(DB), "select count(distinct pre_mirna_id) from mirfam")
dbGetQuery(dbconn(DB), "select count(distinct pre_mirna_id) from pre_mirna")
## so, with all left joins we get as much as possible... it's just important from which table we start...





##*************************************************
##
## Internal functions and methods.
##
##*************************************************
cat("Testing internal stuff...")
mirhostgenes:::tablesByDegree(DB)

mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "array_id", "feature_id"))

mirhostgenes:::.buildQuery(DB, columns=c("pre_mirna_id", "database", "gene_id"))

mirhostgenes:::.buildQuery(DB, columns=c("gene_id", "tx_biotype"))

mirhostgenes:::.buildQuery(DB, columns=c("gene_id", "tx_id"))

mirhostgenes:::.buildQuery(DB, columns=c("gene_id", "tx_id", "gene_name"))

mirhostgenes:::.buildQuery(DB, columns=listColumns(DB, "host_gene"))

## what if we use a GeneidFilter and columns gene_id, tx_id, it should
## filter on tx.gene_id
mirhostgenes:::.buildQuery(DB, columns=c("gene_id", "tx_id"),
                           filter=list(GeneidFilter("a")))
## if we have more attrs from the gene table it should be host_gene.gene_id
mirhostgenes:::.buildQuery(DB, columns=c("gene_id", "tx_id", "gene_name", "gene_biotype"),
                           filter=list(GeneidFilter("a")))


cat("done.\n")





