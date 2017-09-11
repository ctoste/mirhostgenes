## testing basic functionality of the package.
## these function just call the methods, but do not test the
## results.

library("MirhostDb.Hsapiens.v75.v20")
DB <- MirhostDb.Hsapiens.v75.v20


test_matmirnas <- function(){
    ## get all mature miRNAs
    checkTrue(nrow(matmirnas(DB)) > 1)
    ## get mat_mirna and pre_mirna entries for mature miRNA MIMAT0000062
    MM <- matmirnas(DB, columns=unique(c(listColumns(DB, "mat_mirna"),
                             listColumns(DB, "pre_mirna"))),
                    filter=MatMirnaIdFilter("MIMAT0000062"))
    checkEquals(unique(MM$mat_mirna_id), "MIMAT0000062")
    ## Test whether the mature miRNA sequences are correctly calculated.
    MM <- matmirnas(DB, columns=c("mat_mirna_id", "pre_mirna_name",
                                  "sequence"))
    MM <- split(MM, MM$mat_mirna_id)
    Counts <- unlist(lapply(MM, function(z){
        length(unique(z$sequence))
    }))
    checkEquals(length(unique(Counts)), 1)

    ## Check columns, and make sure we're returning just those columns that
    ## are requested.
    ResAll <- matmirnas(DB)
    Res <- matmirnas(DB, columns="mat_mirna_name")
    checkEquals(colnames(Res), "mat_mirna_name")
    ## Nrow of MMall has to be larger, as the same mature miRNA can be encoded
    ## in several genomic loci
    checkTrue(nrow(ResAll) > nrow(Res))
    ## But still all of em have to be present.
    checkEquals(sort(unique(Res$mat_mirna_name)), sort(unique(ResAll$mat_mirna_name)))
    ## And now for columns that are not in mat_mirna
    Res <- matmirnas(DB, columns="pre_mirna_name")
    checkEquals(colnames(Res), "pre_mirna_name")

    ## matmirnasBy
    ResAll <- matmirnasBy(DB)
    Res <- matmirnasBy(DB, columns="mat_mirna_name")
    ## checkEquals(colnames(Res[[1]]), "mat_mirna_name")
}

test_premirnas <- function(){
    PM <- premirnas(DB)
    PM <- premirnasBy(DB)
    PM <- premirnasBy(DB, by="mirfam", use.names=TRUE)
    PM <- premirnasBy(DB, by="host_tx",
                      filter=list(DatabaseFilter("core")))

    ## Check columns, and make sure we're returning just those columns that
    ## are requested.
    ResAll <- premirnas(DB)
    Res <- premirnas(DB, columns="pre_mirna_name")
    checkEquals(colnames(Res), "pre_mirna_name")
    ## Nrow of MMall has to be larger, as the same mature miRNA can be encoded
    ## in several genomic loci
    checkTrue(nrow(ResAll) > nrow(Res))
    ## But still all of em have to be present.
    checkEquals(sort(unique(Res$pre_mirna_name)), sort(unique(ResAll$pre_mirna_name)))
    ## And now for columns that are not in mat_mirna
    Res <- premirnas(DB, columns="mat_mirna_name")
    checkEquals(colnames(Res), "mat_mirna_name")
    ## premirnasBy
    ResAll <- premirnasBy(DB)
    Res <- premirnasBy(DB, columns="pre_mirna_name")
}

test_hostgenes <- function(){
    HG <- hostgenes(DB)
    ## this fetches only from the host_tx table
    HG <- hostgenes(DB,columns=c("gene_id", "tx_id"))
    hostgenes(DB,columns=c("gene_id", "tx_id"),
              filter=list(DatabaseFilter("core")))
    HG <- hostgenesBy(DB, by="pre_mirna")
    ## that query is amazingly slow... but why???
    ## HG <- hostgenesBy(DB, by="pre_mirna",
    ##                   filter=list(DatabaseFilter("core")))
    ## that too...
    ##HG <- hostgenesBy(DB, by="host_tx",
    ##                  filter=list(DatabaseFilter("core")))

    ## Check columns, and make sure we're returning just those columns that
    ## are requested.
    ResAll <- hostgenes(DB)
    Res <- hostgenes(DB, columns="gene_name")
    checkEquals(colnames(Res), "gene_name")
    checkTrue(nrow(ResAll) > nrow(Res))
    ## But still all of em have to be present.
    checkEquals(sort(unique(Res$gene_name)), sort(unique(ResAll$gene_name)))
    ## And now for columns that are not in host_gene
    Res <- hostgenes(DB, columns="mat_mirna_name")
    checkEquals(colnames(Res), "mat_mirna_name")
    ## hostgenesBy
    ResAll <- hostgenesBy(DB)
    Res <- hostgenesBy(DB, columns="gene_name")

}

test_hosttx <- function(){
    HT <- hosttx(DB)
    HT <- hosttx(DB, filter=list(TxBiotypeFilter("miRNA")))
    HT <- hosttxBy(DB, by="host_gene")

    ## Check columns, and make sure we're returning just those columns that
    ## are requested.
    ResAll <- hosttx(DB)
    Res <- hosttx(DB, columns="tx_biotype")
    checkEquals(colnames(Res), "tx_biotype")
    checkTrue(nrow(ResAll) > nrow(Res))
    ## But still all of em have to be present.
    checkEquals(sort(unique(Res$tx_biotype)), sort(unique(ResAll$tx_biotype)))
    ## And now for columns that are not in host_gene
    Res <- hosttx(DB, columns="mat_mirna_name")
    checkEquals(colnames(Res), "mat_mirna_name")
    Resmir <- matmirnas(DB, columns="mat_mirna_name")
    checkEquals(nrow(Res), nrow(Resmir))
    ## hosttxBy
    ResAll <- hosttxBy(DB)
    Res <- hosttxBy(DB, columns="tx_biotype")

}

test_probesets <- function(){
    if(DB@have_array_features){
        ## makes only sense if we have that functionality...
        PS <- probesets(DB)
        PM <- premirnasBy(DB, by="probeset")
        PS <- probesetsBy(DB, by="pre_mirna")

        ResAll <- probesets(DB)
        Res <- probesets(DB, columns="array_id")
        checkEquals(colnames(Res), "array_id")
        checkTrue(nrow(ResAll) > nrow(Res))
        ## But still all of em have to be present.
        checkEquals(sort(unique(Res$array_id)), sort(unique(ResAll$array_id)))
        ## And now for columns that are not in host_gene
        Res <- probesets(DB, columns="mat_mirna_name")
        checkEquals(colnames(Res), "mat_mirna_name")
        Resmir <- matmirnas(DB, columns="mat_mirna_name")
        checkEquals(nrow(Res), nrow(Resmir))
        ## hosttxBy
        ResAll <- probesetsBy(DB)
        Res <- probesetsBy(DB, columns="array_id")
    }
}

test_confidence_filter <- function(){
    if(DB@have_premirna_confidence){
        highPres <- premirnas(DB, filter=PreMirnaConfidence())
    }
    if(DB@have_matmirna_confidence){
        highMats <- matmirnas(DB, filter=MatMirnaConfidence())
    }
}

test_left_join <- function(){
    ## We're using left joins, thus, a call to premirnas including host tx
    ## will yield different results than a call to hosttx including premirnas.
    fromPre <- premirnas(DB, columns=c("pre_mirna_name", "tx_id"))
    fromTx <- hosttx(DB, columns=c("pre_mirna_name", "tx_id"))
    checkTrue(nrow(fromPre) > nrow(fromTx))
    checkTrue(sum(is.na(fromPre$tx_id)) > 0)
    checkTrue(sum(is.na(fromTx$tx_id)) == 0)
}


