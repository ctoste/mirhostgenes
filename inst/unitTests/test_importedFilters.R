## testing imported filters from the ensembldb package.
## GeneifFilter, GenebiotypeFilter, GenenameFilter, TxIdFilter, TxbiotypeFilter, ExonIdFilter, SeqNameFilter,
## SeqStartFilter, SeqEndFilter, SeqstrandFilter.

library("MirhostDb.Hsapiens.v75.v20")
MHDB <- MirhostDb.Hsapiens.v75.v20


## testing GeneIdFilter
test_GeneIdFilter <- function(){
    Filt <- GeneIdFilter("a")
    ## check if column matches the present database.
    checkEquals(column(Filt, MHDB), "host_tx.gene_id")
    ## check error if value is not as expected.
    checkException(GeneIdFilter("ENSG000001", ">"))
}

test_GeneBiotypeFilter <- function(){
    Filt <- GeneBiotypeFilter("protein_coding")
    checkEquals(column(Filt, MHDB), "host_gene.gene_biotype")
    checkException(GeneBiotypeFilter("protein_coding", ">"))
}

test_GenenameFilter <- function(){
    Filt <- GenenameFilter("genename")
    checkEquals(column(Filt, MHDB), "host_gene.gene_name")
    checkException(GenenameFilter("genename", ">"))
}

test_TxIdFilter <- function(){
    Filt <- TxIdFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.tx_id")
    checkException(TxIdFilter("a", ">"))
}

test_TxBiotypeFilter <- function(){
    Filt <- TxBiotypeFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.tx_biotype")
    checkException(TxBiotypeFilter("a", ">"))
}

test_ExonIdFilter <- function(){
    Filt <- ExonIdFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.exon_id")
    checkException(ExonIdFilter("a", ">"))
}

## SeqNameFilter
test_SeqNameFilter <- function(){
    Filt <- SeqNameFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.seq_name")
    checkException(SeqNameFilter("a", ">"))
}

## SeqstrandFilter
test_SeqstrandFilter <- function(){
    Filt <- SeqStrandFilter("-")
    checkEquals(column(Filt, MHDB), "pre_mirna.seq_strand")
}

## SeqStartFilter, feature
test_SeqStartFilter <- function(){
    Filt <- SeqStartFilter(123, feature="pre_mirna")
    checkEquals(column(Filt, MHDB),
                "pre_mirna.pre_mirna_seq_start")
    Filt <- SeqStartFilter(123, feature="mat_mirna")
    checkEquals(column(Filt, MHDB),
                "mat_mirna.mat_mirna_seq_start")
}

## SeqEndFilter
test_SeqEndFilter <- function(){
    Filt <- SeqEndFilter(123, feature="pre_mirna")
    checkEquals(column(Filt, MHDB),
                "pre_mirna.pre_mirna_seq_end")
    Filt <- SeqEndFilter(123, feature="mat_mirna")
    checkEquals(column(Filt, MHDB),
                "mat_mirna.mat_mirna_seq_end")
}

