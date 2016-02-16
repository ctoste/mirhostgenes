## testing imported filters from the ensembldb package.
## GeneifFilter, GenebiotypeFilter, GenenameFilter, TxidFilter, TxbiotypeFilter, ExonidFilter, SeqnameFilter,
## SeqstartFilter, SeqendFilter, SeqstrandFilter.

library("MirhostDb.Hsapiens.v75.v20")
MHDB <- MirhostDb.Hsapiens.v75.v20


## testing GeneidFilter
test_GeneidFilter <- function(){
    Filt <- GeneidFilter("a")
    ## check if column matches the present database.
    checkEquals(column(Filt, MHDB), "host_tx.gene_id")
    ## check error if value is not as expected.
    checkException(GeneidFilter("ENSG000001", ">"))
}

test_GenebiotypeFilter <- function(){
    Filt <- GenebiotypeFilter("protein_coding")
    checkEquals(column(Filt, MHDB), "host_gene.gene_biotype")
    checkException(GenebiotypeFilter("protein_coding", ">"))
}

test_GenenameFilter <- function(){
    Filt <- GenenameFilter("genename")
    checkEquals(column(Filt, MHDB), "host_gene.gene_name")
    checkException(GenenameFilter("genename", ">"))
}

test_TxidFilter <- function(){
    Filt <- TxidFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.tx_id")
    checkException(TxidFilter("a", ">"))
}

test_TxbiotypeFilter <- function(){
    Filt <- TxbiotypeFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.tx_biotype")
    checkException(TxbiotypeFilter("a", ">"))
}

test_ExonidFilter <- function(){
    Filt <- ExonidFilter("a")
    checkEquals(column(Filt, MHDB), "host_tx.exon_id")
    checkException(ExonidFilter("a", ">"))
}

## SeqnameFilter
test_SeqnameFilter <- function(){
    Filt <- SeqnameFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.seq_name")
    checkException(SeqnameFilter("a", ">"))
}

## SeqstrandFilter
test_SeqstrandFilter <- function(){
    checkException(SeqstrandFilter("a"))
    Filt <- SeqstrandFilter("-")
    checkEquals(column(Filt, MHDB), "pre_mirna.seq_strand")
}

## SeqstartFilter, feature
test_SeqstartFilter <- function(){
    Filt <- SeqstartFilter(123, feature="pre_mirna")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_seq_start")
    Filt <- SeqstartFilter(123, feature="mat_mirna")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_seq_start")
}

## SeqendFilter
test_SeqendFilter <- function(){
    Filt <- SeqendFilter(123, feature="pre_mirna")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_seq_end")
    Filt <- SeqendFilter(123, feature="mat_mirna")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_seq_end")
}

