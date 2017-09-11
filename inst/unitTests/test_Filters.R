## testing Filters defined in the package.
## detach("package:MirhostDb.Hsapiens.v75.v20", unload=TRUE)
## detach("package:mirhostgenes", unload=TRUE)

library("MirhostDb.Hsapiens.v75.v20")
MHDB <- MirhostDb.Hsapiens.v75.v20

## testing GeneIdFilter
test_PositionFilter <- function(){
    Filt <- PositionFilter("exonic")
    ## check if column matches the present database.
    checkEquals(column(Filt, MHDB), "host_tx.in_exon")
    Filt <- PositionFilter("intronic")
    ## check if column matches the present database.
    checkEquals(column(Filt, MHDB), "host_tx.in_intron")
    checkException(PositionFilter("a"))

    ## value
    checkEquals(value(Filt), "intronic")
    Filt <- PositionFilter("exonic")
    checkEquals(value(Filt), "exonic")

    ## something wrong; allow setting the value, but throw an error
    ## if we're calling it for MHDB.
    checkException(PositionFilter("something"))
    res <- where(Filt)
    res <- where(Filt, MHDB)
    ## condition
    checkEquals(condition(Filt), "==")
}

test_PreMirnaFilter <- function(){
    Filt <- PreMirnaFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_name")

    ## condition
    checkException(PreMirnaFilter("a", condition = ">"))
    condition(Filt)

    ## Value
    checkEquals(value(Filt), "a")
    Filt <- PreMirnaFilter(c("a", "b"))
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- PreMirnaFilter(c("a", "b"), "!=")
    checkEquals(condition(Filt), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    where(Filt, MHDB)
}

test_PreMirnaIdFilter <- function(){
    Filt <- PreMirnaIdFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_id")
    value(Filt)
    Filt <- PreMirnaIdFilter(c("a", "b"))
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- PreMirnaIdFilter(c("a", "b"), condition = "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    where(Filt, MHDB)
}

test_MatMirnaFilter <- function(){
    Filt <- MatMirnaFilter("a")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_name")
    checkEquals(value(Filt), "a")
    Filt <- MatMirnaFilter(c("a", "b"))
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- MatMirnaFilter(c("a", "b"), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    checkEquals(where(Filt, MHDB), "mat_mirna.mat_mirna_name not in ('a','b')")
}

test_MatMirnaIdFilter <- function(){
    Filt <- MatMirnaIdFilter("a")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_id")
    checkEquals(value(Filt), "a")
    Filt <- MatMirnaIdFilter(c("a", "b"))
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- MatMirnaIdFilter(c("a", "b"), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    checkEquals(where(Filt, MHDB), "mat_mirna.mat_mirna_id not in ('a','b')")
}

test_MirfamFilter <- function(){
    Filt <- MirfamFilter("a")
    checkEquals(column(Filt, MHDB), "mirfam.mirfam_name")
    checkEquals(value(Filt), "a")
    Filt <- MirfamFilter(c("a", "b"))
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- MirfamFilter(c("a", "b"), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    checkEquals(where(Filt, MHDB), "mirfam.mirfam_name not in ('a','b')")
}

test_MirfamIdFilter <- function(){
    Filt <- MirfamIdFilter("a")
    checkEquals(column(Filt, MHDB), "mirfam.mirfam_id")
    checkEquals(value(Filt), "a")
    Filt <- MirfamIdFilter(c("a", "b"))
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- MirfamIdFilter(c("a", "b"), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    checkEquals(where(Filt, MHDB), "mirfam.mirfam_id not in ('a','b')")
}

test_AlignmentIdFilter <- function(){
    Filt <- AlignmentIdFilter(1:10)
    checkEquals(condition(Filt), "==")
    checkEquals(where(Filt), "in ('1','2','3','4','5','6','7','8','9','10')")
    Filt <- AlignmentIdFilter(1:10, condition="!=")
    checkEquals(condition(Filt), "!=")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_algn_id")
}

## test_ArrayFilter <- function(){
##     Filt <- ArrayFilter("a")
##     checkEquals(column(Filt, MHDB), "array_features.array_id")
## }

## test_ProbesetIdFilter <- function(){
##     Filt <- ProbesetIdFilter("a")
##     checkEquals(column(Filt, MHDB), "array_features.feature_id")
## }

test_DatabaseFilter <- function(){
    Filt <- DatabaseFilter("a")
    checkEquals(column(Filt, MHDB), "host_gene.database")
    checkEquals(value(Filt), "a")
    Filt <- DatabaseFilter(c("a", "b"))
    checkEquals(where(Filt), "in ('a','b')")
    Filt <- DatabaseFilter(c("a", "b"), "!=")
    checkEquals(where(Filt), "not in ('a','b')")

    checkEquals(where(Filt, MHDB), "host_gene.database not in ('a','b')")
}

test_MatMirnaConfidence <- function(){
    Conf <- MatMirnaConfidence(value="other")
    ## we expect an error for where and MirhostDb.
    checkEquals(where(Conf), "='other'")
    checkException(where(Conf, MHDB))
    ## doing it the "right" way:
    Conf <- MatMirnaConfidence(value="high")
    checkEquals(where(Conf, MHDB), "mat_mirna.mat_mirna_confidence =1")

}

test_PreMirnaConfidence <- function(){
    Conf <- PreMirnaConfidence(value="other")
    ## we expect an error for where and MirhostDb.
    checkEquals(where(Conf), "='other'")
    checkException(where(Conf, MHDB))
    ## doing it the "right" way:
    Conf <- PreMirnaConfidence(value="high")
    checkEquals(where(Conf, MHDB), "pre_mirna.pre_mirna_confidence =1")
}

test_ReadCountFilter <- function(){
    rcf <- ReadCountFilter(value=10)
    checkEquals(where(rcf), ">10")
    ## throw an error
    checkException(ReadCountFilter(value="a"))
    checkEquals(where(rcf, MHDB), "pre_mirna.pre_mirna_read_count >10")
    checkEquals(column(rcf, MHDB), "pre_mirna.pre_mirna_read_count")

    ## use a not allowed value for of
    rcf <- ReadCountFilter(value=10, of="fasfd")
    checkException(where(rcf, MHDB))

    ## mat_mirna
    rcf <- ReadCountFilter(value="10", of="mat_mirna")
    checkEquals(where(rcf, MHDB), "mat_mirna.mat_mirna_read_count >10")
    checkEquals(column(rcf, MHDB), "mat_mirna.mat_mirna_read_count")

    checkEquals(value(rcf), 10L)
}


