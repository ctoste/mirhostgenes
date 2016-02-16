## testing Filters defined in the package.
## detach("package:MirhostDb.Hsapiens.v75.v20", unload=TRUE)
## detach("package:mirhostgenes", unload=TRUE)

library("MirhostDb.Hsapiens.v75.v20")
MHDB <- MirhostDb.Hsapiens.v75.v20

## testing GeneidFilter
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
    value(Filt) <- "exonic"
    checkEquals(value(Filt), "exonic")
    ## something wrong; allow setting the value, but throw an error
    ## if we're calling it for MHDB.
    value(Filt) <- "something"
    checkEquals(value(Filt), "something")
    checkException(where(Filt, MHDB))
    ## condition
    condition(Filt)
    condition(Filt) <- "="
}

test_PremirnaFilter <- function(){
    Filt <- PremirnaFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_name")

    ## condition
    checkException(condition(Filt) <- ">")
    condition(Filt)

    ## Value
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_PremirnaidFilter <- function(){
    Filt <- PremirnaidFilter("a")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_id")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_MatmirnaFilter <- function(){
    Filt <- MatmirnaFilter("a")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_name")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_MatmirnaidFilter <- function(){
    Filt <- MatmirnaidFilter("a")
    checkEquals(column(Filt, MHDB), "mat_mirna.mat_mirna_id")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_MirfamFilter <- function(){
    Filt <- MirfamFilter("a")
    checkEquals(column(Filt, MHDB), "mirfam.mirfam_name")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_MirfamidFilter <- function(){
    Filt <- MirfamidFilter("a")
    checkEquals(column(Filt, MHDB), "mirfam.mirfam_id")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_AlignmentidFilter <- function(){
    Filt <- AlignmentidFilter(1:10)
    checkEquals(condition(Filt), "in")
    Filt <- AlignmentidFilter(1:10, condition="!=")
    checkEquals(condition(Filt), "not in")
    checkEquals(column(Filt, MHDB), "pre_mirna.pre_mirna_algn_id")
}

## test_ArrayFilter <- function(){
##     Filt <- ArrayFilter("a")
##     checkEquals(column(Filt, MHDB), "array_features.array_id")
## }

## test_ProbesetidFilter <- function(){
##     Filt <- ProbesetidFilter("a")
##     checkEquals(column(Filt, MHDB), "array_features.feature_id")
## }

test_DatabaseFilter <- function(){
    Filt <- DatabaseFilter("a")
    checkEquals(column(Filt, MHDB), "host_gene.database")
    value(Filt)
    value(Filt) <- c("a", "b")
    checkEquals(condition(Filt), "in")
    condition(Filt) <- "!="
    Filt
    checkEquals(condition(Filt), "not in")

    where(Filt, MHDB)
}

test_MatmirnaConfidence <- function(){
    Conf <- MatmirnaConfidence(value="other")
    ## we expect an error for where and MirhostDb.
    checkEquals(where(Conf), "= 'other'")
    checkException(where(Conf, MHDB))
    ## doing it the "right" way:
    Conf <- MatmirnaConfidence(value="high")
    checkEquals(where(Conf, MHDB), "mat_mirna.mat_mirna_confidence =1")

    value(Conf) <- "agkgnfk"
    checkException(where(Conf, MHDB))
}

test_PremirnaConfidence <- function(){
    Conf <- PremirnaConfidence(value="other")
    ## we expect an error for where and MirhostDb.
    checkEquals(where(Conf), "= 'other'")
    checkException(where(Conf, MHDB))
    ## doing it the "right" way:
    Conf <- PremirnaConfidence(value="high")
    checkEquals(where(Conf, MHDB), "pre_mirna.pre_mirna_confidence =1")
}

test_ReadCountFilter <- function(){
    rcf <- ReadCountFilter(value=10)
    checkEquals(where(rcf), "> 10")
    ## throw an error
    checkException(ReadCountFilter(value="a"))
    checkEquals(where(rcf, MHDB), "pre_mirna.pre_mirna_read_count > 10")
    checkEquals(column(rcf, MHDB), "pre_mirna.pre_mirna_read_count")

    ## use a not allowed value for of
    rcf <- ReadCountFilter(value=10, of="fasfd")
    checkException(where(rcf, MHDB))

    ## mat_mirna
    rcf <- ReadCountFilter(value="10", of="mat_mirna")
    checkEquals(where(rcf, MHDB), "mat_mirna.mat_mirna_read_count > 10")
    checkEquals(column(rcf, MHDB), "mat_mirna.mat_mirna_read_count")

    value(rcf)
    value(rcf) <- 20

    checkException(value(rcf) <- "dfdsf")
}


