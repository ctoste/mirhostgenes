.MATMIRNA.ATTRS <- c("mat_mirna_id", "mat_mirna_name", "pre_mirna_algn_id",
                     "mat_mirna_seq_start", "mat_mirna_seq_end",
                     "mat_mirna_confidence", "mat_mirna_read_count",
                     "mat_mirna_experiment_count")
.PREMIRNA.ATTRS <- c("pre_mirna_algn_id", "pre_mirna_id", "pre_mirna_name",
                     "seq_name", "seq_strand", "pre_mirna_seq_start",
                     "pre_mirna_seq_end", "pre_mirna_confidence",
                     "pre_mirna_read_count", "pre_mirna_experiment_count")
.HOSTTX.ATTRS <- c("pre_mirna_algn_id", "tx_id", "tx_biotype", "in_intron",
                   "in_exon", "is_outside", "exon_id", "gene_id")
.HOSTGENE.ATTRS <- c("gene_id", "gene_biotype", "gene_name", "entrezid", "database")


##***********************************************************************
##
##     MirhostDb
##
##     That's the main class with the connection to the database that
##     provides also all functionality to work with the database.
##
##***********************************************************************
setClass("MirhostDb",
         representation(con="DBIConnection",
                        tables="list",
                        have_pre_sequence="logical",
                        have_mirfam="logical",
                        have_array_features="logical",
                        have_premirna_confidence="logical",
                        have_premirna_readcount="logical",
                        have_matmirna_confidence="logical",
                        have_matmirna_readcount="logical"
                        ),
         prototype=list(con=NULL,
                        tables=list(),
                        have_pre_sequence=FALSE,
                        have_mirfam=FALSE,
                        have_array_features=FALSE,
                        have_premirna_confidence=FALSE,
                        have_premirna_readcount=FALSE,
                        have_matmirna_confidence=FALSE,
                        have_matmirna_readcount=FALSE
                        )
         )



##***********************************************************************
##
##     Filter classes
##
##     All filter classes from ensembldb can be used in this package.
##
##     We're defining there in addition:
##     PositionFilter
##     PremirnaidFilter
##     PremirnanameFilter
##     MatmirnaidFilter
##     MatmirnanameFilter
##
##***********************************************************************
## Table pre_mirna
## filter for the alignment id.
setClass("AlignmentIdFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "pre_mirna_algn_id"
         )
         )
AlignmentIdFilter <- function(value, condition = "==") {
    ## here value could be intronic, exonic, both.
    ## condition will not be considered...
    new("AlignmentIdFilter", condition = condition, value = as.character(value))
}

## Table array_feature
## filter for attribute array_id
setClass("ArrayFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "array_id"
         )
         )
ArrayFilter <- function(value, condition = "=="){
    new("ArrayFilter", condition=condition, value=as.character(value))
}
## filter for attribute array_id
setClass("ProbesetIdFilter", contains="CharacterFilter",
         prototype = list(
             condition = "==",
             value="",
             field = "probeset_id"
         )
         )
ProbesetIdFilter <- function(value, condition = "=="){
    new("ProbesetIdFilter", condition = condition, value=as.character(value))
}

## Table host_gene
## filter for the database.
setClass("DatabaseFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "database"
         )
         )
DatabaseFilter <- function(value, condition="=="){
    new("DatabaseFilter", condition=condition, value=as.character(value))
}


## Table host_tx
## filter for position
setClass("PositionFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = ""
         )
         )
PositionFilter <- function(value, condition="=="){
    ## here value could be intronic, exonic, both.
    ## condition will not be considered...
    value <- match.arg(value, c("intronic", "exonic", "both"))
    new("PositionFilter", condition=condition, value=as.character(value))
}

##***********************************************************************
## Table pre_mirna
## filter for pre_mirnas
setClass("PreMirnaFilter", contains="CharacterFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="==",
             value="",
             field = "pre_mirna_name",
             match.case=TRUE
         )
         )
PreMirnaFilter <- function(value, condition="==", match.case=TRUE){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    new("PreMirnaFilter", match.case = match.case, condition = condition,
               value = as.character(value))
}
setClass("PreMirnaIdFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "pre_mirna_id",
             .valueIsCharacter=TRUE
         )
         )
PreMirnaIdFilter <- function(value, condition="=="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    new("PreMirnaIdFilter", condition=condition, value=as.character(value))
}

##***********************************************************************
## Table mat_mirna
## filter for mat_mirnas
setClass("MatMirnaFilter", contains="CharacterFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="==",
             value="",
             field = "mat_mirna_name",
             match.case=TRUE
         )
         )
MatMirnaFilter <- function(value, condition="==", match.case=TRUE){
    new("MatMirnaFilter", match.case=match.case, condition=condition,
        value=as.character(value))
}
setClass("MatMirnaIdFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "mat_mirna_id",
             .valueIsCharacter=TRUE
         )
         )
MatMirnaIdFilter <- function(value, condition="=="){
    new("MatMirnaIdFilter", condition=condition, value=as.character(value))
}


##***********************************************************************
## Table mat_mirna
## filter for mat_mirnas
setClass("MirfamFilter", contains="CharacterFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="==",
             value="",
             field = "mirfam_name",
             match.case=TRUE
         )
         )
MirfamFilter <- function(value, condition="==", match.case=TRUE){
    new("MirfamFilter", match.case=match.case, condition=condition,
        value=as.character(value))
}
setClass("MirfamIdFilter", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="",
             field = "mirfam_id"
         )
         )
MirfamIdFilter <- function(value, condition="=="){
    new("MirfamIdFilter", condition=condition, value=as.character(value))
}

## special type fo filter: high confidence filter.
setClass("MatMirnaConfidence", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="high",
             field = "mat_mirna_confidence"
         ))
MatMirnaConfidence <- function(value="high", condition="=="){
    if(length(value) > 1){
        value <- value[1]
        warning("MatMirnaConfidence filter does only support a single value. Taking value[1].")
    }
    if(!(condition %in% c("==", "!=")))
        stop("MatMirnaConfidence filter does not support a condition other than '=' or '!='.")
    return(new("MatMirnaConfidence", condition=condition,
               value=as.character(value)))
}

## special type fo filter: high confidence filter.
setClass("PreMirnaConfidence", contains="CharacterFilter",
         prototype=list(
             condition="==",
             value="high",
             field = "pre_mirna_confidence"
         ))
PreMirnaConfidence <- function(value="high", condition="=="){
    if(missing(value)){
        value <- "high"
    }
    if(length(value) > 1){
        value <- value[1]
        warning("PreMirnaConfidence filter does only support a single value. Taking value[1].")
    }
    if(!(condition %in% c("==", "!=")))
        stop("PreMirnaConfidence filter does not support a condition other than '=' or '!='.")
    new("PreMirnaConfidence", condition=condition,
        value=as.character(value))
}

setClass("ReadCountFilter", contains="IntegerFilter",
         representation(
             of="character"
         ),
         prototype=list(
             condition = ">",
             value = 0L,
             of = "pre_mirna"
         ))
ReadCountFilter <- function(value=0, condition=">", of="pre_mirna"){
    if(length(value) > 1){
        value <- value[1]
        warning("ReadCountFilter does only support a single value! Taking value[1].")
    }
    if(is.na(as.numeric(value)))
        stop("value has to be numeric for ReadCountFilter!")
    new("ReadCountFilter", condition=condition, value=as.integer(value), of=of)
}
## do: ProbesetFilter, throws error if no required table there
## do: ArrayFilter, throws error if no required table there

setClass("SeqEndFilter", contains="IntegerFilter",
         representation(
             feature="character"
         ),
         prototype=list(
             condition = ">",
             feature = "gene",
             value = 0L
         ))
SeqEndFilter <- function(value = 0, condition = ">", feature = "pre_mirna"){
    if(length(value) > 1){
        value <- value[1]
        warning("SeqEndFilter does only support a single value! Taking value[1].")
    }
    new("SeqEndFilter", condition = condition, value = as.integer(value),
        feature = feature)
}
setClass("SeqStartFilter", contains="IntegerFilter",
         representation(
             feature="character"
         ),
         prototype=list(
             condition = ">",
             feature = "gene",
             value = 0L
         ))
SeqStartFilter <- function(value = 0, condition = ">", feature = "pre_mirna"){
    if(length(value) > 1){
        value <- value[1]
        warning("SeqStartFilter does only support a single value! Taking value[1].")
    }
    new("SeqStartFilter", condition = condition, value = as.integer(value),
        feature = feature)
}


