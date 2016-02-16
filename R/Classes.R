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
setClass("AlignmentidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
AlignmentidFilter <- function(value, condition="="){
    ## here value could be intronic, exonic, both.
    ## condition will not be considered...
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("AlignmentidFilter", condition=condition, value=as.character(value)))
}

## Table array_feature
## filter for attribute array_id
setClass("ArrayFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
ArrayFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("ArrayFilter", condition=condition, value=as.character(value)))
}
## filter for attribute array_id
setClass("ProbesetidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
ProbesetidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("ProbesetidFilter", condition=condition, value=as.character(value)))
}

## Table host_gene
## filter for the database.
setClass("DatabaseFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
DatabaseFilter <- function(value, condition="="){
    ## here value could be intronic, exonic, both.
    ## condition will not be considered...
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("DatabaseFilter", condition=condition, value=as.character(value)))
}


## Table host_tx
## filter for position
setClass("PositionFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
PositionFilter <- function(value, condition="="){
    ## here value could be intronic, exonic, both.
    ## condition will not be considered...
    value <- match.arg(value, c("intronic", "exonic", "both"))
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("PositionFilter", condition=condition, value=as.character(value)))
}

##***********************************************************************
## Table pre_mirna
## filter for pre_mirnas
setClass("PremirnaFilter", contains="BasicFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE,
             match.case=TRUE
         )
         )
PremirnaFilter <- function(value, condition="=", match.case=TRUE){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("PremirnaFilter", match.case=match.case, condition=condition,
               value=as.character(value)))
}
setClass("PremirnaidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
PremirnaidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("PremirnaidFilter", condition=condition, value=as.character(value)))
}

##***********************************************************************
## Table mat_mirna
## filter for mat_mirnas
setClass("MatmirnaFilter", contains="BasicFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE,
             match.case=TRUE
         )
         )
MatmirnaFilter <- function(value, condition="=", match.case=TRUE){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("MatmirnaFilter", match.case=match.case, condition=condition,
               value=as.character(value)))
}
setClass("MatmirnaidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
MatmirnaidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("MatmirnaidFilter", condition=condition, value=as.character(value)))
}


##***********************************************************************
## Table mat_mirna
## filter for mat_mirnas
setClass("MirfamFilter", contains="BasicFilter",
         representation(
             match.case="logical"
         ),
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE,
             match.case=TRUE
         )
         )
MirfamFilter <- function(value, condition="=", match.case=TRUE){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("MirfamFilter", match.case=match.case, condition=condition,
               value=as.character(value)))
}
setClass("MirfamidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
         )
         )
MirfamidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("MirfamidFilter", condition=condition, value=as.character(value)))
}

## special type fo filter: high confidence filter.
setClass("MatmirnaConfidence", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="high",
             .valueIsCharacter=TRUE
         ))
MatmirnaConfidence <- function(value="high", condition="="){
    if(missing(value)){
        value <- "high"
    }
    if(length(value) > 1){
        value <- value[1]
        warning("MatmirnaConfidence filter does only support a single value. Taking value[1].")
    }
    if(!(condition %in% c("=", "!=")))
        stop("MatmirnaConfidence filter does not support a condition other than '=' or '!='.")
    return(new("MatmirnaConfidence", condition=condition,
               value=as.character(value)))
}

## special type fo filter: high confidence filter.
setClass("PremirnaConfidence", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="high",
             .valueIsCharacter=TRUE
         ))
PremirnaConfidence <- function(value="high", condition="="){
    if(missing(value)){
        value <- "high"
    }
    if(length(value) > 1){
        value <- value[1]
        warning("PremirnaConfidence filter does only support a single value. Taking value[1].")
    }
    if(!(condition %in% c("=", "!=")))
        stop("PremirnaConfidence filter does not support a condition other than '=' or '!='.")
    return(new("PremirnaConfidence", condition=condition,
               value=as.character(value)))
}

setClass("ReadCountFilter", contains="BasicFilter",
         representation(
             of="character"
         ),
         prototype=list(
             condition=">",
             value="0",
             .valueIsCharacter=FALSE,
             of="pre_mirna"
         ))
ReadCountFilter <- function(value=0, condition=">", of="pre_mirna"){
    if(length(value) > 1){
        value <- value[1]
        warning("ReadCountFilter does only support a single value! Taking value[1].")
    }
    if(is.na(as.numeric(value)))
        stop("value has to be numeric for ReadCountFilter!")
    return(new("ReadCountFilter", condition=condition, value=as.character(value), of=of))
}
## do: ProbesetFilter, throws error if no required table there
## do: ArrayFilter, throws error if no required table there


