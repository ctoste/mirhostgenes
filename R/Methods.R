###
## Methods and their implementation for classes defined in Classes.R





##***********************************************************************
##
##     Implementations for MirhostDb
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("dbconn", "MirhostDb", function(x){
    con <- x@con
    return(con)
})

validateMirnahostgenesDb <- function(object){
    ## check if the database contains all required tables...
    if(!is.null(object@con)){
        have <- dbListTables(object@con)
        need <- c("host_gene", "host_tx", "metadata", "mat_mirna", "pre_mirna", "chromosome")
        notthere <- need[ !(need %in% have) ]
        if(length(notthere) > 0){
            return(paste("Required tables", notthere, "not found in the database!"))
        }
        ## check the attributes in the tables:
        ## mat_mirna
        res <- dbGetQuery(object@con, "select * from mat_mirna limit 1")
        notthere <- !(.MATMIRNA.ATTRS %in% colnames(res))
        if(any(notthere)){
            if(all(.MATMIRNA.ATTRS[notthere] %in% c("mat_mirna_confidence",
                                                    "mat_mirna_read_count",
                                                    "mat_mirna_experiment_count"))){
                if(any(.MATMIRNA.ATTRS[notthere] == "mat_mirna_confidence"))
                    warning("No confidence information for mature miRNAs available.")
                if(any(.MATMIRNA.ATTRS[notthere] == "mat_mirna_read_count"))
                    warning("No mature miRNA read count information available.")
            }else{
                return(paste("Required database columns",
                             paste(.MATMIRNA.ATTRS[ notthere ], collapse=", "),
                             "not present in database table mat_mirna!"))
            }
        }
        ## pre_mirna
        res <- dbGetQuery(object@con, "select * from pre_mirna limit 1")
        notthere <- !(.PREMIRNA.ATTRS %in% colnames(res))
        if(any(notthere)){
            if(all(.PREMIRNA.ATTRS[notthere] %in% c("pre_mirna_confidence",
                                                    "pre_mirna_read_count",
                                                    "pre_mirna_experiment_count"))){
                if(any(.PREMIRNA.ATTRS[notthere] == "pre_mirna_confidence"))
                    warning("No confidence information for pre-miRNAs available.")
                if(any(.PREMIRNA.ATTRS[notthere] == "pre_mirna_read_count"))
                    warning("No pre-miRNA read count information available.")
            }else{
                return(paste("Required database columns",
                             paste(.PREMIRNA.ATTRS[ notthere ], collapse=", "),
                             "not present in database table pre_mirna!"))
            }
        }
        ## host_tx
        res <- dbGetQuery(object@con, "select * from host_tx limit 1")
        notthere <- !(.HOSTTX.ATTRS %in% colnames(res))
        if(any(notthere))
            return(paste("Required database columns",
                         paste(.HOSTTX.ATTRS[ notthere ], collapse=", "),
                         "not present in database table host_tx!"))
        ## host_gene
        res <- dbGetQuery(object@con, "select * from host_gene limit 1")
        notthere <- !(.HOSTGENE.ATTRS %in% colnames(res))
        if(any(notthere))
            return(paste("Required database columns",
                         paste(.HOSTGENE.ATTRS[ notthere ], collapse=", "),
                         "not present in database table host_gene!"))

    }
    return(TRUE)
}
setValidity("MirhostDb", validateMirnahostgenesDb)
setMethod("initialize", "MirhostDb", function(.Object,...){
    OK <- validateMirnahostgenesDb(.Object)
    if(class(OK)=="character"){
        stop(OK)
    }
    callNextMethod(.Object, ...)
})

setMethod("seqinfo", "MirhostDb", function(x){
    Chrs <- dbGetQuery(dbconn(x), "select * from chromosome")
    Chr.build <- .getMetaDataValue(dbconn(x), "genome_build")
    SI <- Seqinfo(seqnames=Chrs$seq_name,
                  seqlengths=Chrs$seq_length,
                  isCircular=Chrs$is_circular==1, genome=Chr.build)
    return(SI)
})

setMethod("metadata", "MirhostDb", function(x){
    return(dbGetQuery(dbconn(x), "select * from metadata"))
})

setMethod("ensemblVersion", "MirhostDb", function(x){
    return(version(x, "ensembl"))
})

setMethod("mirbaseVersion", "MirhostDb", function(x){
    return(version(x, "mirbase"))
})

setMethod("version", "MirhostDb", function(object, what="ensembl", ...){
    what <- match.arg(what, c("ensembl", "mirbase"))
    return(.getMetaDataValue(object@con, paste0(what, "_version")))
})

setMethod("show", "MirhostDb", function(object){
    if(is.null(object@con)){
        cat("Dash it! Got an empty thing!\n")
    }else{
        con <- object@con
        info <- dbGetQuery(con, "select * from metadata;")
        cat("MirhostDb:\n")
        for(i in 1:nrow(info)){
            if(!(info[ i, "name" ] %in% c("did_core", "did_cdna", "did_otherfeatures",
                                          "did_vega", "prop_probes", "min_probe_algn",
                                          "max_mm")))
            {
                cat(paste0("| ", info[ i, "name" ], ": ", info[ i, "value" ], "\n"))
            }
        }
        ## info string on the databases
        Str <- "| Queried databases: "
        DBs <- NULL
        if(info[ info$name=="did_core", "value" ]==1)
            DBs <- c(DBs, "core")
        if(info[ info$name=="did_cdna", "value" ]==1)
            DBs <- c(DBs, "cdna")
        if(info[ info$name=="did_otherfeatures", "value" ]==1)
            DBs <- c(DBs, "otherfeatures")
        if(info[ info$name=="did_vega", "value" ]==1)
            DBs <- c(DBs, "vega")
        cat(paste0(Str, paste0(DBs, collapse=", "), "\n"))
        Tables <- listTables(object)
        if(any(names(Tables)=="pre_mirna_sequence"))
            cat("| Have pre-miRNA sequences\n")
        if(any(names(Tables)=="mirfam"))
            cat("| Have miRNA family definitions\n")
        if(any(names(Tables)=="array_feature")){
            cat("| Have microarray features:\n")
            cat(paste0("| - minimum required proportion of probes mapping a transcript: ",
                       info[ info$name=="prop_probes", 2 ], "\n"))
            cat(paste0("| - maximum allowed mismatches in probe alignments: ",
                       info[ info$name=="max_mm", 2 ], "\n"))
            cat(paste0("| - required mininal alignment length for a probe: ",
                       info[ info$name=="min_probe_algn", 2 ], "\n"))
        }
        ## summary on the host genes...
        Res <- dbGetQuery(object@con, "select count(distinct pre_mirna_id) from pre_mirna")
        cat(paste0("| Total number of pre-miRNAs: ",  Res[ 1, 1 ], "\n"))
        Res <- dbGetQuery(object@con,
                          paste0("select count(distinct pre_mirna_id) from pre_mirna join",
                                 " host_tx on (pre_mirna.pre_mirna_algn_id",
                                 "=host_tx.pre_mirna_algn_id) where tx_biotype!='miRNA'"))
        cat(paste0("| Number of pre-miRNAs for which host gene(s) are defined: ",
                   Res[ 1, 1 ], "\n"))
    }
})

### tables
## returns a named list with database table attributes
setMethod("listTables", "MirhostDb", function(x, ...){
    if(length(x@tables)==0){
        tables <- dbListTables(x@con)
        ## read the attributes for these tables.
        Tables <- vector(length=length(tables), "list")
        for(i in 1:length(Tables)){
            Tables[[ i ]] <- colnames(dbGetQuery(x@con, paste0("select * from ",
                                                               tables[ i ], " limit 1")))
        }
        names(Tables) <- tables
        x@tables <- Tables
    }
    Tab <- x@tables
    Tab <- Tab[ tablesByDegree(x, tab=names(Tab)) ]
    return(Tab)
})

setMethod("organism", "MirhostDb", function(object){
    org <- .getMetaDataValue(dbconn(object), "Organism")
    return(.cleanOrganismName(org))
})


### listColumns
## lists all attributes.
setMethod("listColumns", "MirhostDb", function(x, table, skip.keys=TRUE, ...){
    if(length(x@tables)==0){
        tables <- dbListTables(dbconn(x))
        ## read the attributes for these tables.
        Tables <- vector(length=length(tables), "list")
        for(i in 1:length(Tables)){
            Tables[[ i ]] <- colnames(dbGetQuery(dbconn(x),
                                                 paste0("select * from ", tables[ i ],
                                                        " limit 1")))
        }
        names(Tables) <- tables
        x@tables <- Tables
    }
    if(!missing(table)){
        attrs <- x@tables[[ table ]]
    }else{
        attrs <- unlist(x@tables, use.names=FALSE)
    }
    if(skip.keys){
        ## remove everything that has a _pk or _fk...
        idx <- grep(attrs, pattern="_fk$")
        if(length(idx) > 0)
            attrs <- attrs[ -idx ]
        idx <- grep(attrs, pattern="_pk$")
        if(length(idx) > 0)
            attrs <- attrs[ -idx ]
    }
    return(attrs)
})

setMethod("listArrays", "MirhostDb", function(x, ...){
    if(!x@have_array_features)
        stop("The database does not provide microarray features!")
    Res <- dbGetQuery(dbconn(x), "select distinct array_id from array_feature;")
    return(Res[ , "array_id" ])
})

setMethod("listDatabases", "MirhostDb", function(x, ...){
    Res <- dbGetQuery(dbconn(x), "select distinct database from host_gene;")
    return(Res[ , "database" ])
})

setMethod("listGenebiotypes", "MirhostDb", function(x, ...){
    Res <- dbGetQuery(dbconn(x), "select distinct gene_biotype from host_gene;")
    return(Res[, 1])
})
setMethod("listTxbiotypes", "MirhostDb", function(x, ...){
    Res <- dbGetQuery(dbconn(x), "select distinct tx_biotype from host_tx;")
    return(Res[, 1])
})



##************************************************************************
##
##   The main methods
##
##************************************************************************
### matmirnas
## get mature_mirnas from the database
setMethod("matmirnas", "MirhostDb", function(x,
                                             columns=listColumns(x, "mat_mirna"),
                                             filter, order.by="mat_mirna_id",
                                             order.type="asc",
                                             return.type="DataFrame"){
    return.type <- match.arg(return.type, c("DataFrame", "data.frame", "GRanges"))
    if(length(columns) == 0)
        stop("Length of argument 'columns' is 0!")
    ## attrs <- unique(c(columns, "mat_mirna_id"))
    attrs <- columns
    if(missing(filter))
        filter <- list()
    filter <- checkIsBasicFilter(filter)
    ## if we're going to return a GRanges, we have to make sure we get the required
    ## columns.
    if(return.type=="GRanges"){
        attrs <- unique(c(attrs, c("seq_name", "seq_strand", "mat_mirna_seq_start",
                                   "mat_mirna_seq_end")))
    }
    ## silently remove order.by=mat_mirna_id if it's not in attrs...
    if(!any(attrs==order.by))
        order.by <- ""
    if(missing(filter))
        filter <- list()
    if(any(attrs=="sequence")){
        ## if we've got sequence as a column, we want in reality to get the sequence of
        ## the mature miRNA. thus, we need mat_mirna_seq_start, mat_mirna_seq_end,
        ## pre_mirna_seq_start, pre_mirna_seq_end also.
        addAttrs <- c("mat_mirna_seq_start", "mat_mirna_seq_end", "pre_mirna_seq_start",
                      "pre_mirna_seq_end", "seq_strand")
        addAttrs <- addAttrs[!(addAttrs %in% attrs)]
    }else{
        addAttrs <- character()
    }
    Res <- .getWhat(x=x, columns=c(attrs, addAttrs), filter=filter,
                    order.by=order.by, order.type=order.type, join="left join",
                    start.table="mat_mirna")
    ## OK, I need to fix the sequence, if present...
    if(any(attrs=="sequence")){
        ## what sucks is that the pre-miRNA is NOT the genomic miRNA sequence,
        ## but the sequence "as.is".
        PreSeqs <- Res$sequence
        relStart <- Res$mat_mirna_seq_start - Res$pre_mirna_seq_start
        relEnd <- Res$mat_mirna_seq_end - Res$pre_mirna_seq_start
        ## for seq_strand == -1 I have to fix that!
        if(any(Res$seq_strand < 0)){
            minusStr <- Res$seq_strand < 0
            relStart[minusStr] <- Res$pre_mirna_seq_end[minusStr] -
                Res$mat_mirna_seq_end[minusStr]
            relEnd[minusStr] <- Res$pre_mirna_seq_end[minusStr] -
                Res$mat_mirna_seq_start[minusStr]
        }
        MatSeqs <- substring(PreSeqs, first=(relStart+1),
                             last=(relEnd+1))
        Res$sequence <- MatSeqs
        Res <- Res[, !(colnames(Res) %in% addAttrs)]
    }

    if(return.type=="DataFrame"){
        return(DataFrame(Res))
    }
    if(return.type=="GRanges"){
        meta.attrs <- attrs[ !(attrs %in% c("seq_name",
                                            "seq_strand",
                                            "mat_mirna_seq_start",
                                            "mat_mirna_seq_end")) ]
        SI <- seqinfo(x)
        SI <- SI[ unique(Res$seq_name) ]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$mat_mirna_seq_start,
                                     end=Res$mat_mirna_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI,
                      Res[ , meta.attrs, drop=FALSE ]
                      )
        return(GR)
    }
    return(Res)
})


### premirnas
## get pre_mirnas from the database
setMethod("premirnas", "MirhostDb", function(x,
                                             columns=listColumns(x, "pre_mirna"),
                                             filter,
                                             order.by="pre_mirna_id",
                                             order.type="asc",
                                             return.type="DataFrame"){
    return.type <- match.arg(return.type, c("DataFrame", "data.frame", "GRanges"))
    if(length(columns) == 0)
        stop("Length of argument 'columns' is 0!")
    attrs <- unique(columns)
    ## FIXME: Do I really have to add pre_mirna_id?
    ## attrs <- unique(c(columns, "pre_mirna_id"))
    if(missing(filter))
        filter <- list()
    filter <- checkIsBasicFilter(filter)
    ## silently remove order.by=mat_mirna_id if it's not in attrs...
    if(!any(attrs==order.by))
        order.by <- ""
    if(return.type=="GRanges"){
        attrs <- unique(c(attrs, c("seq_name", "seq_strand", "pre_mirna_seq_start",
                                   "pre_mirna_seq_end")))
    }
    Res <- .getWhat(x=x, columns=attrs, filter=filter,
                    order.by=order.by, order.type=order.type,
                    join="left join", start.table="pre_mirna")
    if(return.type=="DataFrame"){
        return(DataFrame(Res))
    }
    if(return.type=="GRanges"){
        meta.attrs <- attrs[ !(attrs %in% c("seq_name",
                                            "seq_strand",
                                            "pre_mirna_seq_start",
                                            "pre_mirna_seq_end")) ]
        SI <- seqinfo(x)
        SI <- SI[ unique(Res$seq_name) ]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$pre_mirna_seq_start,
                                     end=Res$pre_mirna_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI,
                      Res[ , meta.attrs, drop=FALSE ]
                      )
        return(GR)
    }
    return(Res)
})


### hosttx
## get host transcripts from the database
setMethod("hosttx", "MirhostDb", function(x,
                                          columns=listColumns(x, "host_tx"),
                                          filter,
                                          order.by="tx_id",
                                          order.type="asc",
                                          return.type="DataFrame"){
    return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
    if(length(columns) == 0)
        stop("Length of argument 'columns' is 0!")
    attrs <- unique(columns)
    ## attrs <- unique(c(columns, "tx_id"))
    if(missing(filter))
        filter <- list()
    filter <- checkIsBasicFilter(filter)
    ## silently remove order.by=mat_mirna_id if it's not in attrs...
    if(!any(attrs==order.by))
        order.by <- ""
    Res <- .getWhat(x=x, columns=attrs, filter=filter,
                    order.by=order.by, order.type=order.type,
                    join="left join", start.table="host_tx")
    if(return.type=="DataFrame"){
        Res <- DataFrame(Res)
    }
    return(Res)
})


### hostgenes
## get host genes from the database
setMethod("hostgenes", "MirhostDb",
          function(x,
                   columns=listColumns(x, "host_gene"),
                   filter,
                   order.by="gene_id",
                   order.type="asc",
                   return.type="DataFrame"){
              return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              attrs <- unique(columns)
              ## attrs <- unique(c(columns, "gene_id"))
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              ## silently remove order.by=mat_mirna_id if it's not in attrs...
              if(!any(attrs==order.by))
                  order.by <- ""
              Res <- .getWhat(x=x, columns=attrs, filter=filter,
                              order.by=order.by, order.type=order.type,
                              join="left join", start.table="host_gene")
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              return(Res)
          })

### probesets
## get host genes from the database
setMethod("probesets", "MirhostDb",
          function(x,
                   columns=listColumns(x, "array_feature"),
                   filter,
                   order.by="probeset_id",
                   order.type="asc",
                   return.type="DataFrame"){
              ## stop if we don't have the required table!!!
              if(!x@have_array_features)
                  stop("The database does not provide the required tables to perform this call!")
              return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              ## attrs <- unique(c(columns, "probeset_id"))
              attrs <- unique(columns)
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              ## silently remove order.by=mat_mirna_id if it's not in attrs...
              if(!any(attrs==order.by))
                  order.by <- ""
              Res <- .getWhat(x=x, columns=attrs, filter=filter,
                              order.by=order.by, order.type=order.type,
                              join="left join", start.table="array_feature")
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              return(Res)
          })

####
## little helper functions: these functions return NULL if the
## required tables are not available, otherwise they return the
## array name
mirfamtab <- function(x){
    if(x@have_mirfam)
        return("mirfam")
    return(NULL)
}
arrayfeature <- function(x){
    if(x@have_array_features){
        return("probeset")
    }else{
        warning("No array_feature table available!")
    }
    return(NULL)
}



### premirnasBy
## get pre_mirna grouped by mat_mirna, host_tx, host_gene, mirfam
setMethod("premirnasBy", "MirhostDb",
          function(x,
                   by="mat_mirna",
                   columns=listColumns(x, "pre_mirna"),
                   filter,
                   return.type="DataFrame",
                   use.names=FALSE
                   ## ,ifEmptyBy=NA
                   ){
              return.type <- match.arg(return.type, c("DataFrame", "data.frame", "GRanges"))
              by <- match.arg(by, c("mat_mirna", "host_tx", "host_gene",
                                    mirfamtab(x), "database", arrayfeature(x)))
              split.by <- paste0(by, "_id")
              if(by=="host_tx"){
                  split.by <- "tx_id"
                  if(use.names){
                      warning("Using transcript IDs instead of transcript names,",
                              " as there are no transcript names in the database!")
                  }
              }
              if(by=="pre_mirna_algn"){
                  by <- "pre_mirna"
              }
              if(by=="host_gene"){
                  split.by <- "gene_id"
                  if(use.names)
                      split.by <- "gene_name"
              }
              if(by=="database"){
                  split.by <- "database"
                  by <- "host_gene"
              }
              if(by=="mat_mirna"){
                  if(use.names)
                      split.by <- "mat_mirna_name"
              }
              if(by=="mirfam"){
                  if(use.names)
                      split.by <- "mirfam_name"
              }
              if(by=="probeset"){
                  by <- "array_feature"
              }
              ## we require here at least the pre_mirna_id and the split.by!!!
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              ## attrs <- unique(c(columns, split.by, "pre_mirna_id"))
              attrs <- unique(c(columns, split.by))
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              order.by=""
              order.type="asc"
              if(return.type=="GRanges"){
                  attrs <- unique(c(attrs, c("seq_name", "seq_strand",
                                             "pre_mirna_seq_start",
                                             "pre_mirna_seq_end")))
              }
              ## elements for which column "by" does not contain a value will not be returned
              ## by the call below
              Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
                              order.type=order.type, join="left join", start.table=by)
              ## check if we have empty ones...
              Empties <- is.na(Res[ , split.by ]) | Res[ , split.by ]==""
              if(any(Empties)){
                  Warnstring <- paste(sum(Empties), split.by, "values are empty!")
                  ## if(!is.na(ifEmptyBy)){
                  ##     Res[ Empties, split.by ] <- ifEmptyBy
                  ##     Warnstring <- paste(Warnstring, "These are returned in the list element named", ifEmptyBy, ".")
                  ## }else{
                  ##     Warnstring <- paste(Warnstring, "These will not be returned.")
                  ## }
                  Warnstring <- paste(Warnstring, "These will not be returned.")
                  warning(Warnstring)
              }
              SplitFac <- factor(Res[ , split.by ])
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              if(return.type=="GRanges"){
                  attrs.metadata <- attrs[ !(attrs %in% c("seq_name", "seq_strand",
                                                          "pre_mirna_seq_start",
                                                          "pre_mirna_seq_end")) ]
                  SI <- seqinfo(x)
                  SI <- SI[ unique(Res$seq_name) ]
                  Res <- GRanges(seqnames=Rle(Res$seq_name),
                                 strand=Rle(Res$seq_strand),
                                 ranges=IRanges(start=Res$pre_mirna_seq_start,
                                                end=Res$pre_mirna_seq_end),
                                 seqinfo=SI,
                                 Res[ , attrs.metadata, drop=FALSE ]
                                 )
              }
              return(split(Res, f=SplitFac))
          })

### matmirnasBy
## get mat_mirna grouped by pre_mirna, host_tx, host_gene, mirfam
setMethod("matmirnasBy", "MirhostDb",
          function(x,
                   by="pre_mirna_algn",
                   columns=listColumns(x, "mat_mirna"),
                   filter,
                   return.type="DataFrame",
                   use.names=FALSE
                   ## ,ifEmptyBy=NA
                   ){
              return.type <- match.arg(return.type,
                                       c("DataFrame", "data.frame", "GRanges"))
              by <- match.arg(by, c("pre_mirna_algn", "pre_mirna", "host_tx",
                                    "host_gene", mirfamtab(x), "database",
                                    arrayfeature(x)))
              split.by <- paste0(by, "_id")
              if(by=="host_tx"){
                  split.by <- "tx_id"
                  if(use.names){
                      warning("Using transcript IDs instead of transcript names,",
                              " as there are no transcript names in the database!")
                  }
              }
              if(by=="pre_mirna_algn"){
                  by <- "pre_mirna"
              }
              if(by=="host_gene"){
                  split.by <- "gene_id"
                  if(use.names){
                      split.by <- "gene_name"
                  }
              }
              if(by=="database"){
                  split.by <- "database"
                  by <- "host_gene"
              }
              if(by=="pre_mirna"){
                  if(use.names)
                      split.by <- "pre_mirna_name"
              }
              if(by=="mirfam"){
                  if(use.names)
                      split.by <- "mirfam_name"
              }
              if(by=="probeset"){
                  by <- "array_feature"
              }
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              attrs <- unique(c(columns, split.by))
              ## we require here at least the pre_mirna_id and the split.by!!!
              ## attrs <- unique(c(columns, split.by, "mat_mirna_id"))
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              order.by <- ""
              order.type <- "asc"
              if(return.type=="GRanges"){
                  attrs <- unique(c(attrs, c("seq_name", "seq_strand", "mat_mirna_seq_start",
                                             "mat_mirna_seq_end")))
              }
              if(any(attrs=="sequence")){
                  ## if we've got sequence as a column, we want in reality to get the sequence of
                  ## the mature miRNA. thus, we need mat_mirna_seq_start, mat_mirna_seq_end,
                  ## pre_mirna_seq_start, pre_mirna_seq_end also.
                  addAttrs <- c("mat_mirna_seq_start", "mat_mirna_seq_end",
                                "pre_mirna_seq_start",
                                "pre_mirna_seq_end", "seq_strand")
                  addAttrs <- addAttrs[!(addAttrs %in% attrs)]
              }else{
                  addAttrs <- character()
              }
              Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
                              order.type=order.type, join="left join", start.table=by)
              ## OK, I need to fix the sequence, if present...
              if(any(attrs=="sequence")){
                  ## what sucks is that the pre-miRNA is NOT the genomic miRNA sequence,
                  ## but the sequence "as.is".
                  PreSeqs <- Res$sequence
                  relStart <- Res$mat_mirna_seq_start - Res$pre_mirna_seq_start
                  relEnd <- Res$mat_mirna_seq_end - Res$pre_mirna_seq_start
                  ## for seq_strand == -1 I have to fix that!
                  if(any(Res$seq_strand < 0)){
                      minusStr <- Res$seq_strand < 0
                      relStart[minusStr] <- Res$pre_mirna_seq_end[minusStr] -
                          Res$mat_mirna_seq_end[minusStr]
                      relEnd[minusStr] <- Res$pre_mirna_seq_end[minusStr] -
                          Res$mat_mirna_seq_start[minusStr]
                  }
                  MatSeqs <- substring(PreSeqs, first=(relStart+1),
                                       last=(relEnd+1))
                  Res$sequence <- MatSeqs
                  Res <- Res[, !(colnames(Res) %in% addAttrs)]
              }

              ## check if we have empty ones...
              Empties <- is.na(Res[ , split.by ]) | Res[ , split.by ]==""
              if(any(Empties)){
                  Warnstring <- paste(sum(Empties), split.by, "values are empty!")
                  ## if(!is.na(ifEmptyBy)){
                  ##     Res[ Empties, split.by ] <- ifEmptyBy
                  ##     Warnstring <- paste(Warnstring, "These are returned in the list element named", ifEmptyBy, ".")
                  ## }else{
                  ##     Warnstring <- paste(Warnstring, "These will not be returned.")
                  ## }
                  Warnstring <- paste(Warnstring, "These will not be returned.")
                  warning(Warnstring)
              }
              SplitFac <- factor(Res[ , split.by ])
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              if(return.type=="GRanges"){
                  attrs.metadata <- attrs[ !(attrs %in% c("seq_name", "seq_strand",
                                                          "mat_mirna_seq_start",
                                                          "mat_mirna_seq_end")) ]
                  SI <- seqinfo(x)
                  SI <- SI[ unique(Res$seq_name) ]
                  Res <- GRanges(seqnames=Rle(Res$seq_name),
                                 strand=Rle(Res$seq_strand),
                                 ranges=IRanges(start=Res$mat_mirna_seq_start,
                                                end=Res$mat_mirna_seq_end),
                                 seqinfo=SI,
                                 Res[ , attrs.metadata, drop=FALSE ]
                                 )
              }
              return(split(Res, f=SplitFac))
          })


### hosttxBy
## it's a little more tricky since we could have empty values for split by.
## get host_tx grouped by host_gene, mat_mirna, pre_mirna, mirfam
setMethod("hosttxBy", "MirhostDb",
          function(x,
                   by="pre_mirna_algn",
                   columns=listColumns(x, "host_tx"),
                   filter,
                   return.type="DataFrame",
                   drop.empty=TRUE,
                   use.names=FALSE
                   ## ,ifEmptyBy=NA
                   ){
              return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
              by <- match.arg(by, c("pre_mirna_algn", "pre_mirna", "host_gene",
                                    "mat_mirna", mirfamtab(x), "database",
                                    arrayfeature(x)))
              split.by <- paste0(by, "_id")
              if(by=="host_tx"){
                  split.by <- "tx_id"
                  if(use.names){
                      warning("Using transcript IDs instead of transcript names,",
                              " as there are no transcript names in the database!")
                  }
              }
              if(by=="pre_mirna_algn"){
                  by <- "pre_mirna"
              }
              if(by=="host_gene"){
                  split.by <- "gene_id"
                  if(use.names){
                      split.by <- "gene_name"
                  }
              }
              if(by=="database"){
                  split.by <- "database"
                  by <- "host_gene"
              }
              if(by=="mirfam"){
                  if(use.names)
                      split.by <- "mirfam_name"
              }
              if(by=="pre_mirna"){
                  if(use.names)
                      split.by <- "pre_mirna_name"
              }
              if(by=="mat_mirna"){
                  if(use.names)
                      split.by <- "mat_mirna_name"
              }
              if(by=="probeset"){
                  by <- "array_feature"
              }
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              attrs <- unique(c(columns, split.by))
              ## attrs <- unique(c(columns, split.by, "tx_id"))
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              order.by <- ""
              order.type <- ""
              Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
                              order.type=order.type, join="left join", start.table=by)
              ## check if we have empty ones...
              Empties <- is.na(Res[ , split.by ]) | Res[ , split.by ]==""
              if(any(Empties)){
                  Warnstring <- paste(sum(Empties), split.by, "values are empty!")
                  ## if(!is.na(ifEmptyBy)){
                  ##     Res[ Empties, split.by ] <- ifEmptyBy
                  ##     Warnstring <- paste(Warnstring, "These are returned in the list element named", ifEmptyBy, ".")
                  ## }else{
                  ##     Warnstring <- paste(Warnstring, "These will not be returned.")
                  ## }
                  Warnstring <- paste(Warnstring, "These will not be returned.")
                  warning(Warnstring)
              }
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              if(drop.empty){
                  if(any(columns == "tx_id")){
                      Res <- Res[ !is.na(Res[ , "tx_id" ]), ]
                  }else{
                      warning("Did not drop empty elements, as required column 'tx_id' was not requested.")
                  }
              }
              return(split(Res, f=factor(Res[ , split.by ])))
          })

### hostgeneBy
## get host_gene grouped by mat_mirna, pre_mirna, mirfam
setMethod("hostgenesBy", "MirhostDb",
          function(x,
                   by="pre_mirna_algn",
                   columns=listColumns(x, "host_gene"),
                   filter,
                   return.type="DataFrame",
                   drop.empty=TRUE,
                   use.names=FALSE
                   ## ,ifEmptyBy=NA
                   ){
              return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
              by <- match.arg(by, c("pre_mirna_algn", "pre_mirna", "host_tx",
                                    "mat_mirna", mirfamtab(x), "database",
                                    arrayfeature(x)))
              split.by <- paste0(by, "_id")
              if(by=="host_tx"){
                  split.by <- "tx_id"
                  if(use.names){
                      warning("Using transcript IDs instead of transcript names,",
                              " as there are no transcript names in the database!")
                  }
              }
              if(by=="pre_mirna_algn"){
                  by <- "pre_mirna"
              }
              if(by=="host_gene"){
                  split.by <- "gene_id"
                  if(use.names){
                      split.by <- "gene_name"
                  }
              }
              if(by=="database"){
                  split.by <- "database"
                  by <- "host_gene"
              }
              if(by=="mirfam"){
                  if(use.names)
                      split.by <- "mirfam_name"
              }
              if(by=="pre_mirna"){
                  if(use.names)
                      split.by <- "pre_mirna_name"
              }
              if(by=="mat_mirna"){
                  if(use.names)
                      split.by <- "mat_mirna_name"
              }
              if(by=="probeset"){
                  by <- "array_feature"
              }
              order.by <- ""
              order.type <- ""
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              attrs <- unique(c(columns, split.by))
              ## attrs <- unique(c(columns, split.by, "gene_id"))
              Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
                              order.type=order.type, join="left join", start.table=by)
              ## Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
              ##                 order.type=order.type, join="left join", start.table="host_gene")
              ## check if we have empty ones...
              Empties <- is.na(Res[ , split.by ]) | Res[ , split.by ]==""
              if(any(Empties)){
                  Warnstring <- paste(sum(Empties), split.by, "values are empty!")
                  ## if(!is.na(ifEmptyBy)){
                  ##     Res[ Empties, split.by ] <- ifEmptyBy
                  ##     Warnstring <- paste(Warnstring, "These are returned in the list element named", ifEmptyBy, ".")
                  ## }else{
                  ##     Warnstring <- paste(Warnstring, "These will not be returned.")
                  ## }
                  Warnstring <- paste(Warnstring, "These will not be returned.")
                  warning(Warnstring)
              }
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              if(drop.empty){
                  if(any(columns == "gene_id")){
                      Res <- Res[ !is.na(Res[ , "gene_id" ]), ]
                  }else{
                      warning("Did not drop empty elements, as required column 'gene_id' was not requested.")
                  }
              }
              return(split(Res, f=factor(Res[ , split.by ])))
          })


### probesetsBy
## get probe sets grouped by something.
setMethod("probesetsBy", "MirhostDb",
          function(x,
                   by="pre_mirna_algn",
                   columns=listColumns(x, "array_feature"),
                   filter,
                   return.type="DataFrame",
                   drop.empty=TRUE,
                   use.names=FALSE
                   ## ,ifEmptyBy=NA
                   ){
              if(!x@have_array_features)
                  stop("The database does not provide the required tables to perform this call!")
              return.type <- match.arg(return.type, c("DataFrame", "data.frame"))
              by <- match.arg(by, c("pre_mirna_algn", "pre_mirna", "host_tx",
                                    "mat_mirna", mirfamtab(x), "database",
                                    "host_gene"))
              split.by <- paste0(by, "_id")
              if(by=="host_tx"){
                  split.by <- "tx_id"
                  if(use.names){
                      warning("Using transcript IDs instead of transcript names,",
                              " as there are no transcript names in the database!")
                  }
              }
              if(by=="pre_mirna_algn"){
                  by <- "pre_mirna"
              }
              if(by=="host_gene"){
                  split.by <- "gene_id"
                  if(use.names){
                      split.by <- "gene_name"
                  }
              }
              if(by=="database"){
                  split.by <- "database"
                  by <- "host_gene"
              }
              if(by=="mirfam"){
                  if(use.names)
                      split.by <- "mirfam_name"
              }
              if(by=="pre_mirna"){
                  if(use.names)
                      split.by <- "pre_mirna_name"
              }
              if(by=="mat_mirna"){
                  if(use.names)
                      split.by <- "mat_mirna_name"
              }
              order.by <- ""
              order.type <- ""
              if(missing(filter))
                  filter <- list()
              filter <- checkIsBasicFilter(filter)
              if(length(columns) == 0)
                  stop("Length of argument 'columns' is 0!")
              attrs <- unique(c(columns, split.by))
              ## attrs <- unique(c(columns, split.by, "probeset_id"))
              Res <- .getWhat(x=x, columns=attrs, filter=filter, order.by=order.by,
                              order.type=order.type, join="left join", start.table=by)
              ## check if we have empty ones...
              Empties <- is.na(Res[ , split.by ]) | Res[ , split.by ]==""
              if(any(Empties)){
                  Warnstring <- paste(sum(Empties), split.by, "values are empty!")
                  ## if(!is.na(ifEmptyBy)){
                  ##     Res[ Empties, split.by ] <- ifEmptyBy
                  ##     Warnstring <- paste(Warnstring, "These are returned in the list element named", ifEmptyBy, ".")
                  ## }else{
                  ##     Warnstring <- paste(Warnstring, "These will not be returned.")
                  ## }
                  Warnstring <- paste(Warnstring, "These will not be returned.")
                  warning(Warnstring)
              }
              if(return.type=="DataFrame"){
                  Res <- DataFrame(Res)
              }
              if(drop.empty){
                  if(any(columns == "probeset_id")){
                      Res <- Res[ !is.na(Res[ , "probeset_id" ]), ]
                  }else{
                      warning("Did not drop empty elements, as required column 'probeset_id' was not requested.")
                  }
              }
              return(split(Res, f=factor(Res[ , split.by ])))
          })


##***********************************************************************
##
##     summary methods
##
##***********************************************************************
setMethod("matmirnasInMultiplePremirnas",
          "MirhostDb",
          function(x,
                   columns=c(listColumns(x, "mat_mirna"), "pre_mirna_id", "pre_mirna_name"),
                   filter=list(),
                   return.type="DataFrame"){
              columns <- unique(c("mat_mirna_name", "mat_mirna_id", columns,
                                  "pre_mirna_name", "pre_mirna_id"))
              if(return.type == "GRanges"){
                  columns <- unique(c(columns, "seq_name", "seq_strand",
                                      "mat_mirna_seq_start", "mat_mirna_seq_end"))
              }
              All <- matmirnas(x, columns=columns, filter=filter,
                               return.type="data.frame")
              All <- split(All, All$mat_mirna_id)
              multis <- unlist(lapply(All, function(z){
                  z <- unique(z[, c("mat_mirna_id", "pre_mirna_id")])
                  nrow(z) > 1
              }))
              All <- do.call(rbind, All[multis])
              if(return.type == "DataFrame")
                  All <- DataFrame(All)
              if(return.type=="GRanges"){
                  meta.attrs <- columns[ !(columns %in% c("seq_name",
                                                          "seq_strand",
                                                          "mat_mirna_seq_start",
                                                          "mat_mirna_seq_end")) ]
                  SI <- seqinfo(x)
                  SI <- SI[unique(All$seq_name)]
                  GR <- GRanges(seqnames=Rle(All$seq_name),
                                ranges=IRanges(start=All$mat_mirna_seq_start,
                                               end=All$mat_mirna_seq_end),
                                strand=Rle(All$seq_strand),
                                seqinfo=SI,
                                All[ , meta.attrs, drop=FALSE ]
                                )
                  return(GR)
              }
              rownames(All) <- NULL
              return(All)
          })

setMethod("premirnasWithMultipleAlignments", "MirhostDb",
          function(x,
                   columns=listColumns(x, "pre_mirna"),
                   filter=list(),
                   return.type="DataFrame"){
              columns <- unique(c(columns, "pre_mirna_algn_id",
                                  "pre_mirna_name", "pre_mirna_id"))
              if(return.type == "GRanges"){
                  columns <- unique(c(columns, "seq_name", "seq_strand",
                                      "pre_mirna_seq_start", "pre_mirna_seq_end"))
              }
              All <- premirnas(x, columns=columns, filter=filter,
                               return.type="data.frame")
              All <- split(All, All$pre_mirna_id)
              multis <- unlist(lapply(All, function(z){
                  z <- unique(z[, c("pre_mirna_id", "pre_mirna_algn_id")])
                  nrow(z) > 1
              }))
              All <- do.call(rbind, All[multis])
              if(return.type == "DataFrame")
                  All <- DataFrame(All)
              if(return.type=="GRanges"){
                  meta.attrs <- columns[ !(columns %in% c("seq_name",
                                                          "seq_strand",
                                                          "pre_mirna_seq_start",
                                                          "pre_mirna_seq_end")) ]
                  SI <- seqinfo(x)
                  SI <- SI[unique(All$seq_name)]
                  GR <- GRanges(seqnames=Rle(All$seq_name),
                                ranges=IRanges(start=All$pre_mirna_seq_start,
                                               end=All$pre_mirna_seq_end),
                                strand=Rle(All$seq_strand),
                                seqinfo=SI,
                                All[ , meta.attrs, drop=FALSE ]
                                )
                  return(GR)
              }
              rownames(All) <- NULL
              return(All)
          })


##***********************************************************************
##
##     methods to match values between host genes and mature mirnas
##
##***********************************************************************
## matchValues
## 1) define the mapping between the names of x and names of y.
## 2) call doMatchValues to match the values.
setMethod("pairData", signature(x="numeric", y="numeric", object="MirhostDb",
                                chooseFunX="function", chooseFunY="function",
                                xNamesAre="character", yNamesAre="character"),
          function(x, y, object, chooseFunX=chooseAll, chooseFunY=chooseAll,
                   xNamesAre="mat_mirna_name", yNamesAre="probeset_id"){
              ## first, check xNamesAre, yNamesAre
              suppNames <- c("pre_mirna_id", "pre_mirna_name", "tx_id",
                             "mat_mirna_id", "mat_mirna_name",
                             "gene_id", "gene_name", "probeset_id")
              ## define the table name for each...
              names(suppNames) <- c("pre_mirna", "pre_mirna", "host_tx",
                                    "mat_mirna", "mat_mirna",
                                    "host_gene", "host_gene", "array_feature")
              xNamesAre <- match.arg(xNamesAre, suppNames)
              yNamesAre <- match.arg(yNamesAre, suppNames)
              ## just check if we do have probesets available... if needed...
              if(xNamesAre == "probeset_id" | yNamesAre == "probeset_id"){
                  if(!object@have_array_features)
                      stop("No probe set definitions available in the database!",
                           " Thus 'probeset_id' can not be used for arguments 'xNamesAre' or 'yNamesAre'.")
              }
              ## use .getWhat (without filter), but need to define start.table...
              allData <- .getWhat(object, columns=c(xNamesAre, yNamesAre),
                                  start.table=names(suppNames)[suppNames == xNamesAre])
          })


## setMethod("transferValues", signature(x="numeric", object="MirhostDb", xNamesAre="character",
##                                       toNames="character", solveFun="missing",
##                                       filter="missing", na.rm="missing"),
##           function(x, object, xNamesAre="mat_mirna_name", toNames="pre_mirna_name"){
##               transferValues(x=x, object=object, xNamesAre=xNamesAre, toNames=toNames,
##                              solveFun=chooseAll,
##                              list(DatabaseFilter("core"),
##                                GenebiotypeFilter("miRNA", condition="!="),
##                                ArrayFilter("HG-U133_Plus_2")),
##                    na.rm=FALSE)
##           })

## function to:
## transfer values, map values.
## input is some values for specific ids, either mature_mirnas, pre_mirnas, probeset_ids,
## gene_ids, transcript_ids.
## idea is:
## 1) map the original ids to the result ids.
## 2) transfer the values.
## sounds simple, but isn't:
## a) one mature miRNA can be encoded in more than one pre-miRNA: assign one value to
##    potentially more than one pre-miRNA
## b) each pre-miRNA however encodes two mature miRNAs; so, map two values to one pre-miRNA.
## c) one host transcript can encode 1 to x pre-miRNAs, thus, map x pre-miRNAs to one host
##    transcriot.
## d) each host genes/transcript can be detected by more than one probe set.
##
## filter: don't change filter.
## setMethod("transferValues", signature(x="numeric", object="MirhostDb", xNamesAre="character",
##                                       toNames="character", solveFun="function", filter="ANY",
##                                       na.rm="logical"),
setMethod("transferValues", signature(x="numeric", object="MirhostDb"),
          function(x, object, xNamesAre="mat_mirna_name", toNames="pre_mirna_name",
                   solveFun=chooseAll,
                   filter=list(),
                   na.rm=FALSE){
              if(is.null(names(x)))
                  stop("'x' has to be a named numeric vector!")
              if(class(filter)!="list")
                  filter <- list(filter)
              ## 1) check xNamesAre, based on that the toNames.
              suppNames <- c("pre_mirna_id", "pre_mirna_name", "tx_id", "mat_mirna_id",
                             "mat_mirna_name",
                             "gene_id", "gene_name", "probeset_id")
              ## define the table name for each...
              names(suppNames) <- c("pre_mirna", "pre_mirna", "host_tx", "mat_mirna",
                                    "mat_mirna",
                                    "host_gene", "host_gene", "array_feature")
              xNamesAre <- match.arg(xNamesAre, suppNames)
              ## OK, now I've got xNames, let's see to what we want to transfer that.
              toNames <- match.arg(toNames, suppNames)
              if(toNames == xNamesAre)
                  stop("Value for 'xNamesAre' should be different from 'toNames'!")
              ## just check if we do have probesets available... if needed...
              if(xNamesAre == "probeset_id" | toNames == "probeset_id"){
                  if(!object@have_array_features)
                      stop("No probe set definitions available in the database!",
                           " Thus 'probeset_id' can not be used for arguments",
                           " 'xNamesAre' or 'yNamesAre'.")
              }
              ## check, if we're staying within the miRNA world we don't need the filters:
              ## if(length(grep(xNamesAre, pattern="mirna")) > 0 & length(grep(toNames, pattern="mirna"))){
              ##     filter <- list()
              ## }
              if(toNames != "probeset_id"){
                  ## get rid of the array filter!
                  af <- unlist(lapply(filter, function(f){
                      return(class(f)=="ArrayFilter")
                  }))
                  if(any(af))
                      filter <- filter[!af]
              }
              ## 2) map between ids
              allData <- .getWhat(object, columns=c(xNamesAre, toNames),
                                  start.table=names(suppNames)[suppNames == xNamesAre],
                                  filter=filter)
              if(nrow(allData) == 0){
                  ## something bad happened.
                  stop("Could not get any mapping data! Please re-evaluate",
                       " eventually specified filters!")
              }
              ## 3) subset the input values.
              xNoMapIdx <- which(!(names(x) %in% allData[, xNamesAre]))
              if(length(xNoMapIdx) > 0){
                  warning(paste0("Got no mapping for ids: ",
                                 paste(names(x)[xNoMapIdx], collapse=", ")))
                  xNoMap <- x[xNoMapIdx]
                  x <- x[-xNoMapIdx]
                  if(length(x) == 0)
                      stop("Could not map any of the names of 'x' to ", toNames, "!")
              }
              ## 4) assign values, solve multi-mapping using solveFun.
              allData <- allData[allData[, xNamesAre] %in% names(x), ]
              xOrig <- x    ## keep the original input
              x <- x[allData[, xNamesAre]]
              names(x) <- allData[, toNames]
              ## OK, now choose the value,
              x <- doSelectData(x, chooseFunX=solveFun)
              ## get the original name matching the toNames and the value:
              xvals <- cbind(x, origName=allData[x$x.idx, xNamesAre])
              xvals <- xvals[, c("x", "name", "origName")]
              ## NA values???
              if(!na.rm & length(xNoMapIdx) > 0){
                  empties <- data.frame(x=xNoMap, name=rep(NA, length(xNoMap)),
                                        origName=names(xNoMap),
                                        stringsAsFactors=FALSE)
                  xvals <- rbind(xvals, empties)
              }
              names(xvals) <- c("x", toNames, xNamesAre)
              rownames(xvals) <- NULL
              return(xvals)
          })



##***********************************************************************
##
##     private methods
##***********************************************************************
### tablesForColumns
## returns the tables for the specified columns.
setMethod("tablesForColumns", "MirhostDb", function(x, columns, ...){
    if(missing(columns))
        stop("No columns submitted!")
    bm <- unlist(lapply(listTables(x), function(z){
        return(any(z %in% columns))
    }))
    if(!any(bm))
        return(NULL)
    Tables <- names(bm)[ bm ]
    Tables <- Tables[ !(Tables %in% c("chromosome", "metadata")) ]
    return(Tables)
})

### cleanColumns
## checks the provided columns and removes all that are not in the database.
setMethod("cleanColumns", "MirhostDb", function(x, columns){
    if(missing(columns))
        stop("No columns provided!")
    TableExcludes <- "metadata"
    Tables <- listTables(x)
    Allowed <- unlist(Tables[ !(names(Tables) %in% TableExcludes) ],
                      use.names=FALSE)
    bm <- columns %in% Allowed
    removed <- columns[ !bm ]
    if(length(removed) > 0){
        warning("The columns ", paste(removed, collapse=","),
                " have been removed as they do not match any column names in the database.")
    }
    columns <- columns[ bm ]
    if(length(columns)==0)
        stop("Columns are not known to the database!")
    return(columns)
})
## returns the table names ordered by degree, i.e. edges to other tables
setMethod("tablesByDegree", "MirhostDb",
          function(x, tab=names(listTables(x)), ...){
              ## ## to do this with a graph:
              ## DBgraph <- graphNEL(nodes=c("pre_mirna", "host_tx", "mat_mirna",
              ##                     "host_gene", "mirfam", "pre_mirna_seq", "information"),
              ##                  edgeL=list(pre_mirna=c("host_tx", "mat_mirna", "mirfam", "pre_mirna_seq"),
              ##                      host_tx=c("pre_mirna", "host_gene"),
              ##                      mat_mirna="pre_mirna",
              ##                      host_gene="host_tx",
              ##                      mirfam="pre_mirna",
              ##                      pre_mirna_seq="pre_mirna"
              ##                          ))
              ## Tab <- names(sort(degree(DBgraph), decreasing=TRUE))
              ##Table.order <- c(pre_mirna=1, mat_mirna=2, host_tx=3, host_gene=4, mirfam=5, pre_mirna_sequence=6)
              Table.order <- c(pre_mirna=1, host_tx=2, mat_mirna=3,
                               host_gene=4, mirfam=5, pre_mirna_sequence=6,
                               array_feature=7, metadata=8,
                               chromosome=9)
              Tab <- tab[ order(Table.order[ tab ]) ]
              return(Tab)
          })




