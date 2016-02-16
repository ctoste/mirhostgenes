## simple function that loads the database and creates the
## MirhostDb object.
MirhostDb <- function(x){
    options(useFancyQuotes=FALSE)
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname=x, flags=SQLITE_RO)
    tables <- dbListTables(con)
    ## read the attributes for these tables.
    Tables <- vector(length=length(tables), "list")
    for(i in 1:length(Tables)){
        Tables[[ i ]] <- colnames(dbGetQuery(con,
                                             paste0("select * from ",
                                                    tables[ i ], " limit 1")))
    }
    names(Tables) <- tables
    have_mirfam <- FALSE
    have_pre_sequence <- FALSE
    have_array_features <- FALSE
    if(any(tables=="mirfam"))
        have_mirfam <- TRUE
    if(any(tables=="pre_mirna_sequence"))
        have_pre_sequence <- TRUE
    if(any(tables=="array_feature"))
        have_array_features <- TRUE
    ## check if we've got some additional stuff:
    havePreConf <- FALSE
    havePreRead <- FALSE
    haveMatConf <- FALSE
    haveMatRead <- FALSE
    PreConfMeta <- .getMetaDataValue(con, "have_premirna_confidence")
    if(!is.na(PreConfMeta))
        havePreConf <- ifelse(PreConfMeta == "TRUE", yes=TRUE, no=FALSE)
    PreReadMeta <- .getMetaDataValue(con, "have_premirna_readcount")
    if(!is.na(PreReadMeta))
        havePreRead <- ifelse(PreReadMeta == "TRUE", yes=TRUE, no=FALSE)
    MatConfMeta <- .getMetaDataValue(con, "have_matmirna_confidence")
    if(!is.na(MatConfMeta))
        haveMatConf <- ifelse(MatConfMeta == "TRUE", yes=TRUE, no=FALSE)
    MatReadMeta <- .getMetaDataValue(con, "have_matmirna_readcount")
    if(!is.na(MatReadMeta))
        haveMatRead <- ifelse(MatReadMeta == "TRUE", yes=TRUE, no=FALSE)
    MDB <- new("MirhostDb", con=con, tables=Tables,
               have_mirfam=have_mirfam,
               have_pre_sequence=have_pre_sequence,
               have_array_features=have_array_features,
               have_premirna_confidence=havePreConf,
               have_premirna_readcount=havePreRead,
               have_matmirna_confidence=haveMatConf,
               have_matmirna_readcount=haveMatRead
               )
    return(MDB)
}

## x is the connection to the database, name is the name of the entry to fetch
.getMetaDataValue <- function(x, name){
    return(dbGetQuery(x, paste0("select value from metadata where name='",
                                name, "'"))[ 1, 1])
}



## filter is a list of Filters; if there is any filter that has a slotName match.case and
## the value of that slot is FALSE, then the function will return TRUE.
.doNocase <- function(filter){
    if(length(filter)==0)
        return(FALSE)
    Res <- unlist(lapply(filter, function(x){
        if(any(slotNames(x)=="match.case")){
            return(x@match.case)
        }else{
            return(TRUE)
        }
    }))
    return(any(Res)==FALSE)
}

####
## Note: that's the central function that checks which tables are needed for the
## least expensive join!!! The names of the tables should then also be submitted
## to any other method that calls prefixColumns (e.g. where of the Filter classes)
##
## this function checks:
## a) for multi-table attributes, selects the table with the highest degree (i.e.
##    the table connected to most other tables).
## b) pre-pend (inverse of append ;)) the table name to the attribute name.
## returns a list, names being the tables and the values being the attributes
## named: <table name>.<attribute name>
## clean: whether a cleanColumns should be called on the submitted attributes.
## with.tables: force the prefix to be specifically on the submitted tables.
prefixColumns <- function(x, columns, clean=TRUE, with.tables){
    if(missing(columns))
        stop("columns is empty! No columns provided!")
    ## first get to the tables that contain these columns
    Tab <- listTables(x)   ## returns the tables by degree!
    if(!missing(with.tables)){
        with.tables <- with.tables[ with.tables %in% names(Tab) ]
        if(length(with.tables) > 0)
            Tab <- Tab[ with.tables ]
        if(length(Tab)==0)
            stop("None of the tables submitted with with.tables is present in the database!")
    }
    if(clean)
        columns <- cleanColumns(x, columns)
    if(length(columns)==0){
        return(NULL)
    }
    ## group the columns by table.
    columns.bytable <- sapply(Tab, function(z){
        return(z[ z %in% columns ])
    }, simplify=FALSE, USE.NAMES=TRUE)
    ## kick out empty tables...
    columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable, function(z){
        return(length(z) > 0)
    })) ]
    if(length(columns.bytable)==0)
        stop("No columns available!")
    have.columns <- NULL
    ## new approach! order the tables by number of elements, and after that, re-order them.
    columns.bytable <- columns.bytable[ order(unlist(lapply(columns.bytable, length)),
                                              decreasing=TRUE) ]
    ## has to be a for loop!!!
    ## loop throught the columns by table and sequentially kick out columns
    ## for the current table if they where already
    ## in a previous (more relevant) table
    for(i in 1:length(columns.bytable)){
        bm <- columns.bytable[[ i ]] %in% have.columns
        keepvals <- columns.bytable[[ i  ]][ !bm ]   ## keep those
        if(length(keepvals) > 0){
            have.columns <- c(have.columns, keepvals)
        }
        if(length(keepvals) > 0){
            columns.bytable[[ i ]] <- paste(names(columns.bytable)[ i ], keepvals, sep=".")
        }else{
            columns.bytable[[ i ]] <- keepvals
        }
    }
    ## kick out those tables with no elements left...
    columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable, function(z){
        return(length(z) > 0)
    })) ]
    return(columns.bytable)
}


## define a function to create a join query based on columns
## this function has to first get all tables that contain the columns,
## and then select, for columns present in more than one
## table (i.e. pre_mirna_id, gene_id) select the table with the higher
## degree (since that one is better to join).
## x... MirhostDb
## columns... the columns
joinQueryOnColumns <- function(x, columns, join="join", start.table){
    columns.bytable <- prefixColumns(x, columns)
    ## based on that we can build the query based on the tables we've got.
    ## Note that the function internally
    ## adds tables that might be needed for the join.
    Query <- joinQueryOnTables(x, names(columns.bytable), join=join, start.table=start.table)
    return(Query)
}


## only list direct joins!!!
.JOINS <- rbind(
    c("pre_mirna", "mat_mirna", "on (pre_mirna.pre_mirna_algn_id=mat_mirna.pre_mirna_algn_id)"),
    c("pre_mirna", "host_tx", "on (pre_mirna.pre_mirna_algn_id=host_tx.pre_mirna_algn_id)"),
    c("pre_mirna", "mirfam", "on (pre_mirna.pre_mirna_id=mirfam.pre_mirna_id)"),
    c("pre_mirna", "pre_mirna_sequence",
      "on (pre_mirna.pre_mirna_id=pre_mirna_sequence.pre_mirna_id)"),
    c("host_tx", "host_gene", "on (host_tx.gene_id=host_gene.gene_id)"),
    c("host_tx", "array_feature", "on (host_tx.tx_id=array_feature.tx_id)")
    )


###  JOINS:
##   pre_mirna - mat_mirna ON pre_mirna_pk=pre_mirna_fk
##   pre_mirna - host_tx ON pre_mirna_pk=pre_mirna_fk
##   host_tx - host_gene ON gene_id=gene_id
##   pre_mirna_sequence - pre_mirna ON pre_mirna_id=pre_mirna_id
##   mirfam - pre_mirna ON pre_mirna_id=pre_mirna_id
## x... MirhostDb object.
## tab... character vector of table names to be joined.
## join... how the tables should be joined. allowed are join and left join.
## start.table... optionally specify the "most left" table, i.e. the start table from
## which the values should be fetched first.
joinQueryOnTables <- function(x, tab, join="join", start.table){
    ## just to be on the save side: evaluate whether we have all required tables to join;
    ## this will also ensure that the order is by degree.
    tab <- addRequiredTables(x, tab)
    if(!missing(start.table)){
        ## check if we do have this table...
        if(any(tab==start.table)){
            ## re-order the tables such that the start table is the first...
            tab <- c(start.table, tab[ tab!=start.table ])
        }
    }
    Query <- tab[ 1 ]
    previous.table <- tab[ 1 ]
    remaining.table <- tab[ -1 ]
    ## now, have to check what joins I can have starting from this first one and based
    ## on all tables I have.
    while(length(previous.table)!=length(tab)){
        ## repeat until I don't have joined all tables!
        ## check which joins are available for the last table(s)
        previous.idx <- which(apply(.JOINS[ , 1:2, drop=FALSE ], MARGIN=1,
                                    function(z){ any(z %in% previous.table) }))
        if(length(previous.idx)==0)
            stop("Argh, something weird happened...")
        ## now check if in the remaining tables there is one that could be joined to this one.
        tmp <- .JOINS[ previous.idx, , drop=FALSE ]
        remaining.idx <- which(apply(tmp[ , 1:2, drop=FALSE ], MARGIN=1,
                                     function(z){ any(z %in% remaining.table) }))
        ## add this table to the previous.table vector
        previous.table <- unique(c(previous.table,
                                   tmp[ remaining.idx[ 1 ], 1:2 ][ tmp[ remaining.idx[ 1 ], 1:2 ]!=previous.table[ length(previous.table) ] ]))
        remaining.table <- remaining.table[ remaining.table!=previous.table[ length(previous.table) ] ]
        Query <- paste(Query, join, previous.table[ length(previous.table) ],
                       tmp[ remaining.idx[ 1 ], 3 ])
    }
    return(Query)
}


###
## add additional tables in case the submitted tables are not directly connected
## and can thus not be joined.
addRequiredTables <- function(x, tab){
    ## dash it, as long as I can't find a way to get connected objects in a
    ## graph I'll do it manually...
    ## so, if I have host_gene in it, and any other table, I need host_tx
    if(any(tab=="host_gene") & length(tab) > 1){
        tab <- unique(c(tab, "host_tx"))
    }
    if(any(tab=="array_feature") & length(tab) > 1){
        tab <- unique(c(tab, "host_tx"))
    }
    ## if we have host_tx and any table other than host_gene, we need pre_mirna
    if(any(tab=="host_tx") &
       any(tab %in% c("mat_mirna", "mirfam", "pre_mirna_sequence"))){
        tab <- unique(c(tab, "pre_mirna"))
    }
    ## if we have mat_mirna, pre_mirna_sequence or mir_fam and any other table, we need pre_mirna
    if(any(tab %in% c("mat_mirna", "mirfam", "pre_mirna_sequence")) & length(tab) > 1){
        tab <- unique(c(tab, "pre_mirna"))
    }
    return(tablesByDegree(x, tab))
}


.buildQuery <- function(x, columns, filter=list(), order.by="", order.type="asc",
                        group.by, skip.order.check=FALSE, join="join", start.table, return.all.columns=TRUE){
    resultcolumns <- columns    ## just to remember what we really want to give back
    join <- match.arg(join, c("join", "left join", "left outer join"))
    ## 1) get all column names from the filters also removing the prefix.
    if(class(filter)!="list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    if(length(filter) > 0){
        ## check filter!
        ## add the columns needed for the filter
        filtercolumns <- unlist(lapply(filter, column, x))
        ## remove the prefix (column name for these)
        filtercolumns <- sapply(filtercolumns, removePrefix, USE.NAMES=FALSE)
        columns <- unique(c(columns, filtercolumns))
    }
    ## 2) get all column names for the order.by:
    if(order.by!=""){
        ## if we have skip.order.check set we use the order.by as is.
        if(!skip.order.check){
            order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
            order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
            ## allow only order.by that are also in the columns.
            order.by.nocolumns <- order.by[ !(order.by %in% columns) ]
            order.by <- order.by[ order.by %in% columns ]
            if(length(order.by.nocolumns) > 0){
                warning("columns provided in order.by (",
                        paste(order.by.nocolumns, collapse=","),
                        ") are not in columns and were thus removed.")
            }
            if(length(order.by)==0){
                order.by <- ""
            }
        }
    }else{
        order.by <- ""
    }
    ## Note: order by is now a vector!!!
    ## columns are now all columns that we want to fetch or that we need to
    ## filter or to sort.
    ## 3) check which tables we need for all of these columns:
    need.tables <- names(prefixColumns(x, columns))
    ##
    ## Now we can begin to build the query parts!
    ## a) the query part that joins all required tables.
    joinquery <- joinQueryOnColumns(x, columns=columns, join=join,
                                    start.table=start.table)
    ## b) the filter part of the query
    if(length(filter) > 0){
        filterquery <- paste(" where",
                             paste(unlist(lapply(filter, where, x,
                                                 with.tables=need.tables)),
                                   collapse=" and "))
    }else{
        filterquery <- ""
    }
    ## c) the order part of the query
    if(order.by!=""){
        if(!skip.order.check){
            order.by <- paste(unlist(prefixColumns(x=x, columns=order.by,
                                                   with.tables=need.tables),
                                     use.names=FALSE), collapse=",")
        }
        orderquery <- paste(" order by", order.by, order.type)
    }else{
        orderquery <- ""
    }
    ## And finally build the final query
    if(return.all.columns){
        resultcolumns <- columns
    }
    finalquery <- paste0("select distinct ",
                         paste(unlist(prefixColumns(x,
                                                    resultcolumns,
                                                    with.tables=need.tables),
                                      use.names=FALSE), collapse=","),
                         " from ",
                         joinquery,
                         filterquery,
                         orderquery
                         )
    return(finalquery)
}


## the buildQuery function; basically create the SQL query.
## x... MirhostDb
## columns... the columns to retrieve.
## filter... list of filters
## join: either "join" or "left join"
## start.table: specify from which table we want to start the join; useful for "left join"
## WARNING!!!! DOES NOT WORK!!!
.buildQueryOld <- function(x, columns, filter, order.by="", order.type="asc",
                           group.by, skip.order.check=FALSE, join="join", start.table, return.all.columns=TRUE){
    resultcolumns <- columns    ## just to remember what we really want to give back
    join <- match.arg(join, c("join", "left join", "left outer join"))
    ## first checking the filters and eventually add required columns to the columns:
    if(missing(filter)){
        filter <- list()
    }
    if(class(filter)!="list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    if(length(filter) > 0){
        ## check filter!
        ## add the columns needed for the filter
        filtercolumns <- unlist(lapply(filter, column, x))   ## this returns the full named columns!
        filtercolumns <- sapply(filtercolumns, removePrefix, USE.NAMES=FALSE)
        columns <- unique(c(columns, filtercolumns))
        ## next we're building the where query.
        wherequery <- paste(" where",
                            paste(unlist(lapply(filter, where, x)), collapse=" and "))
    }else{
        wherequery <- ""
    }
    ## should we do an order.by?
    if(!missing(order.by) & order.by!=""){
        ## if we have skip.order.check set we use the order.by as is.
        if(!skip.order.check){
            order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
            order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
            ## allow only order.by that are also in the columns.
            order.by.nocolumns <- order.by[ !(order.by %in% columns) ]
            order.by <- order.by[ order.by %in% columns ]
            if(length(order.by.nocolumns) > 0){
                warning("columns provided in order.by (",
                        paste(order.by.nocolumns, collapse=","),
                        ") are not in columns and were thus removed.")
            }
            if(length(order.by)==0){
                order.by <- ""
            }else{
                ## have to pre-pend the database table name...
                order.by <- paste(unlist(prefixColumns(x=x, columns=order.by),
                                         use.names=FALSE), collapse=",")
            }
        }
    }else{
        order.by <- ""
    }
    ## OK, build that query.
    if(order.by!=""){
        orderquery <- paste(" order by", order.by, order.type)
    }else{
        orderquery <- ""
    }
    ## now build the join query that joins all required tables.
    joinquery <- joinQueryOnColumns(x, columns=columns, join=join, start.table=start.table)
    if(return.all.columns){
        resultcolumns <- columns
    }
    finalquery <- paste0("select distinct ",
                         paste(unlist(prefixColumns(x, resultcolumns),
                                      use.names=FALSE), collapse=","),
                         " from ",
                         joinquery,
                         wherequery,
                         orderquery
                         )
    return(finalquery)
}

## just to add another layer; basically just calls buildQuery and executes the query
## return.all.columns: returns also columns added because of the filters.
.getWhat <- function(x, columns, filter=list(), order.by="", order.type="asc",
                     group.by, skip.order.check=FALSE, join="join", start.table, return.all.columns=TRUE){
    Q <- .buildQuery(x=x, columns=columns, filter=filter, order.by=order.by,
                     order.type=order.type, group.by=group.by,
                     skip.order.check=skip.order.check, join=join,
                     start.table=start.table, return.all.columns=return.all.columns)
    return(dbGetQuery(x@con, Q))
}

## remove the prefix again...
removePrefix <- function(x, split=".", fixed=TRUE){
    return(sapply(x, function(z){
        tmp <- unlist(strsplit(z, split=split, fixed=fixed))
        return(tmp[ length(tmp) ])
    }))
}


####
## Note: functions we might want to import from ensembldb instead:
## removePrefix
## .buildQuery
## joinQueryOnTables
## joinQueryOnColumns
## prefixColumns

