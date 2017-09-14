tableNames <- function(x){
    return(names(listTables(x)))
}
##***********************************************************************
##
##     Implementations for list
##
##***********************************************************************
setMethod("where", signature(object="list",db="MirhostDb"),
          function(object, db, ...){
              wherequery <- paste(" where", paste(unlist(lapply(object, where, db)),
                                                  collapse=" and "))
              return(wherequery)
          })
setMethod("value", "BasicFilter", function(x) {
    x@value
})
setReplaceMethod("value", "BasicFilter", function(x, value) {
    x@value <- as.character(value)
    x
})
setReplaceMethod("condition", "BasicFilter", function(x, value) {
    x@condition <- as.character(value)
    x
})



##***********************************************************************
##
##     Implementations for PositionFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="PositionFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PositionFilter, db=missing,",
                      "with.tables=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="PositionFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PositionFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              ## ok, in database MirhostDb, the columns for position could be
              ## in_intron, or in_exon.
              Val <- object@value
              if(!(Val %in% c("exonic", "intronic", "both")))
                  stop("A value of '", Val, "' is not supported for MirhostDb!",
                       " Only 'exonic', 'intronic' or 'both' are allowed.")
              if(Val == "exonic")
                  return("host_tx.in_exon")
              if(Val == "intronic")
                  return("host_tx.in_intron")
              if(Val == "both")
                  return(c("host_tx.in_exon", "host_tx.in_intron"))
              return(NA)
          })
## that's the version without a database specified.
setMethod("where", signature(object="PositionFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          }
          )
setMethod("where", signature(object="PositionFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PositionFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              ## note: value can be exonic, intronic or both. We're going to
              ## translate that string into the corresponding where clause that
              ## checks if in_intron or in_exon (or both) are > 0.
              Str <- column(object, db, with.tables=with.tables)
              Str <- sapply(Str, function(x){
                  paste(x, "> 0")
              })
              return(paste(Str, collapse=" and "))
          })


##***********************************************************************
##
##     Implementations for AlignmentidFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="AlignmentidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=AlignmentidFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="AlignmentidFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="AlignmentidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_algn_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="AlignmentidFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="AlignmentidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="AlignmentidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for ArrayFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="ArrayFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=ArrayFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="ArrayFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ArrayFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  ## means we don't have that data available...
                  warning("Can not use a ArrayFilter, as the database doesn't",
                          " have the required table.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "array_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="ArrayFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="ArrayFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ArrayFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  ## means we don't have that data available...
                  warning("Can not use a ArrayFilter, as the database doesn't",
                          " have the required table.\n")
                  return(NULL)
              }
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })



##***********************************************************************
##
##     Implementations for DatabaseFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="DatabaseFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=DatabaseFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="DatabaseFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="DatabaseFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "database", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="DatabaseFilter", db="missing", with.tables="missing"),
          function(object, db, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="DatabaseFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="DatabaseFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })



##***********************************************************************
##
##     Implementations for ProbesetidFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="ProbesetidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=ProbesetidFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="ProbesetidFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ProbesetidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  warning("Can not use a ProbesetidFilter, as the database",
                          " doesn't have the required table.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "probeset_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="ProbesetidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="ProbesetidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ProbesetidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  ## means we don't have that data available...
                  warning("Can not use a ProbesetidFilter, as the database",
                          " doesn't have the required table.\n")
                  return(NULL)
              }
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for MatmirnaConfidence
##
##
##***********************************************************************
setMethod("column", signature(object="MatmirnaConfidence", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatmirnaConfidence, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="MatmirnaConfidence", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MatmirnaConfidence", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_matmirna_confidence){
                  warning("Can not use MatmirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "mat_mirna_confidence", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatmirnaConfidence", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatmirnaConfidence", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatmirnaConfidence", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_matmirna_confidence){
                  ## means we don't have that data available...
                  warning("Can not use MatmirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              ## check manually what values we've got. Thus far we support only
              ## confidence of value "high"
              if(object@value!="high")
                  stop("MatmirnaConfidence can only take the value 'high' for",
                       " db being a MirhostDb database.")
              if(object@value == "high"){
                  suff <- paste0(object@condition, 1)
              }
              ## suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for PremirnaConfidence
##
##
##***********************************************************************
setMethod("column", signature(object="PremirnaConfidence", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PremirnaConfidence,",
                      " db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="PremirnaConfidence", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PremirnaConfidence", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_premirna_confidence){
                  warning("Can not use PremirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "pre_mirna_confidence", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PremirnaConfidence", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PremirnaConfidence", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PremirnaConfidence", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_premirna_confidence){
                  ## means we don't have that data available...
                  warning("Can not use PremirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              ## check manually what values we've got. Thus far we support only
              ## confidence of value "high"
              if(object@value!="high")
                  stop("PremirnaConfidence can only take the value 'high' for",
                       " db being a MirhostDb database.")
              if(object@value == "high"){
                  suff <- paste0(object@condition, 1)
              }
              ## suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for ReadCountFilter
##
##
##***********************************************************************
setMethod("column", signature(object="ReadCountFilter", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=ReadCountFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="ReadCountFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ReadCountFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!(object@of %in% c("mat_mirna", "pre_mirna")))
                  stop("Parameter 'of' of ReadCountFilter for MirhostDb has",
                       " to be either 'mat_mirna' or 'pre_mirna'!")
              if(object@of == "mat_mirna"){
                  if(!db@have_matmirna_readcount){
                      warning("Can not use ReadCountFilter, as the database",
                              " doesn't have the required mat_mirna read count information.")
                      return(NULL)
                  }
              }
              if(object@of == "pre_mirna"){
                  if(!db@have_premirna_readcount){
                      warning("Can not use ReadCountFilter, as the database",
                              " doesn't have the required pre_mirna read count information.")
                      return(NULL)
                  }
              }
              return(unlist(prefixColumns(db, paste0(object@of, "_read_count"),
                                          with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="ReadCountFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="ReadCountFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ReadCountFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!(object@of %in% c("mat_mirna", "pre_mirna")))
                  stop("Parameter 'of' of ReadCountFilter for MirhostDb has to be",
                       " either 'mat_mirna' or 'pre_mirna'!")
              if(object@of == "mat_mirna"){
                  if(!db@have_matmirna_readcount){
                      warning("Can not use ReadCountFilter, as the database doesn't",
                              " have the required mat_mirna read count information.")
                      return(NULL)
                  }
              }
              if(object@of == "pre_mirna"){
                  if(!db@have_premirna_readcount){
                      warning("Can not use ReadCountFilter, as the database doesn't",
                              " have the required pre_mirna read count information.")
                      return(NULL)
                  }
              }
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setReplaceMethod("value", "ReadCountFilter", function(x, value){
    ## Complain if value is not numeric.
    if(!is.numeric(value)){
        suppressWarnings(
            numVal <- as.numeric(value)
        )
        if(is.na(numVal))
            stop("ReadCountFilter only supports numeric values!")
    }
    x@value <- as.character(value)
    return(x)
})

##***********************************************************************
##
##     Implementations for PremirnaFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="PremirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PremirnaFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="PremirnaFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PremirnaFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_name", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PremirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PremirnaFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PremirnaFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object) ## that returns the result from the BasicFilter...
              nc <- NULL
              if(!object@match.case)
                  nc <- "collate nocase"
              return(paste(c(column(object, db, with.tables=with.tables), suff, nc),
                           collapse=" "))
          })


##***********************************************************************
##
##     Implementations for PremirnaidFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="PremirnaidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PremirnaidFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="PremirnaidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PremirnaidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PremirnaidFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PremirnaidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PremirnaidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object)
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for MatmirnaFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MatmirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatmirnaFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MatmirnaFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MatmirnaFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mat_mirna_name", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatmirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatmirnaFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatmirnaFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object)
              nc <- NULL
              if(!object@match.case)
                  nc <- "collate nocase"
              return(paste(c(column(object, db, with.tables=with.tables), suff, nc),
                           collapse=" "))
          })



##***********************************************************************
##
##     Implementations for MatmirnaidFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MatmirnaidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatmirnaidFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MatmirnaidFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, tn))
          })
setMethod("column", signature(object="MatmirnaidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mat_mirna_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatmirnaidFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatmirnaidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatmirnaidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object)
              return(paste(column(object, db, with.tables=with.tables), suff))
          })



##***********************************************************************
##
##     Implementations for MirfamFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MirfamFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MirfamFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MirfamFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MirfamFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_mirfam){
                  ## means we don't have that data available...
                  warning("Can not use a MirfamFilter, as the database doesn't",
                          " have the required table.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "mirfam_name", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MirfamFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MirfamFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MirfamFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_mirfam){
                  ## means we don't have that data available...
                  warning("Can not use a MirfamFilter, as the database doesn't",
                          " have the required table.\n")
                  return(NULL)
              }
              suff <- where(object)
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for MirfamidFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MirfamidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MirfamidFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MirfamidFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MirfamidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mirfam_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MirfamidFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MirfamidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MirfamidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object)
              return(paste(column(object, db, with.tables=with.tables), suff))
          })




##***********************************************************************
##
##     Classes imported from ensembldb
##
##***********************************************************************
setMethod("where", signature(object="BasicFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## just call the plain method without database.
              return(ensembldb:::.where(object))
          })
setMethod("where", signature(object="BasicFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              ## just call the plain method without database.
              return(ensembldb:::.where(object))
          })


##***********************************************************************
##
##     Implementations for GeneidFilter.
##
##     overwriting/implementation of methods for GeneidFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="GeneidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GeneidFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="GeneidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GeneidFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for GenenameFilter.
##
##     overwriting/implementation of methods for GenenameFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="GenenameFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GenenameFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="GenenameFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables,...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GenenameFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables,...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for GenebiotypeFilter.
##
##     overwriting/implementation of methods for GenebiotypeFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="GenebiotypeFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GenebiotypeFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="GenebiotypeFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GenebiotypeFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })



##***********************************************************************
##
##     Implementations for TxidFilter.
##
##     overwriting/implementation of methods for TxidFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="TxidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="TxidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for TxbiotypeFilter.
##
##     overwriting/implementation of methods for TxbiotypeFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="TxbiotypeFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxbiotypeFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="TxbiotypeFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxbiotypeFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for ExonidFilter.
##
##     overwriting/implementation of methods for ExonidFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="ExonidFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables,  ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ExonidFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables,  ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="ExonidFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ExonidFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqnameFilter.
##
##     overwriting/implementation of methods for SeqnameFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqnameFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqnameFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqnameFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqnameFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqstrandFilter.
##
##     overwriting/implementation of methods for SeqstrandFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqstrandFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqstrandFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqstrandFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqstrandFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqstartFilter.
##
##     overwriting/implementation of methods for SeqstartFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqstartFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqstartFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              feat <- object@feature
              feat <- match.arg(feat, c("mat_mirna", "pre_mirna"))
              return(unlist(prefixColumns(db,
                                          paste0(feat, "_seq_start"),
                                          with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqstartFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqstartFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              return(paste0(column(object, db, with.tables=with.tables),
                            object@condition, object@value))
          })


##***********************************************************************
##
##     Implementations for SeqendFilter.
##
##     overwriting/implementation of methods for SeqendFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqendFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqendFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              feat <- object@feature
              feat <- match.arg(feat, c("mat_mirna", "pre_mirna"))
              return(unlist(prefixColumns(db,
                                          paste0(feat, "_seq_end"),
                                          with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqendFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqendFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              return(paste0(column(object, db, with.tables=with.tables),
                            object@condition, object@value))
          })




