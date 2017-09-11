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
##     Implementations for AlignmentIdFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="AlignmentIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=AlignmentIdFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="AlignmentIdFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="AlignmentIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_algn_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="AlignmentIdFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="AlignmentIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="AlignmentIdFilter", db="MirhostDb",
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
##     Implementations for ProbesetIdFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="ProbesetIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=ProbesetIdFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="ProbesetIdFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ProbesetIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  warning("Can not use a ProbesetIdFilter, as the database",
                          " doesn't have the required table.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "probeset_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="ProbesetIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="ProbesetIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ProbesetIdFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_array_features){
                  ## means we don't have that data available...
                  warning("Can not use a ProbesetIdFilter, as the database",
                          " doesn't have the required table.\n")
                  return(NULL)
              }
              suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for MatMirnaConfidence
##
##
##***********************************************************************
setMethod("column", signature(object="MatMirnaConfidence", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatMirnaConfidence, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="MatMirnaConfidence", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MatMirnaConfidence", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_matmirna_confidence){
                  warning("Can not use MatMirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "mat_mirna_confidence", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatMirnaConfidence", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatMirnaConfidence", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatMirnaConfidence", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_matmirna_confidence){
                  ## means we don't have that data available...
                  warning("Can not use MatMirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              ## check manually what values we've got. Thus far we support only
              ## confidence of value "high"
              if(object@value!="high")
                  stop("MatMirnaConfidence can only take the value 'high' for",
                       " db being a MirhostDb database.")
              if(object@value == "high"){
                  suff <- paste0(.condition(object), 1)
              }
              ## suff <- where(object) ## that returns the result from the BasicFilter...
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for PreMirnaConfidence
##
##
##***********************************************************************
setMethod("column", signature(object="PreMirnaConfidence", db="missing",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PreMirnaConfidence,",
                      " db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="PreMirnaConfidence", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PreMirnaConfidence", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_premirna_confidence){
                  warning("Can not use PreMirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              return(unlist(prefixColumns(db, "pre_mirna_confidence", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PreMirnaConfidence", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PreMirnaConfidence", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PreMirnaConfidence", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              if(!db@have_premirna_confidence){
                  ## means we don't have that data available...
                  warning("Can not use PreMirnaConfidence, as the database",
                          " doesn't have the required information.\n")
                  return(NULL)
              }
              ## check manually what values we've got. Thus far we support only
              ## confidence of value "high"
              if(object@value!="high")
                  stop("PreMirnaConfidence can only take the value 'high' for",
                       " db being a MirhostDb database.")
              if(object@value == "high"){
                  suff <- paste0(.condition(object), 1)
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
## setReplaceMethod("value", "ReadCountFilter", function(x, value){
##     ## Complain if value is not numeric.
##     if(!is.numeric(value)){
##         suppressWarnings(
##             numVal <- as.numeric(value)
##         )
##         if(is.na(numVal))
##             stop("ReadCountFilter only supports numeric values!")
##     }
##     x@value <- as.character(value)
##     return(x)
## })

##***********************************************************************
##
##     Implementations for PreMirnaFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="PreMirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PreMirnaFilter, db=missing not implemented!")
              return(NULL)
          })
setMethod("column", signature(object="PreMirnaFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PreMirnaFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_name", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PreMirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## without a database we're just calling the where of BasicFilter
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PreMirnaFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PreMirnaFilter", db="MirhostDb", with.tables="character"),
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
##     Implementations for PreMirnaIdFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="PreMirnaIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=PreMirnaIdFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="PreMirnaIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="PreMirnaIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "pre_mirna_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="PreMirnaIdFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="PreMirnaIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="PreMirnaIdFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- where(object)
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for MatMirnaFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MatMirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatMirnaFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MatMirnaFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MatMirnaFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mat_mirna_name", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatMirnaFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatMirnaFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatMirnaFilter", db="MirhostDb",
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
##     Implementations for MatMirnaIdFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MatMirnaIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MatMirnaIdFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MatMirnaIdFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, tn))
          })
setMethod("column", signature(object="MatMirnaIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mat_mirna_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MatMirnaIdFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MatMirnaIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MatMirnaIdFilter", db="MirhostDb",
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
##     Implementations for MirfamIdFilter.
##
##***********************************************************************
##
##     public methods
##***********************************************************************
setMethod("column", signature(object="MirfamIdFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              warning("Method column for object=MirfamIdFilter, db=missing not implemented!")
              return(NA)
          })
setMethod("column", signature(object="MirfamIdFilter", db="MirhostDb",
                              with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="MirfamIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, "mirfam_id", with.tables=with.tables),
                            use.names=FALSE))
          })
## that's the version without a database specified.
setMethod("where", signature(object="MirfamIdFilter", db="missing",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              return(callNextMethod())
          }
          )
setMethod("where", signature(object="MirfamIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="MirfamIdFilter", db="MirhostDb",
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
setMethod("where", signature(object="AnnotationFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## just call the plain method without database.
              .where(object)
          })
setMethod("where", signature(object="AnnotationFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              ## just call the plain method without database.
              .where(object)
          })
setMethod("where", signature(object = "AnnotationFilter", db = "missing",
                             with.tables = "missing"),
          function(object, db, with.tabbles, ...) {
              .where(object)
          })
setMethod("column", signature(object = "AnnotationFilter", db = "missing",
                              with.tables = "missing"),
          function(object, db, with.tables, ...){
              field(object)
          })

.where <- function(x) {
    ## field condition value
    paste0(.condition(x), .value(x))
}
.condition <- function(x) {
    cond <- condition(x)
    if (length(unique(value(x))) > 1) {
        if (cond == "==")
            cond <- "in "
        if (cond == "!=")
            cond <- "not in "
    }
    if (cond == "==")
        cond <- "="
    if (cond %in% c("startsWith", "endsWith", "contains"))
        cond <- "like "
    cond    
}
.value <- function(x) {
    vals <- unique(value(x))
    if (is(x, "CharacterFilter")) {
        vals <- sQuote(gsub(unique(vals), pattern = "'", replacement = "''"))
    }
    if (length(vals) > 1)
        vals <- paste0("(",  paste0(vals, collapse = ","), ")")
    ## Process the like/startsWith/endsWith
    if (condition(x) == "startsWith")
        vals <- paste0("'", unique(x@value), "%'")
    if (condition(x) == "endsWith")
        vals <- paste0("'%", unique(x@value), "'")
    if (condition(x) == "contains")
        vals <- paste0("'%", unique(x@value), "%'")
    vals
}

##***********************************************************************
##
##     Implementations for GeneIdFilter.
##
##     overwriting/implementation of methods for GeneIdFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="GeneIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GeneIdFilter", db="MirhostDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="GeneIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GeneIdFilter", db="MirhostDb", with.tables="character"),
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
setMethod("column", signature(object = "GenenameFilter", db = "missing",
                              with.tables = "missing"),
          function(object, db, with.tables, ...)
              "gene_name"
          )
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
##     Implementations for GeneBiotypeFilter.
##
##     overwriting/implementation of methods for GeneBiotypeFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="GeneBiotypeFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GeneBiotypeFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="GeneBiotypeFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GeneBiotypeFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })



##***********************************************************************
##
##     Implementations for TxIdFilter.
##
##     overwriting/implementation of methods for TxIdFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="TxIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="TxIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxIdFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for TxBiotypeFilter.
##
##     overwriting/implementation of methods for TxBiotypeFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="TxBiotypeFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxBiotypeFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="TxBiotypeFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxBiotypeFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for ExonIdFilter.
##
##     overwriting/implementation of methods for ExonIdFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="ExonIdFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables,  ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ExonIdFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables,  ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="ExonIdFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ExonIdFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqNameFilter.
##
##     overwriting/implementation of methods for SeqNameFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqNameFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqNameFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqNameFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqNameFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqStrandFilter.
##
##     overwriting/implementation of methods for SeqStrandFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqStrandFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqStrandFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("where", signature(object="SeqStrandFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqStrandFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })


##***********************************************************************
##
##     Implementations for SeqStartFilter.
##
##     overwriting/implementation of methods for SeqStartFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqStartFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqStartFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              feat <- object@feature
              feat <- match.arg(feat, c("mat_mirna", "pre_mirna"))
              return(unlist(prefixColumns(db,
                                          paste0(feat, "_seq_start"),
                                          with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("column", signature(object = "SeqStartFilter", db = "missing",
                              with.tables = "missing"),
          function(object, db, with.tables, ...) {
              paste0(object@feature, "_seq_start")
          })
setMethod("where", signature(object="SeqStartFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqStartFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              return(paste0(column(object, db, with.tables=with.tables),
                            object@condition, object@value))
          })


##***********************************************************************
##
##     Implementations for SeqEndFilter.
##
##     overwriting/implementation of methods for SeqEndFilter defined
##     in makeEnsemblDb to be applicable in mirhostgenesdb
##
##***********************************************************************
setMethod("column", signature(object="SeqEndFilter", db="MirhostDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(column(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqEndFilter", db="MirhostDb",
                              with.tables="character"),
          function(object, db, with.tables, ...){
              feat <- object@feature
              feat <- match.arg(feat, c("mat_mirna", "pre_mirna"))
              return(unlist(prefixColumns(db,
                                          paste0(feat, "_seq_end"),
                                          with.tables=with.tables),
                            use.names=FALSE))
          })
setMethod("column", signature(object = "SeqEndFilter", db = "missing",
                              with.tables = "missing"),
          function(object, db, with.tables, ...) {
              paste0(object@feature, "_seq_end")
          })
setMethod("where", signature(object="SeqEndFilter", db="MirhostDb",
                             with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- tableNames(db)
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqEndFilter", db="MirhostDb",
                             with.tables="character"),
          function(object, db, with.tables, ...){
              return(paste0(column(object, db, with.tables=with.tables),
                            object@condition, object@value))
          })




