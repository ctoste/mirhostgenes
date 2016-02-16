\name{MirhostDb-class}
\Rdversion{1.1}
\docType{class}
\alias{MirhostDb-class}
\alias{dbconn}
\alias{dbconn,MirhostDb-method}
\alias{listArrays}
\alias{listArrays,MirhostDb-method}
\alias{listColumns}
\alias{listColumns,MirhostDb-method}
\alias{listDatabases}
\alias{listDatabases,MirhostDb-method}
\alias{listGenebiotypes}
\alias{listGenebiotypes,MirhostDb-method}
\alias{listTxbiotypes}
\alias{listTxbiotypes,MirhostDb-method}
\alias{metadata,MirhostDb-method}
\alias{mirbaseVersion}
\alias{mirbaseVersion,MirhostDb-method}
\alias{show}
\alias{show,MirhostDb-method}
\alias{organism}
\alias{organism,MirhostDb-method}
\alias{listTables}
\alias{listTables,MirhostDb-method}
\alias{version}
\alias{version,MirhostDb-method}
\alias{MirhostDb}

\title{miRNA host gene definitions}
\description{
  Basic functionality for miRNA host gene packages. Such database
  packages can be created using the
  \code{\link{makeMirhostgenesPackage}} function (see the
  corresponding help page for more information).

  For methods and examples to retrieve host gene definitions and similar
  from the database, refer to the \code{\link{hostgenes}} man page.

}
\section{Objects from the Class}{
  A connection to the respective annotation database is created upon
  loading of an annotation package created with the
  \code{\link{makeMirhostgenesPackage}} function.
  See examples below for basic usage.
}
\usage{

\S4method{dbconn}{MirhostDb}(x)

\S4method{listArrays}{MirhostDb}(x, ...)

\S4method{listColumns}{MirhostDb}(x, table, skip.keys=TRUE, ...)

\S4method{listTables}{MirhostDb}(x, ...)

\S4method{listDatabases}{MirhostDb}(x, ...)

\S4method{listGenebiotypes}{MirhostDb}(x, ...)

\S4method{listTxbiotypes}{MirhostDb}(x, ...)

\S4method{metadata}{MirhostDb}(x)

\S4method{mirbaseVersion}{MirhostDb}(x)

MirhostDb(x)

\S4method{organism}{MirhostDb}(object)

\S4method{version}{MirhostDb}(object, what="ensembl", ...)

}
\arguments{

  \item{...}{
    Additional arguments. For not used for any methods.
  }

  \item{object}{
    For \code{version} and \code{organism}: the \code{MirhostDb} instance for which
    the version of the miRBase or Ensembl database should be returned.
  }

  \item{skip.keys}{
    For \code{listColumns}: whether primary and foreign keys should be
    returned.
  }

  \item{table}{
    For \code{listColumns}: specify the database table name for which
    the columns should be listed.
  }

  \item{what}{
    For \code{version}: either \code{"ensembl"} or \code{"mirbase"}
    specifying the resource for which the version used to build the
    database should be returned.
  }

  \item{x}{
    For \code{MirhostDb} the file name of the SQLite database
    file with the host gene data (as created by the
    \code{\link{makeHostgeneSQLiteFromTables}} function).
    For all other methods the \code{MirhostDb} instance.
  }

}
\section{Slots}{
  \describe{
    \item{con}{
      Object of class \code{"DBIConnection"}: the
      connection to the database.
    }

    \item{tables}{
      named list of database table columns with the names being the
      database table names. The order in which the tables are returned
      are the by their degree, i.e. the number of other tables they can
      be connected to (joined with).
    }

    \item{have_pre_sequence}{
      Logical specifying whether a database table with pre-miRNA
      sequences in available.
    }
    \item{have_mirfam}{
      Logical specifying whether a database table with miRNA family
      definitions is available.
    }
    \item{have_array_features}{
      Logical specifying whether a database table with microarray
      features (probe sets) is available.
    }
    \item{have_premirna_confidence}{
      Logical specifying whether pre-miRNA confidence information, as
      provided by the miRBase database, is available.
    }
    \item{have_premirna_readcount}{
      Logical specifying whether the read counts information for
      pre-miRNAs as provided by the miRBase database is available.
    }
    \item{have_matmirna_confidence}{
      Logical specifying whether mature miRNA confidence information, as
      provided by the miRBase database, is available.
    }
    \item{have_matmirna_readcount}{
      Logical specifying whether the read counts information for
      mature miRNAs as provided by the miRBase database is available.
    }
  }

}
\section{Basic usage}{

  (in alphabetical order)

  \describe{

    \item{dbconn}{
      Returns the \code{SQLiteConnection} connection to the internal database.
    }

    \item{listArrays}{
      Returns a character vector with the names of the microarrays for
      which probe set to host transcript mappings are available.
    }

    \item{listColumns}{
      Returns a character vector of all columns stored in a database,
      or in a single database table if \code{table} was specified.
    }

    \item{listDatabases}{
      Returns a character vector of all databases from which
      gene/transcript models for host gene predictions were retrieved.
      In most instances the databases \code{"core"} (Ensembl core),
      \code{"otherfeatures"} (Ensembl otherfeatures containing CCDs and
      RefSeq genes) and \code{"vega"} (Ensembl Vega database with
      curated genes) are queried.
    }

    \item{listGenebiotypes}{
      Returns a character vector of all available gene biotypes from the
      database. These values can be used in a
      \code{\link{GenebiotypeFilter}} to filter results queried from the
      database.
    }

    \item{listTables}{
      Returns a named list of database table columns (names of the
      list being the database table names).
    }

    \item{listTxbiotypes}{
      Returns all available transcript biotypes from the database. These
      values can be used in a \code{\link{TxbiotypeFilter}} to filter
      results queried from the database.
    }

    \item{metadata}{
      Retrieves the metadata of the database. Returns a
      \code{data.frame} with the key-value pairs.
    }

    \item{mirbaseVersion}{
      Returns the version of the miRBase database the \code{MirhostDb}
      database bases on as a character string.
    }

    \item{MirhostDb}{
      The constructor function. Returns a \code{MirhostDb} instance.
    }

    \item{organism}{
      Returns the organism (e.g. \code{"Homo sapiens"}).
    }

    \item{show}{
      Displays some informations from the database.
    }

    \item{version}{
      Returns the version of the resource (either Ensembl or miRBase)
      that was used to build the database.
    }

  }
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{hostgenes}}, \code{\link{premirnas}}, \code{\link{makeMirhostgenesPackage}}, \code{\linkS4class{PositionFilter}}
}
\examples{

library(MirhostDb.Hsapiens.v75.v20)
## make a shortcut...
mhg <- MirhostDb.Hsapiens.v75.v20

## show some information of the database
mhg


## the version of the Ensembl database
ensemblVersion(mhg)

## the version of the mirbase database
mirbaseVersion(mhg)

## list all tables along with their columns
listTables(mhg)

## list all columns stored in the database for database table host_tx
listColumns(mhg, table="host_tx")

## list all databases on which gene/transcript models host genes were
## defined
listDatabases(mhg)

## for what organism is this package
organism(mhg)

## list all available gene biotypes
listGenebiotypes(mhg)

## list all available transcript biotypes
listTxbiotypes(mhg)

## list all microarrays for which probe sets have been indentified
## targeting the host genes/transcripts
listArrays(mhg)

## get an overview of the full metadata
metadata(mhg)

## get some information about the genome
genome(mhg)

## or the chromosomes on which the miRNA precursors (pre-miRNAs)
## where defined
seqinfo(mhg)

}
\keyword{classes}

