\name{GeneidFilter-class}
\Rdversion{1.1}
\docType{class}
% defined in ensembldb
% column, Filter, MirhostDb, character
\alias{column,GeneidFilter,MirhostDb,character-method}
\alias{column,GenenameFilter,MirhostDb,character-method}
\alias{column,GenebiotypeFilter,MirhostDb,character-method}
\alias{column,EntrezgeneFilter,MirhostDb,character-method}
\alias{column,TxidFilter,MirhostDb,character-method}
\alias{column,TxbiotypeFilter,MirhostDb,character-method}
\alias{column,ExonidFilter,MirhostDb,character-method}
\alias{column,SeqnameFilter,MirhostDb,character-method}
\alias{column,SeqstrandFilter,MirhostDb,character-method}
\alias{column,SeqstartFilter,MirhostDb,character-method}
\alias{column,SeqendFilter,MirhostDb,character-method}
% where, Filter, MirhostDb, character
\alias{where,GeneidFilter,MirhostDb,character-method}
\alias{where,GenenameFilter,MirhostDb,character-method}
\alias{where,GenebiotypeFilter,MirhostDb,character-method}
\alias{where,EntrezgeneFilter,MirhostDb,character-method}
\alias{where,TxidFilter,MirhostDb,character-method}
\alias{where,TxbiotypeFilter,MirhostDb,character-method}
\alias{where,ExonidFilter,MirhostDb,character-method}
\alias{where,SeqnameFilter,MirhostDb,character-method}
\alias{where,SeqstrandFilter,MirhostDb,character-method}
\alias{where,SeqstartFilter,MirhostDb,character-method}
\alias{where,SeqendFilter,MirhostDb,character-method}
\alias{where,BasicFilter,MirhostDb,character-method}
% column, Filter, MirhostDb, missing
\alias{column,GeneidFilter,MirhostDb,missing-method}
\alias{column,GenenameFilter,MirhostDb,missing-method}
\alias{column,GenebiotypeFilter,MirhostDb,missing-method}
\alias{column,EntrezgeneFilter,MirhostDb,missing-method}
\alias{column,TxidFilter,MirhostDb,missing-method}
\alias{column,TxbiotypeFilter,MirhostDb,missing-method}
\alias{column,ExonidFilter,MirhostDb,missing-method}
\alias{column,SeqnameFilter,MirhostDb,missing-method}
\alias{column,SeqstrandFilter,MirhostDb,missing-method}
\alias{column,SeqstartFilter,MirhostDb,missing-method}
\alias{column,SeqendFilter,MirhostDb,missing-method}
% where, Filter, MirhostDb, missing
\alias{where,GeneidFilter,MirhostDb,missing-method}
\alias{where,GenenameFilter,MirhostDb,missing-method}
\alias{where,GenebiotypeFilter,MirhostDb,missing-method}
\alias{where,EntrezgeneFilter,MirhostDb,missing-method}
\alias{where,TxidFilter,MirhostDb,missing-method}
\alias{where,TxbiotypeFilter,MirhostDb,missing-method}
\alias{where,ExonidFilter,MirhostDb,missing-method}
\alias{where,SeqnameFilter,MirhostDb,missing-method}
\alias{where,SeqstrandFilter,MirhostDb,missing-method}
\alias{where,SeqstartFilter,MirhostDb,missing-method}
\alias{where,SeqendFilter,MirhostDb,missing-method}
\alias{where,BasicFilter,MirhostDb,missing-method}
% column, Filter, missing, missing
\alias{column,GeneidFilter,missing,missing-method}
\alias{column,GenenameFilter,missing,missing-method}
\alias{column,GenebiotypeFilter,missing,missing-method}
\alias{column,EntrezgeneFilter,missing,missing-method}
\alias{column,TxidFilter,missing,missing-method}
\alias{column,TxbiotypeFilter,missing,missing-method}
\alias{column,ExonidFilter,missing,missing-method}
\alias{column,SeqnameFilter,missing,missing-method}
\alias{column,SeqstrandFilter,missing,missing-method}
\alias{column,SeqstartFilter,missing,missing-method}
\alias{column,SeqendFilter,missing,missing-method}
% where, Filter, missing, missing
\alias{where,GeneidFilter,missing,missing-method}
\alias{where,GenenameFilter,missing,missing-method}
\alias{where,GenebiotypeFilter,missing,missing-method}
\alias{where,EntrezgeneFilter,missing,missing-method}
\alias{where,TxidFilter,missing,missing-method}
\alias{where,TxbiotypeFilter,missing,missing-method}
\alias{where,ExonidFilter,missing,missing-method}
\alias{where,SeqnameFilter,missing,missing-method}
\alias{where,SeqstrandFilter,missing,missing-method}
\alias{where,SeqstartFilter,missing,missing-method}
\alias{where,SeqendFilter,missing,missing-method}
\alias{where,BasicFilter,missing,missing-method}
%
\alias{GeneidFilter}
\alias{GenenameFilter}
\alias{GenebiotypeFilter}
\alias{EntrezidFilter}
\alias{TxidFilter}
\alias{TxbiotypeFilter}
\alias{ExonidFilter}
\alias{SeqnameFilter}
\alias{SeqstrandFilter}
\alias{SeqstartFilter}
\alias{SeqendFilter}

\title{Filters imported from ensembldb}
\description{
  These filters are imported from the \code{ensembldb} package and can
  also be used to specify which entries (i.e. genes, transcripts)
  should be retrieved from the database.
}
\section{Objects from the Class}{
  While objects can be created by calls e.g. of the form
  \code{new("GeneidFilter", ...)} users are strongly encouraged to use
  the specific functions: \code{\link[ensembldb]{GeneidFilter}},
  \code{\link[ensembldb]{GenenameFilter}},
  \code{\link[ensembldb]{GenebiotypeFilter}},
  \code{\link[ensembldb]{EntrezidFilter}},
  \code{\link[ensembldb]{TxidFilter}},
  \code{\link[ensembldb]{ExonidFilter}},
  \code{\link[ensembldb]{SeqnameFilter}},
  \code{\link[ensembldb]{SeqstrandFilter}},
  \code{\link[ensembldb]{SeqstartFilter}} and
  \code{\link[ensembldb]{SeqendFilter}}.
  See examples below for usage.
}
\section{Slots}{
  \describe{
    \item{\code{condition}:}{
      Object of class \code{"character"}: can be
      either \code{"="}, \code{"in"} or \code{"like"} to filter on character values
      (e.g. gene id, gene biotype, seqname etc), or \code{"="}, \code{">"}
      or \code{"<"} for numerical values (chromosome/seq
      coordinates). Note that for \code{"like"} \code{value} should be a
      SQL pattern (e.g. \code{"ENS\%"}).
    }

    \item{\code{value}:}{
      Object of class \code{"character"}: the value
      to be used for filtering.
    }

  }
}
\section{Extends}{
  Class \code{"BasicFilter"} defined in the \code{ensembldb} package, directly.
}
\section{Methods}{
  \describe{

    \item{column}{
      \code{signature(object = "GeneidFilter", db = "MirhostDb",
      with.tables="character")}: returns the column (table column name
      in the miRNA host genes database) to be used for the
      filtering. The method can be applied to any filter class defined
      in the \code{mirhostgenes} or \code{ensembldb} packages. Argument
      \code{with.tables} allows to restrict the filter to be applied to
      certain database tables. By default the filter is applied to all
      database tables.

      The method returns a character string with the name of the column.
    }

    \item{where}{
      \code{signature(object = "GeneidFilter", db = "MirhostDb",
      with.tables="character")}: returns the where clause for the SQL
      call to the miRNA host genes database. The method can be applied
      to any filter class defined in the \code{mirhostgenes} or
      \code{ensembldb} packages. Argument \code{with.tables} allows to
      restrict the filter to be applied to certain database tables. By
      default the filter is applied to all database tables.

      The method returns the SQL \emph{where} clause that will be used
      to restrict the results returned by the database. This string is
      dependent on the type of database.
    }

  }
}
\details{
  Description of filters (in alphabetical order).

  \describe{
    \item{\code{GenebiotypeFilter}}{
      Allows to filter on the biotype of the host gene. Use method
      \code{\link{listGenebiotypes}} to get an overview of the available
      choices.
    }

    \item{\code{GeneidFilter}}{
      Allows to filter results based on the gene ID of the host gene.
    }

    \item{\code{GenenameFilter}}{
      Allows to filter results based on the gene name (symbol) of the
      host gene.
    }

    \item{\code{EntrezidFilter}}{
      Allows to filter on the NCBI Entrezgene ID of the host gene.
    }

    \item{\code{TxbiotypeFilter}}{
      Allows to filter on the biotype of the host transcript. Use the
      method \code{\link{listTxbiotypes}} to get an overview of
      available choices.
    }

    \item{\code{TxidFilter}}{
      Allows to filter on the transcript ID of the host gene.
    }

    \item{\code{ExonidFilter}}{
      Filter results based on the ID of the exon in which a pre-miRNA is
      encoded (only for exonic miRNAs).
    }

    \item{\code{SeqendFilter}}{
      Filter results based on the chromosomal end coordinate of mature
      miRNAs (if \code{feature="mat_mirna"}) or pre-miRNAs (if
      \code{feature="pre_mirna"}).
    }

    \item{\code{SeqnameFilter}}{
      Filter results based on the seqname (chromosome name) on which the
      miRNAs are encoded.
    }

    \item{\code{SeqstrandFilter}}{
      Filter results based on the strand on which the miRNAs are encoded.
    }

    \item{\code{SeqstartFilter}}{
      Filter results based on the chromosomal start coordinate of mature
      miRNAs (if \code{feature="mat_mirna"}) or pre-miRNAs (if
      \code{feature="pre_mirna"}).
    }
  }
}
\value{
  See method descriptions above for the return values.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link[ensembldb]{GeneidFilter}}, \code{\link{PremirnaFilter}}
}
\examples{


library(MirhostDb.Hsapiens.v75.v20)

## create a filter that could be used to retrieve all informations for
## the respective gene.
Gif <- GeneidFilter("ENSG00000012817")
Gif
## returns the where clause of the SQL querys specific for the
## mirhostgenes database
where(Gif, MirhostDb.Hsapiens.v75.v20)

## create a filter for a chromosomal end position of a gene
Sef <- SeqendFilter(10000, condition=">", "gene")
Sef

## in the miRNA host genes database we are rather interested to get
## pre-miRNAs from certain genomic positions, thus we define "pre_mirna"
## as the feature and call where by supplying also the MirhostDb
## object
Sef <- SeqendFilter(1000, condition=">", "pre_mirna")
where(Sef, MirhostDb.Hsapiens.v75.v20)


}
\keyword{classes}

