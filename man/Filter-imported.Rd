\name{GeneIdFilter-class}
\Rdversion{1.1}
\docType{class}
% defined in ensembldb
% column, Filter, MirhostDb, character
\alias{column,GeneIdFilter,MirhostDb,character-method}
\alias{column,GenenameFilter,MirhostDb,character-method}
\alias{column,GeneBiotypeFilter,MirhostDb,character-method}
\alias{column,EntrezgeneFilter,MirhostDb,character-method}
\alias{column,TxIdFilter,MirhostDb,character-method}
\alias{column,TxBiotypeFilter,MirhostDb,character-method}
\alias{column,ExonIdFilter,MirhostDb,character-method}
\alias{column,SeqNameFilter,MirhostDb,character-method}
\alias{column,SeqStrandFilter,MirhostDb,character-method}
\alias{column,SeqStartFilter,MirhostDb,character-method}
\alias{column,SeqEndFilter,MirhostDb,character-method}
% where, Filter, MirhostDb, character
\alias{where,GeneIdFilter,MirhostDb,character-method}
\alias{where,GenenameFilter,MirhostDb,character-method}
\alias{where,GeneBiotypeFilter,MirhostDb,character-method}
\alias{where,EntrezgeneFilter,MirhostDb,character-method}
\alias{where,TxIdFilter,MirhostDb,character-method}
\alias{where,TxBiotypeFilter,MirhostDb,character-method}
\alias{where,ExonIdFilter,MirhostDb,character-method}
\alias{where,SeqNameFilter,MirhostDb,character-method}
\alias{where,SeqStrandFilter,MirhostDb,character-method}
\alias{where,SeqStartFilter,MirhostDb,character-method}
\alias{where,SeqEndFilter,MirhostDb,character-method}
\alias{where,AnnotationFilter,MirhostDb,character-method}
% column, Filter, MirhostDb, missing
\alias{column,GeneIdFilter,MirhostDb,missing-method}
\alias{column,GenenameFilter,MirhostDb,missing-method}
\alias{column,GeneBiotypeFilter,MirhostDb,missing-method}
\alias{column,EntrezgeneFilter,MirhostDb,missing-method}
\alias{column,TxIdFilter,MirhostDb,missing-method}
\alias{column,TxBiotypeFilter,MirhostDb,missing-method}
\alias{column,ExonIdFilter,MirhostDb,missing-method}
\alias{column,SeqNameFilter,MirhostDb,missing-method}
\alias{column,SeqStrandFilter,MirhostDb,missing-method}
\alias{column,SeqStartFilter,MirhostDb,missing-method}
\alias{column,SeqEndFilter,MirhostDb,missing-method}
% where, Filter, MirhostDb, missing
\alias{where,GeneIdFilter,MirhostDb,missing-method}
\alias{where,GenenameFilter,MirhostDb,missing-method}
\alias{where,GeneBiotypeFilter,MirhostDb,missing-method}
\alias{where,EntrezgeneFilter,MirhostDb,missing-method}
\alias{where,TxIdFilter,MirhostDb,missing-method}
\alias{where,TxBiotypeFilter,MirhostDb,missing-method}
\alias{where,ExonIdFilter,MirhostDb,missing-method}
\alias{where,SeqNameFilter,MirhostDb,missing-method}
\alias{where,SeqStrandFilter,MirhostDb,missing-method}
\alias{where,SeqStartFilter,MirhostDb,missing-method}
\alias{where,SeqEndFilter,MirhostDb,missing-method}
\alias{where,AnnotationFilter,MirhostDb,missing-method}
% column, Filter, missing, missing
\alias{column,GeneIdFilter,missing,missing-method}
\alias{column,GenenameFilter,missing,missing-method}
\alias{column,GeneBiotypeFilter,missing,missing-method}
\alias{column,EntrezgeneFilter,missing,missing-method}
\alias{column,TxIdFilter,missing,missing-method}
\alias{column,TxBiotypeFilter,missing,missing-method}
\alias{column,ExonIdFilter,missing,missing-method}
\alias{column,SeqNameFilter,missing,missing-method}
\alias{column,SeqStrandFilter,missing,missing-method}
\alias{column,SeqStartFilter,missing,missing-method}
\alias{column,SeqEndFilter,missing,missing-method}
\alias{column,AnnotationFilter,missing,missing-method}
% where, Filter, missing, missing
\alias{where,GeneIdFilter,missing,missing-method}
\alias{where,GenenameFilter,missing,missing-method}
\alias{where,GeneBiotypeFilter,missing,missing-method}
\alias{where,EntrezgeneFilter,missing,missing-method}
\alias{where,TxIdFilter,missing,missing-method}
\alias{where,TxBiotypeFilter,missing,missing-method}
\alias{where,ExonIdFilter,missing,missing-method}
\alias{where,SeqNameFilter,missing,missing-method}
\alias{where,SeqStrandFilter,missing,missing-method}
\alias{where,SeqStartFilter,missing,missing-method}
\alias{where,SeqEndFilter,missing,missing-method}
\alias{where,AnnotationFilter,missing,missing-method}
%
\alias{GeneIdFilter}
\alias{GenenameFilter}
\alias{GeneBiotypeFilter}
\alias{EntrezFilter}
\alias{TxIdFilter}
\alias{TxBiotypeFilter}
\alias{ExonIdFilter}
\alias{SeqNameFilter}
\alias{SeqStrandFilter}
\alias{SeqStartFilter}
\alias{SeqEndFilter}

\title{Filters imported from ensembldb}
\description{
  These filters are imported from the \code{ensembldb} package and can
  also be used to specify which entries (i.e. genes, transcripts)
  should be retrieved from the database.
}
\section{Objects from the Class}{
  While objects can be created by calls e.g. of the form
  \code{new("GeneIdFilter", ...)} users are strongly encouraged to use
  the specific functions: \code{\link[AnnotationFilter]{GeneIdFilter}},
  \code{\link[AnnotationFilter]{GenenameFilter}},
  \code{\link[AnnotationFilter]{GeneBiotypeFilter}},
  \code{\link[AnnotationFilter]{EntrezFilter}},
  \code{\link[AnnotationFilter]{TxIdFilter}},
  \code{\link[AnnotationFilter]{ExonIdFilter}},
  \code{\link[AnnotationFilter]{SeqNameFilter}},
  \code{\link[AnnotationFilter]{SeqStrandFilter}}
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
  Class \code{"AnnotationFilter"} defined in the \code{AnnotationFilter} package.
}
\section{Methods}{
  \describe{

    \item{column}{
      \code{signature(object = "GeneIdFilter", db = "MirhostDb",
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
      \code{signature(object = "GeneIdFilter", db = "MirhostDb",
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
    \item{\code{GeneBiotypeFilter}}{
      Allows to filter on the biotype of the host gene. Use method
      \code{\link{listGenebiotypes}} to get an overview of the available
      choices.
    }

    \item{\code{GeneIdFilter}}{
      Allows to filter results based on the gene ID of the host gene.
    }

    \item{\code{GenenameFilter}}{
      Allows to filter results based on the gene name (symbol) of the
      host gene.
    }

    \item{\code{EntrezFilter}}{
      Allows to filter on the NCBI Entrezgene ID of the host gene.
    }

    \item{\code{TxBiotypeFilter}}{
      Allows to filter on the biotype of the host transcript. Use the
      method \code{\link{listTxbiotypes}} to get an overview of
      available choices.
    }

    \item{\code{TxIdFilter}}{
      Allows to filter on the transcript ID of the host gene.
    }

    \item{\code{ExonIdFilter}}{
      Filter results based on the ID of the exon in which a pre-miRNA is
      encoded (only for exonic miRNAs).
    }

    \item{\code{SeqEndFilter}}{
      Filter results based on the chromosomal end coordinate of mature
      miRNAs (if \code{feature="mat_mirna"}) or pre-miRNAs (if
      \code{feature="pre_mirna"}).
    }

    \item{\code{SeqNameFilter}}{
      Filter results based on the seqname (chromosome name) on which the
      miRNAs are encoded.
    }

    \item{\code{SeqStrandFilter}}{
      Filter results based on the strand on which the miRNAs are encoded.
    }

    \item{\code{SeqStartFilter}}{
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
  \code{\link[AnnotationFilter]{GeneIdFilter}}, \code{\link{PreMirnaFilter}}
}
\examples{


library(MirhostDb.Hsapiens.v75.v20)

## create a filter that could be used to retrieve all informations for
## the respective gene.
Gif <- GeneIdFilter("ENSG00000012817")
Gif
## returns the where clause of the SQL querys specific for the
## mirhostgenes database
where(Gif, MirhostDb.Hsapiens.v75.v20)

## create a filter for a chromosomal end position of a gene
Sef <- SeqEndFilter(10000, condition=">", "gene")
Sef

## in the miRNA host genes database we are rather interested to get
## pre-miRNAs from certain genomic positions, thus we define "pre_mirna"
## as the feature and call where by supplying also the MirhostDb
## object
Sef <- SeqEndFilter(1000, condition=">", "pre_mirna")
where(Sef, MirhostDb.Hsapiens.v75.v20)


}
\keyword{classes}

