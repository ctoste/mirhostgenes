\name{PositionFilter}
\alias{AlignmentIdFilter}
\alias{ArrayFilter}
\alias{DatabaseFilter}
\alias{PositionFilter}
\alias{PreMirnaConfidence}
\alias{PreMirnaFilter}
\alias{PreMirnaIdFilter}
\alias{ProbesetIdFilter}
\alias{MatMirnaConfidence}
\alias{MatMirnaFilter}
\alias{MatMirnaIdFilter}
\alias{MirfamFilter}
\alias{MirfamIdFilter}
\alias{ReadCountFilter}
\title{
  Constructor functions for filter objects
}
\description{
  These functions allow to create filter objects that can be used to
  retrieve specific elements from the annotation database.
}
\usage{

AlignmentIdFilter(value, condition = "==")

ArrayFilter(value, condition = "==")

DatabaseFilter(value, condition = "==")

MatMirnaConfidence(value="high", condition="==")

MatMirnaFilter(value, condition = "==", match.case=TRUE)

MatMirnaIdFilter(value, condition = "==")

MirfamFilter(value, condition = "==", match.case=TRUE)

MirfamIdFilter(value, condition = "==")

PositionFilter(value, condition = "==")

PreMirnaConfidence(value="high", condition="==")

PreMirnaFilter(value, condition = "==", match.case=TRUE)

PreMirnaIdFilter(value, condition = "==")

ProbesetIdFilter(value, condition = "==")

ReadCountFilter(value=0, condition=">", of="pre_mirna")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{condition}{
    The condition to be used in the comparison. For character values
    \code{"="}, \code{"in"} and \code{"like"} are allowed, for numeric values
    (\code{SeqstartFilter} and \code{SeqendFilter}) \code{"="},
    \code{">"}, \code{">="}, \code{"<"} and \code{"<="}. Note that for
    \code{"like"} \code{value} should be a SQL pattern (e.g. \code{"ENS\%"}).
  }

  \item{match.case}{
    For \code{PreMirnaFilter}, \code{MatMirnaFilter}
    and \code{MirfamFilter}: whether the query should be performed case
    sensitive or insensitive.
  }

  \item{of}{
    For \code{ReadCountFilter}: whether the read count filter should be
    applied to pre-miRNAs (\code{"pre_mirna"}) or mature miRNAs
    (\code{"mat_mirna"}).
  }

  \item{value}{
    The filter value, e.g.,
    for \code{GeneIdFilter} the id of the feature for which the data
    should be retrieved. For \code{PositionFilter} either
    \code{"exonic"}, \code{"intronic"} or \code{"both"} are allowed,
    specifying that the (pre-)miRNA should be located in an exon, intron
    or within a host transcript.  Can be any string for
    \code{PreMirnaFilter}, \code{PreMirnaIdFilter},
    \code{MatMirnaFilter} and \code{MatMirnaIdFilter}.
    For \code{MatMirnaConfidence} and \code{PreMirnaConfidence} at
    present only a value of \code{"high"} is allowed, which allows to
    retrieve pre-miRNAs or mature miRNAs specified as \emph{high confidence}
    in the miRBase.
  }

}


\value{
  Depending on the function called, an instance of:
  \code{\linkS4class{ArrayFilter}},
  \code{\linkS4class{AlignmentIdFilter}},
  \code{\linkS4class{DatabaseFilter}},
  \code{\linkS4class{MatMirnaConfidence}},
  \code{\linkS4class{MatMirnaFilter}},
  \code{\linkS4class{MatMirnaIdFilter}},
  \code{\linkS4class{MirfamFilter}},
  \code{\linkS4class{MirfamIdFilter}},
  \code{\linkS4class{PositionFilter}},
  \code{\linkS4class{PreMirnaConfidence}},
  \code{\linkS4class{PreMirnaFilter}},
  \code{\linkS4class{PreMirnaIdFilter}},
  \code{\linkS4class{ProbesetIdFilter}},
  \code{\linkS4class{ReadCountFilter}}
}

\author{ Johannes Rainer }

\seealso{ \code{\link{hostgenes}}, \code{\link{hosttx}},
  \code{\linkS4class{ArrayFilter}},
  \code{\linkS4class{AlignmentIdFilter}},
  \code{\linkS4class{DatabaseFilter}},
  \code{\linkS4class{PositionFilter}},
  \code{\linkS4class{PreMirnaFilter}},
  \code{\linkS4class{PreMirnaIdFilter}},
  \code{\linkS4class{ProbesetIdFilter}},
  \code{\linkS4class{MatMirnaFilter}},
  \code{\linkS4class{MatMirnaIdFilter}}
  \code{\link[AnnotationFilter]{GeneIdFilter}},
  \code{\link[AnnotationFilter]{GenenameFilter}},
  \code{\link[AnnotationFilter]{GeneBiotypeFilter}},
  \code{\link[AnnotationFilter]{TxIdFilter}},
  \code{\link[AnnotationFilter]{TxBiotypeFilter}},
  \code{\link[AnnotationFilter]{ExonIdFilter}},
  \code{\link[AnnotationFilter]{SeqNameFilter}},
  \code{\link[AnnotationFilter]{SeqStrandFilter}}
}
\examples{

## Note: for additional examples see the help of the
## hostgenes, hosttx etc. functions.

library(MirhostDb.Hsapiens.v75.v20)

## create a filter that could be used to retrieve all host genes with
## a pre-miRNA encoded in one of their exons.
Pos <- PositionFilter("exonic")

where(Pos, MirhostDb.Hsapiens.v75.v20)
column(Pos, MirhostDb.Hsapiens.v75.v20)

## Next we define a SeqstartFilter (defined in the AnnotationFilter
## package)
Start <- SeqStartFilter(12333)

## the default where and column calls will return us the settings
## required for the makeEnsemblDb package
where(Start)
column(Start)

## however, that would not work in a mirhostgenes package since we're
## interested in the genomic alignment of miRNAs, not genes. Thus,
## submitting in addition the instance of the mirhostgenes database
## returns the correct columns and where clause:
Start <- SeqStartFilter(12333, feature="pre_mirna")
where(Start, MirhostDb.Hsapiens.v75.v20)
column(Start, MirhostDb.Hsapiens.v75.v20)

}
\keyword{data}
