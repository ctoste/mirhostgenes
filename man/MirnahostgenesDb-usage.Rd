\name{hostgenes}
\Rdversion{1.1}
\alias{hostgenes}
\alias{hostgenes,MirhostDb-method}
\alias{hostgenesBy}
\alias{hostgenesBy,MirhostDb-method}
\alias{hosttx}
\alias{hosttx,MirhostDb-method}
\alias{hosttxBy}
\alias{hosttxBy,MirhostDb-method}
\alias{matmirnas}
\alias{matmirnas,MirhostDb-method}
\alias{matmirnasBy}
\alias{matmirnasBy,MirhostDb-method}
\alias{matmirnasInMultiplePremirnas}
\alias{matmirnasInMultiplePremirnas,MirhostDb-method}
\alias{premirnas}
\alias{premirnas,MirhostDb-method}
\alias{premirnasBy}
\alias{premirnasBy,MirhostDb-method}
\alias{premirnasWithMultipleAlignments}
\alias{premirnasWithMultipleAlignments,MirhostDb-method}
\alias{probesets}
\alias{probesets,MirhostDb-method}
\alias{probesetsBy}
\alias{probesetsBy,MirhostDb-method}


\title{Retrieving miRNA information and host gene definitions}
\description{
  Use and retrieve miRNA host gene definitions stored in a corresponding
  database. Such database packages can be created using the
  \code{\link{makeMirhostgenesPackage}} function (see the
  corresponding help page for more information).
  For some basic usage of the database and package see the
  \code{\linkS4class{MirhostDb}}.
}
\usage{

\S4method{hostgenes}{MirhostDb}(x, columns=listColumns(x, "host_gene"),
                                filter, order.by="gene_id",
                                order.type="asc", return.type="DataFrame")

\S4method{hostgenesBy}{MirhostDb}(x, by="pre_mirna_algn",
                                  columns=listColumns(x, "host_gene"), filter,
                                  return.type="DataFrame", drop.empty=TRUE,
                                  use.names=FALSE)

\S4method{hosttx}{MirhostDb}(x, columns=listColumns(x, "host_tx"), filter,
                             order.by="tx_id", order.type="asc",
                             return.type="DataFrame")

\S4method{hosttxBy}{MirhostDb}(x, by="pre_mirna_algn",
                               columns=listColumns(x, "host_tx"), filter,
                               return.type="DataFrame", drop.empty=TRUE,
                               use.names=FALSE)

\S4method{matmirnas}{MirhostDb}(x, columns=listColumns(x, "mat_mirna"),
                                filter, order.by="mat_mirna_id",
                                order.type="asc", return.type="DataFrame")

\S4method{matmirnasBy}{MirhostDb}(x, by="pre_mirna_algn",
                                  columns=listColumns(x, "mat_mirna"), filter,
                                  return.type="DataFrame", use.names=FALSE)

\S4method{matmirnasInMultiplePremirnas}{MirhostDb}(x,columns=c(listColumns(x, "mat_mirna"),
                                                   "pre_mirna_id", "pre_mirna_name"),
                                                   filter=list(),
                                                   return.type="DataFrame")

\S4method{premirnas}{MirhostDb}(x, columns=listColumns(x, "pre_mirna"), filter,
                                order.by="pre_mirna_id", order.type="asc",
                                return.type="DataFrame")

\S4method{premirnasBy}{MirhostDb}(x, by="mat_mirna",
                                  columns=listColumns(x, "pre_mirna"),
                                  filter, return.type="DataFrame",
                                  use.names=FALSE)

\S4method{premirnasWithMultipleAlignments}{MirhostDb}(x,
                                                      columns=listColumns(x, "pre_mirna"),
                                                      filter=list(),
                                                      return.type="DataFrame")

\S4method{probesets}{MirhostDb}(x, columns=listColumns(x, "array_feature"), filter,
                                order.by="probeset_id", order.type="asc",
                                return.type="DataFrame")

\S4method{probesetsBy}{MirhostDb}(x, by="pre_mirna_algn",
                                  columns=listColumns(x, "array_feature"),
                                  filter, return.type="DataFrame",
                                  drop.empty=TRUE,
                                  use.names=FALSE)

}
\arguments{
  (in alphabetical order)

  \item{by}{
    For \code{hostgenesBy}, \code{hosttxBy}, \code{matmirnasBy} and
    \code{premirnasBy}: by what the entries should be grouped
    (\code{"pre_mirna_algn"}, \code{"pre_mirna"}, \code{"mat_mirna"},
    \code{"mirfam"}, \code{"host_gene"}, \code{"host_tx"},
    \code{"database"} or \code{"probeset"} to group results by the alingmnet ID of the
    pre-miRNA, by pre-miRNA, mature miRNA, miRNA family, host gene, host
    transcript, database or microarray probe set associated with the
    host transcript, respectively). The default for all methods
    except \code{premirnasBy} is \code{pre_mirna_algn}, thus it returns
    the entries grouped by the unique pre-miRNA alignments.
    \code{by="database"} causes the entries being grouped by
    the database in which the transcript/gene model of the host
    transcript/gene was defined. To get a list of all databases use the
    \code{link{listDatabases}} method.
  }

  \item{columns}{
    Character vector of columns (attributes) to return from the database. For a
    complete list of available columns use the methods
    \code{\link{listTables}} or \code{\link{listColumns}}.
  }

  \item{drop.empty}{
    For \code{hostgenesBy} and \code{hosttxBy}: whether empty list
    elements should be dropped (default). Otherwise all elements are
    returned, also empty ones, e.g. representing pre-miRNAs for which no
    host gene or transcript was defined (if \code{by="pre_mirna"}).
  }

  \item{filter}{
    A single filter instance or \code{list} of filter instances to be
    used to fetch specific elements from the database. See help for
    \code{\linkS4class{PositionFilter}} or \code{AnnotationFilter} in package
    \code{ensembldb} for information on filter objects and their use.
  }

  %% \item{ifEmptyBy}{
  %%   For \code{matmirnasBy}, \code{premirnasBy}, \code{hostgenesBy} and
  %%   \code{hosttxBy}: allows to return database entries for which the attribute
  %%   specified with argument \code{by} is empty. Values will be combined
  %%   into the list element named according to \code{ifEmptyBy}, or not
  %%   returned (if \code{ofEmptyBy=NA}, which is the default). See notes
  %%   below for more information.
  %% }

  \item{order.by}{
    The column by which the result should be ordered. Can also be a
    string with multiple columns, separated by a \code{","}.
  }

  \item{order.type}{
    Either \code{"asc"} or \code{"desc"} depending on whether the
    results should be returned in ascending or descending order.
  }

  \item{return.type}{
    Allows to specify the class of the result object. Allowed values are
    \code{"data.frame"} or \code{"DataFrame"} (the
    default). Additionally, for methods
    \code{matmirnas}, \code{matmirnasBy}, \code{premirnas} and
    \code{premirnasBy}, \code{return.type="GRanges"} can be
    specified which returns a \code{GRanges} object for the mature miRNA
    or pre-miRNA (i.e. representing its genomic alignment) with all
    additional annotations added as metadata columns.

    Note that methods
    \code{premirnasBy}, \code{matmirnasBy}, \code{hostgenesBy}
    \code{hosttxBy} and \code{probesetsBy} split the resulty by the argument \code{by}, thus,
    \code{return.type} specifies the class of the elements in the \code{list} (for
    \code{return.type="data.frame"}) or \code{SplitDataFrameList} (for
    \code{return.type="DataFrame"}) that is returned.
  }

  \item{use.names}{
    Uses, if available, the names instead of the IDs to group elements
    (e.g. the pre-miRNA name instead of the pre-miRNA ID). Note, that
    the gene name (symbol) might be empty for some genes, thus, all
    entries for genes without a name would be grouped together.
  }

  \item{x}{
    The \code{MirhostDb} instance from which the data should be retrieved.
  }

}
\section{Retrieving miRNA information}{
  These methods allow to access all miRNA related information from the
  database (i.e. get mature miRNAs and pre-miRNAs).

  \describe{
    \item{matmirnas}{
      Returns all mature miRNAs from the database along with optional
      additional columns from other database tables (which can be
      empty for some mature miRNAs). Note that column \code{"sequence"}
      returns the actual RNA sequence of the mature miRNA, not the
      genomic DNA defined by the columns \code{"mat_mirna_seq_start"}
      and \code{"mat_mirna_seq_end"}.
      Also, be aware that mature miRNAs encoded in several pre-miRNAs or in
      pre-miRNAs with multiple genomic alignments are listed in
      multiple rows of the results table (as their start and end
      coordinates differ).
      To get only a unique list of miRNAs \code{columns} should be set to
      \code{c("mat_mirna_id", "mat_mirna_name")}.

      The method returns a \code{DataFrame}, \code{data.frame} or
      \code{GRanges} depending on the value of the \code{return.type}
      argument (\code{"DataFrame"}, \code{"data.frame"} or
      \code{"GRanges"}, respectively). Entries in the returned object
      are ordered according to the parameter \code{order.by}, NOT by any
      ordering of values in eventually submitted filter objects.
    }

    \item{matmirnasBy}{
      Returns a \code{CompressedSplitDataFrameList} of \code{DataFrame}s
      or a \code{list} of \code{data.frame}s with the names of the list
      being the ids by which the mature miRNAs are grouped (e.g. pre-miRNA
      ids) and the elements of the list being the host gene
      entries. Similar to \code{matmirnas}, column \code{"sequence"} in
      the result object contains the RNA sequence of the mature miRNA.

      The method returns a \code{SplitDataFrameList} (list of
      \code{DataFrame}s), a \code{list} of \code{data.frame}s or a
      \code{GRangesList}, depending on the value of the parameter
      \code{return.type} (\code{"DataFrame"}, \code{"data.frame"} or
      \code{"GRanges"}, respectively). The results are ordered by the
      value of the \code{by} parameter.
    }

    \item{matmirnasInMultiplePremirnas}{
      Returns mature miRNAs which are encoded in more than one
      pre-miRNA. The return object is the same than for
      \code{matrmirnas}.
    }

    \item{premirnas}{
      Returns pre-miRNAs defined by the miRBase along with optional
      additional columns from other database tables (which can be
      \code{NA} for some pre-miRNAs). Note that column \code{"sequence"}
      returns the actual RNA sequence of the pre-miRNA, not the
      genomic DNA defined by the columns \code{"pre_mirna_seq_start"}
      and \code{"pre_mirna_seq_end"}.
      Also, some pre-miRNAs might have multiple genomic alignments
      and might thus be listed multiple times in the returned object.

      The method returns a \code{DataFrame}, \code{data.frame} or
      \code{GRanges} depending on the value of the \code{return.type}
      argument (\code{"DataFrame"}, \code{"data.frame"} or
      \code{"GRanges"}, respectively). Entries in the returned object
      are ordered according to the parameter \code{order.by}, NOT by any
      ordering of values in eventually submitted filter objects.
    }

    \item{premirnasBy}{
      Returns a \code{CompressedSplitDataFrameList} of \code{DataFrame}s
      or a \code{list} of \code{data.frame}s with the names of the list
      being the ids by which the pre-miRNAs are grouped (e.g. mature miRNA
      ids) and the elements of the list being the host gene entries.

      The method returns a \code{SplitDataFrameList} (list of
      \code{DataFrame}s), a \code{list} of \code{data.frame}s or a
      \code{GRangesList}, depending on the value of the parameter
      \code{return.type} (\code{"DataFrame"}, \code{"data.frame"} or
      \code{"GRanges"}, respectively). The results are ordered by the
      value of the \code{by} parameter.
    }

    \item{premirnasWithMultipleAlignments}{
      Returns pre-miRNAs wich are encoded in several genomic loci. The
      return object is the same than for \code{premirnas}.
    }

  }
}

\section{Retrieving host genes and transcripts}{
  These methods allow to retrieve host genes and transcripts as well as
  microarray features (probe sets) targeting these.
  \describe{
    \item{hostgenes}{
      Returns all predicted host genes from the database along with
      optional additional columns from other database tables.
      Host genes with \code{gene_biotype} equal to \code{"miRNA"} should
      be taken with care, as they represent the actual
      pre-miRNAs. Ensembl defines genes for some of the pre-miRNAs
      defined in the miRBase. The column/attribute \code{database}
      specifies in which database the gene is defined (\code{"core"},
      \code{"otherfeatures"} and \code{"vega"} indicating the Ensembl
      core database with all known genes, the Ensembl otherfeatures
      database and the manually curated Ensembl vega database).

      The method returns a \code{DataFrame} or \code{data.frame}
      depending on the value of the \code{return.type} argument
      (\code{"DataFrame"} or \code{"data.frame"}). Entries in the returned
      object are ordered according to the parameter \code{order.by}, NOT
      by any ordering of values in eventually submitted filter objects.
    }

    \item{hostgenesBy}{
      Returns a \code{CompressedSplitDataFrameList} of \code{DataFrame}s
      or a \code{list} of \code{data.frame}s with the names of the list
      being the ids by which the host genes are grouped (e.g. pre-miRNA
      ids) and the elements of the list being the host gene entries.
      Note that by default empty elements are dropped (see parameter
      \code{drop.empty}).

      The method returns a \code{SplitDataFrameList} (list of
      \code{DataFrame}s) or a \code{list} of \code{data.frame}s
      depending on the value of the parameter \code{return.type}
      (\code{"DataFrame"} or \code{"data.frame"}). The results are
      ordered by the value of the \code{by} parameter.
    }

    \item{hosttx}{
      Returns all predicted host transcripts from the database along
      with optional additional columns from other database tables.
      Note that for host transcripts being the host for several
      pre-miRNAs multiple rows are present in the result table (one for
      each pre-miRNA). To get a unique list of host transcripts, the
      \code{columns} parameter should be restricted to \code{c("tx_id",
	"tx_biotype", "gene_id")}.
      The columns \code{in_intron} and \code{in_exon} specify in which
      intron or exon of the transcript the pre-miRNA is encoded (0 for
      not in intron or exon), \code{exon_id} indicates the exon id for
      exonic pre-miRNAs and the column \code{is_outside} indicates
      whether the pre-miRNA is only partially inside the transcript.
      See the package's vignette for a detailed description.

      The method returns a \code{DataFrame} or \code{data.frame}
      depending on the value of the \code{return.type} argument
      (\code{"DataFrame"} or \code{"data.frame"}). Entries in the returned
      object are ordered according to the parameter \code{order.by}, NOT
      by any ordering of values in eventually submitted filter objects.
    }

    \item{hosttxBy}{
      Returns a \code{CompressedSplitDataFrameList} of \code{DataFrame}s
      or a \code{list} of \code{data.frame}s with the names of the list
      being the ids by which the host transcripts are grouped (e.g. pre-miRNA
      ids) and the elements of the list being the host gene entries.
      Note that by default empty elements are dropped (see parameter
      \code{drop.empty}).

      The method returns a \code{SplitDataFrameList} (list of
      \code{DataFrame}s) or a \code{list} of \code{data.frame}s
      depending on the value of the parameter \code{return.type}
      (\code{"DataFrame"} or \code{"data.frame"}). The results are
      ordered by the value of the \code{by} parameter.
    }

    \item{probesets}{
      Returns microarray probe sets which where found to target the host
      transcripts. Note that in the database probe sets for different
      microarrays can be stored, thus it might be advisable to use a
      \code{\link{ArrayFilter}} to restrict to probe sets for one
      specific microarray (use \code{\link{listArrays}} to get an
      overview of all microarrays for which probe sets are available).

      The method returns a \code{DataFrame} or \code{data.frame}
      depending on the value of the \code{return.type} argument
      (\code{"DataFrame"} or \code{"data.frame"}). Entries in the returned
      object are ordered according to the parameter \code{order.by}, NOT
      by any ordering of values in eventually submitted filter objects.
    }

    \item{probesetsBy}{
      Returns microarray probe sets grouped by the column specified with
      the argument \code{by}.

      The method returns a \code{SplitDataFrameList} (list of
      \code{DataFrame}s) or a \code{list} of \code{data.frame}s
      depending on the value of the parameter \code{return.type}
      (\code{"DataFrame"} or \code{"data.frame"}). The results are
      ordered by the value of the \code{by} parameter.
    }
  }
}
\note{
  The default grouping of transcripts or genes for \code{hosttxBy} and
  \code{hostgenesBy} is by the \code{pre_mirna_algn} (i.e. the alignment
  ID of the pre-miRNA), since pre-miRNAs might have multiple genomic
  alignments and the thus returned, grouped, transcripts or genes
  might be encoded on different chromosomes.

  For the \code{matmirnas},\code{premirnas}, \code{hostgenes} and
  \code{hosttx} methods the internal SQL call uses a left join starting
  from the respective table (e.g. \code{"mature_mirna"} for
  \code{matmirnas}), thus returning all entries from that table, but
  eventually \code{NA}s for columns from other tables if no value from
  that table is linked to any of the entries in the first table.
  As a result, a call to \code{premirnas} with columns set to
  \code{"pre_mirna_name"} and \code{"tx_id"} will return the IDs of all
  pre-miRNAs and the ID of their respective putative host transcripts,
  or \code{NA} if none was defined. A call to \code{hosttx} with the
  same columns will however return less results from the database, as
  IDs of pre-miRNAs without a specified host transcripts are not
  returned (see example below).

  %% Along that lines, \code{matmirnasBy}, \code{premirnasBy}, \code{hostgenesBy} and
  %% \code{hosttxBy} will not return entries (e.g. pre-miRNAs)
  %% for which the column specified by \code{by} is empty (e.g. pre-miRNAs
  %% for which no host transcript was identified when was set to
  %% \code{by="host_tx"}). The same holds true for host genes or
  %% transcripts for which no microarray probe sets have been indentified
  %% (i.e. \code{hostgenesBy(..., by="probeset")} or
  %% \code{hosttxBy(..., by="probeset")}).
  %% However, it is still possible to
  %% return such entries by setting \code{ifEmptyBy} different from
  %% \code{NA}. In that case, all entries for which column \code{by} is
  %% empty will be combined into a list element named according to
  %% \code{ifEmptyBy} (see examples below).

  In functions \code{matmirnasBy}, \code{premirnasBy},
  \code{hostgenesBy} and \code{hosttxBy}, the internal left join starts
  from the database table in which the attribute (column) specified with
  the \code{by} argument is defined. As a consequence, entries for which
  the column specified by \code{by} is empty are NOT returned.
  To get all entries from the database, the methods \code{matmirnas},
  \code{premirnas}, \code{hostgenes} and \code{hosttx} can be used
  instead, adding additional column names to the \code{columns}
  argument.

}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{MirhostDb}}, \code{\link{listColumns}}, \code{\link{listTables}} \code{\link{makeMirhostgenesPackage}}, \code{\linkS4class{PositionFilter}}
}
\examples{



library(MirhostDb.Hsapiens.v75.v20)

## define a "shortcut" to the database
Mhdb <- MirhostDb.Hsapiens.v75.v20

##***************************************
##
##  mature miRNAs
##
##***************************************

## Simply get all mature miRNAs; the result is however not a unique list of miRNAs,
## since miRNAs from pre-miRNAs with multiple genomic alignments are listed in
## mulitple rows.
MatMir <- matmirnas(Mhdb)
MatMir
length(unique(MatMir$mat_mirna_id))

## Get mat_mirna and pre_mirna entries for mature miRNA MIMAT0000062.
MatMir <- matmirnas(Mhdb,
                    columns=unique(c(listColumns(Mhdb, "mat_mirna"),
                        listColumns(Mhdb, "pre_mirna"))),
                    filter=list(MatMirnaIdFilter("MIMAT0000062")))
MatMir
## The same mature miRNA is encoded in 3 different pre-miRNAs.

## Get all mature miRNAs along with their pre-miRNAs in which they are encoded
## and their sequence.
MatMir <- matmirnas(Mhdb, columns=c("mat_mirna_id", "mat_mirna_name",
                              "pre_mirna_name", "seq_name", "sequence"))
MatMir
length(unique(MatMir$mat_mirna_id))
length(unique(MatMir$pre_mirna_name))

## Get all mature miRNAs along with the potential host gene in which they are encoded.
MatMir <- matmirnas(Mhdb, columns=c("mat_mirna_id", "mat_mirna_name",
                              "seq_name", "gene_id", "gene_name", "gene_biotype"))
MatMir
## The mature miRNAs present in host genes.
MatMir.inhg <- MatMir[ !is.na(MatMir$gene_id), ]
MatMir.nohg <- MatMir[ is.na(MatMir$gene_id), ]

MatMir.inhg
## However, a considerable number of "host genes" are actually the pre-miRNAs, which some of them
## are stored in the Ensembl database as "gene" with the biotype "miRNA".
table(MatMir.inhg$gene_biotype)

## Now, get all mature miRNAs for which the gene_biotype!=miRNA.
MatMir <- matmirnas(Mhdb, columns=c("mat_mirna_id", "mat_mirna_name",
                              "seq_name", "gene_id", "gene_name", "gene_biotype"),
                    filter=list(GeneBiotypeFilter("miRNA", condition="!=")))
MatMir
sum(is.na(MatMir$gene_biotype))
table(MatMir$gene_biotype)

## Get all mature miRNAs as GRanges.
matmirnas(Mhdb, return.type="GRanges")

## Get all mature miRNAs that are encoded in more than one pre-miRNA.
matmirnasInMultiplePremirnas(Mhdb)

##***************************
## matmirnasBy
## Get all mature miRNAs grouped by pre-miRNA.
matmirnasBy(Mhdb, by="pre_mirna")

## Get all mature miRNAs groped by mirfam as GRanges.
matmirnasBy(Mhdb, by="mirfam", return.type="GRanges")

## Get mature miRNAs for pre-miRNA miR-16-1 and miR-16-2.
matmirnasBy(Mhdb,
            filter=list(PreMirnaFilter(c("hsa-mir-16-2", "hsa-mir-16-1"))))



##***************************************
##
##  pre-miRNAs
##
##***************************************

## Get all pre-miRNAs.
PreMir <- premirnas(Mhdb)
PreMir
length(unique(PreMir$pre_mirna_name))

## Get all pre-miRNAs as GRanges.
premirnas(Mhdb, return.type="GRanges")

## Get all pre-miRNAs along with their miRNA family and their sequence.
## Since we don't ask for the pre_mirna_seq_start and end we get a
## unique table of pre-miRNAs.
PreMir <- premirnas(Mhdb, columns=c("pre_mirna_name", "mirfam_name",
                              "sequence"))
PreMir

## We have some pre-miRNAs without family
sum(is.na(PreMir$mirfam_name))
## but none without sequence.
sum(is.na(PreMir$sequence))

## Get all pre-miRNAs with multiple genomic alignments.
premirnasWithMultipleAlignments(Mhdb)

##***************************
## premirnasBy
## Get the pre-miRNAs by the mature_mirna.
PB <- premirnasBy(Mhdb, by="mat_mirna")

## Add also additional stuff and fetch all pre-miRNAs for host gene SMC4:
premirnasBy(Mhdb, columns=c("pre_mirna_name", "sequence", "mirfam_name",
                      "mat_mirna_name"), filter=list(GenenameFilter("SMC4")))

## Get all pre-miRNAs by host_gene SMC4.
premirnasBy(Mhdb, by="host_gene", filter=list(GenenameFilter("SMC4")))


## Get all pre-miRNAs by host_gene SMC4 as GRanges.
premirnasBy(Mhdb, by="host_gene", filter=list(GenenameFilter("SMC4")),
            return.type="GRanges")


##***************************************
##
##  host transcripts
##
##***************************************

## Get all host transcripts from the database.
HT <- hosttx(Mhdb)
HT
nrow(HT)
## The same host_tx might be the host for multiple miRNAs, thus we do have non-unique tx_ids.
length(unique(HT$tx_id))

## Get a unique table of host transcripts.
HT <- hosttx(Mhdb, columns=c("tx_id", "tx_biotype", "gene_id"))
HT
nrow(HT)
length(unique(HT$tx_id))

## Get the host transcripts along with the corresponding gene.
HT <- hosttx(Mhdb, columns=c("tx_id", "in_intron", "in_exon", "gene_id",
                       "gene_name", "entrezid", "database"))
HT
## In what databases are these transcripts defined?
table(HT$database)
nrow(HT)

## Note that the information from the various databases is redundant
## (e.g. the same gene can be defined in the Ensembl code database as
## well as in the NCBI RefSeq database which genes are provided through
## the Ensembl otherfeatures database.
## To avoid getting redundant entries it is possible to use a
## DatabaseFilter:
HT <- hosttx(Mhdb, columns=c("tx_id", "in_intron", "in_exon", "gene_id",
                       "gene_name", "entrezid", "database"),
             filter=list(DatabaseFilter("core")))
HT
nrow(HT)



## Include now also the pre_mirna ids.
HT <- hosttx(Mhdb, columns=c("tx_id", "in_intron", "in_exon", "gene_id",
                       "gene_name", "entrezid", "database",
                       "pre_mirna_id", "pre_mirna_name"))
HT
nrow(HT)
## We have now more rows, since different pre-miRNAs might be
## associated with the same host_tx.
length(unique(HT$tx_id))


##***************************
## hosttxBy
## Get the host transcripts by the pre-miRNA
## this will drop automatically empty entries, i.e. pre-miRNAs for which
## no host transcript was defined.
HT <- hosttxBy(Mhdb, by="pre_mirna", columns=c("tx_id", "tx_biotype",
                                         "in_intron", "in_exon",
                                         "pre_mirna_name"))
HT

## To get all of them we scan set drop.empty=FALSE.
HT <- hosttxBy(Mhdb, by="pre_mirna",
               columns=c("tx_id", "tx_biotype", "in_intron", "in_exon",
                   "pre_mirna_name"), drop.empty=FALSE)
HT

## There are however also some without any entries:
empties <- unlist(lapply(HT, function(z){ return(all(is.na(z$tx_id))) }))
sum(empties)
HT[ empties ]

## Host transcripts by gene.
HT <- hosttxBy(Mhdb, by="host_gene")
HT



##***************************************
##
##  host genes
##
##***************************************

## With the host genes it is just the same as above.
HG <- hostgenes(Mhdb)
HG
length(unique(HG$gene_id))
nrow(HG)


##***************************
## hostgenesBy
## Get the host genes by the pre-miRNA.
HG <- hostgenesBy(Mhdb, by="pre_mirna")
HG

## Get host genes by mirfam.
HG <- hostgenesBy(Mhdb, by="mirfam",
                  columns=c("gene_id", "gene_name", "mirfam_name"))
HG



##***************************************
##
##  probe sets
##
##***************************************

## First get a list of microarrays for which probe sets are available.
listArrays(Mhdb)

AF <- ArrayFilter("HG-U133_Plus_2")

## Get all probe sets from the database along with the gene name and
## the pre-miRNA name.
PS <- probesets(Mhdb, columns=c(listColumns(Mhdb, "array_feature" ),
                          "gene_name", "pre_mirna_name"), filter=list(AF))
PS

## Get all probe sets grouped by pre-miRNA name.
PS <- probesetsBy(Mhdb, by="pre_mirna", use.names=TRUE, filter=list(AF))
PS


##***************************************
##
##  The effect of the left join
##
##***************************************
## Get all pre-miRNAs and the ID of the host transcript.
fromPre <- premirnas(Mhdb, columns=c("pre_mirna_name", "tx_id"))
## Get the same columns, but starting from table "host_tx"
fromTx <- hosttx(Mhdb, columns=c("pre_mirna_name", "tx_id"))
## We have less rows for the latter query.
nrow(fromPre)
nrow(fromTx)

## The reason being, that pre-miRNAs without host transcript are not returned
## by the second query, while they are for the first.
sum(is.na(fromPre$tx_id))
sum(is.na(fromTx$tx_id))

}
\keyword{classes}

