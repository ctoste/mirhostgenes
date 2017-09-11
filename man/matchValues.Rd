\name{tranferValues}
\Rdversion{1.1}
\alias{transferValues}
\alias{chooseAll}
\alias{chooseClosestForEach}
\alias{chooseOrderedValue}
\alias{doSelectData}
\alias{transferValues,numeric,MirhostDb-method}
\alias{doPairData}

\title{Pairing host gene and mature miRNA expression data}
\description{
  These functions allow to pair/match values by name between any two
  named vectors. The \code{transferValues} method allows to map/transfer
  values from one type of id (e.g. mature miRNA name) to any other id
  defined in the \code{MirhostDb} database. Combined, these functions
  enable to map between expression/regulation values of miRNAs and their
  corresponding host genes.
}
\usage{

chooseAll(x, y)

chooseClosestForEach(x, y, abs=TRUE)

chooseOrderedValue(x, y, orderFun=function(z){order(z, decreasing=TRUE)})

doPairData(x, y, chooseFunX=chooseAll, chooseFunY=chooseAll)

doSelectData(x, chooseFunX=chooseAll)

\S4method{transferValues}{numeric,MirhostDb}(x, object,
                          xNamesAre="mat_mirna_name",
                          toNames="pre_mirna_name",
                          solveFun=chooseAll,
                          filter=list(),
                          na.rm=FALSE)


}
\arguments{
  (in alphabetical order)

  \item{abs}{
    For \code{chooseClosestBetween}: if the values should be ordered by
    the absolute value of their differences.
  }

  \item{chooseFunX}{
    The function applied to input vector \code{x} to select a (or all)
    value(s) among all values with the same name. See description below
    for more details.
  }
  \item{chooseFunY}{
    The function to be applied to input vector \code{y}.
  }

  \item{filter}{
    For \code{transferValues}: filters to be used in the query to
    retrieve the mapping between ids. To map values from miRNAs to host
    genes or transcripts it is advisable to use a \code{DatabaseFilter}
    to restrict to genes defined in one of the databases (\code{"core"},
    \code{"vega"} or \code{"otherfeatures"}), otherwise the values are
    mapped to a redundant list of genes (different ids for the same
    genes defined in the different databases).
  }

  \item{na.rm}{
    For \code{transferValues}: for \code{na.rm=FALSE} (default) the
    function returns also input values for which no mapping between
    ids can be established. Otherwise, only values for which the names
    can be mapped to ids of the type \code{toNames} are returned.
  }

  \item{object}{
    For \code{transferValues}: the \code{MirhostDb} object.
  }

  \item{orderFun}{
    For \code{chooseOrderedValue}: the function how values should be
    ordered.
  }

  \item{solveFun}{
    For \code{transferValues}: the function that should be applied to
    resolve multi-mappings between ids. See description and examples
    below for more details.
  }

  \item{toNames}{
    Defines to which ids the values should be transfered/mapped. Allowed values
    are \code{"pre_mirna_id"}, \code{"pre_mirna_name"}, \code{"tx_id"},
    \code{"mat_mirna_id"}, \code{"mat_mirna_name"}, \code{"gene_id"},
    \code{"gene_name"} or \code{"probeset_id"}.
  }

  \item{x}{
    A named numerical vector.
  }

  \item{xNamesAre}{
    Defines which type of ids the names of \code{x} are. See
    \code{toNames} for supported type of ids.
  }

  \item{y}{
    A named numerical vector.
  }

}
\section{Simple utilities to select values in a numeric vector}{
  These functions allow to select specific, or all, values for
  same-named entities in a numeric vector.

  \describe{
    \item{chooseAll}{
      This function returns a numeric sequence from 1 to
      \code{length(x)}.
    }
    \item{chooseOrderedValue}{
      This function returns (by default) the index of the largest value
      in \code{x}.
    }
    \item{chooseClosestForEach}{
      This function returns for each value in \code{y}, the index of the
      closest value in \code{x}. The closest value is defined as the value with
      the smallest abolute difference between values. The length of the
      returned numeric vector corresponds to the length of \code{y}.
    }
    \item{doSelectData}{
      Selects values for same named entities in a submitted numeric
      vector. The function uses one of the \emph{choose} functions
      defined above to select the value and returns a \code{data.frame}
      containing the selected values in column \code{"x"}, their index
      in the original input argument \code{x} in column \code{"x.idx"}
      and their corresponding name in column \code{"name"}.
    }
  }
}

\section{Pairing data}{

  \describe{
    \item{doPairData}{
      This function takes two named numeric vectors and matches/pairs values
      between them based on their names.
      The length of the two vectors do not have to match, names not
      present in one or the other vector are added.
      The aim of the function is to pair entries in input vector \code{x} with
      values in input vector \code{y} based on their names, i.e. values
      with the same name are matched with each other and returned in the
      same row in the result table. If \code{x} and/or \code{y} each
      contains entries with the same name, values are repeated such that
      each value in \code{x} is paired with each value in \code{y} with the
      same name (i.e. if \code{x} contains 2 values with the name
      \code{"a"}, and \code{y} 3 with that name, the function pairs
      each of the 2 values in \code{x} with each of the 3 values in
      \code{y} resulting in 6 rows in the result table for name
      \code{"a"}). This default behaviour can be changed by specifying a
      function different than \code{chooseAll} for argument \code{chooseFunX}
      or \code{chooseFunY}. Note that \code{chooseFunX} is applied
      first, before \code{chooseFunY}, thus resulting eventually in
      different results when arguments \code{x} and \code{y} are
      swapped. See examples and note below for more information.


      The function returns a \code{data.frame} with columns
      \code{"name"} (the common name for the values in \code{x} and
      \code{y}), \code{"x"} (the values from input vector \code{x}),
      \code{"x.idx"} (an integer value representing the original index of
      the value in the input vector \code{x}), \code{"y"} (the values from
      the input vector \code{y}, paired to x) and \code{"y.idx"} (the
      index of the value in the input vector \code{y}).
    }


    \item{transferValues}{
      This function takes a named numeric vector as input and maps the
      names to the new type of ids specified by \code{toNames}. Only
      mappings between ids defined in the \code{MirhostDb} are supported
      (i.e. \code{"pre_mirna_id"}, \code{"pre_mirna_name"}, \code{"tx_id"},
      \code{"mat_mirna_id"}, \code{"mat_mirna_name"}, \code{"gene_id"},
      \code{"gene_name"} or \code{"probeset_id"}).

      Multi-mappings between ids (e.g. \code{"mat_mirna_name"} with
      \code{"pre_mirna_name"} as each pre-miRNA encodes two mature
      miRNAs) can be resolved with \code{solveFun}. The same functions
      than for \code{chooseFunX} of the \code{pairData} are supported,
      i.e. \code{chooseAll} (selects all ids eventually repeating the
      value) or \code{chooseOrderedValue} (selects the value of the first
      id ordered by the sorting function). See below for examples.

      The method returns a \code{data.frame} with 3 columns, the first
      (\code{"x"}) with the numeric values from the input argument
      \code{x}, the second (named according to argument \code{toName})
      with the ids to which the values have been mapped and the third
      (named according to argument \code{xNamesAre}) with the original
      names in the input argument \code{x}.
    }
  }
}
\note{
      The \code{chooseFunX} and \code{chooseFunY} can be used to modify
      the way values are selected for same named entries within the
      input vector \code{x} or \code{y}. The \code{chooseAll} function
      simply selects all values. \code{chooseOrderedValue} chooses the
      first value of all values with the same name ordered according to
      the specified scheme (selecting by default the largest
      value). \code{chooseClosestBetween} chooses the one value in its
      input argument \code{x} which is closest (most similar value) to
      any of the values in its input argument \code{y}.

      Any function returning indices can be defined and used with
      arguments \code{chooseFunX} and \code{chooseFunY}.

      Possible scenarios for pairing of same-name values:
      \describe{
	\item{Pair the largest values in each vector}{
	  Specify \code{chooseOrderedValue} for both, \code{chooseFunX}
	  and \code{chooseFunY}.
	}
	\item{Get the largest value in \code{x} and pair that with the
	  most similar in \code{y}}{
	  Specify \code{chooseOrderedValue} for \code{chooseFunX} and
	  \code{chooseClosestForEach} for \code{chooseFunY}.
	}
      }
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{MirhostDb}},
  \code{\link{listColumns}},
  \code{\link{hostgenes}},
  \code{\link{hosttx}}
}
\examples{

########################
##
##  Simple value pairing between two numeric vectors.
##
## A and B are named numeric vectors, names corresponding e.g. to
## gene names and values to their expression. Each name can occur
## several times.
##
A <- 1:10
names(A) <- c("b", "a", "c", "b", "a", "a", "e", "f", "a", "c")

B <- c(2, 5, 3, 2, 5, 16, 4, 20)
names(B) <- c("a", "b", "a", "c", "c", "c", "n", "f")

## In the most generic case we want to pair each value in A with each
## value in B with the same name. Each value in A will be repeated as
## many times as there are values with the same name in B and vice
## versa.
doPairData(A, B)

## Next we want to pair all values in A to the value in B with the
## largest value (for the same name). This returns each value in A
## only once.
doPairData(A, B, chooseFunY=chooseOrderedValue)

## Next we select the largest value for all same-named entries for
## both input vectors.
doPairData(A, B, chooseFunX=chooseOrderedValue,
           chooseFunY=chooseOrderedValue)

## Next we first select the largest among all values with the same name
## in x and select from y the value that is closest to that (most similar)
doPairData(A, B, chooseFunX=chooseOrderedValue,
           chooseFunY=chooseClosestForEach)

## Next we reverse the order, i.e. choose first for each value in B the closest
## value in A and next select the largest balue in B.
doPairData(A, B, chooseFunX=chooseClosestForEach,
           chooseFunY=chooseOrderedValue)
## As we see, the result is different.

## At last we select all pairs of most similar values.
doPairData(A, B, chooseFunX=chooseClosestForEach,
              chooseFunY=chooseClosestForEach)
## The function does the following:
## for name b: A x(1, 4), B y(5):
##             choose x: (4) (closest to 5)
##             choose y: (5) (closest to 4)
## for name a: A x(2, 5, 6, 9), B y(2, 3)
##             choose x: (2, 2)
##             choose y: (2, 2)
## for name c: A x(3, 10), B y(2, 5, 8)
##             choose x: (3, 3, 10)
##             choose y: (2, 2, 8)
## for name f: A x(8), B y(20)
##             choose x: (8)
##             choose y: (20)
## for name e: A x(7), B y(NA)
## for name n: A x(NA), B y(4)
## At last unique is called on the resulting data.frame that reduced the
## values for a and c.

## another example:
C <- c(3, 8, 7, 1)
names(C) <- rep("a", length(C))
D <- c(9, 2, 4, 20, 3)
names(D) <- rep("a", length(D))
doPairData(C, D, chooseFunX=chooseClosestForEach,
           chooseFunY=chooseClosestForEach)


##################
##
##  simple examples for the "chooseFuns"
##
chooseAll(1:5)

chooseOrderedValue(1:5)

chooseClosestForEach(1:10, 5:1)

someVals <- 1:10
names(someVals) <- c("a", "c", "a", "b", "c", "a", "a", "b", "f", "d")
## just return all of them.
doSelectData(someVals)
## select for each the largest value
doSelectData(someVals, chooseFunX=chooseOrderedValue)



##################
##
##  mapping/transfering values
##
## Map mature miRNA expression values to pre-miRNA ids.
## Each pair of mature miRNAs is encoded in a single pre-miRNA. In addition,
## miR-16-5p is encoded in two different pre-miRNAs.
library(MirhostDb.Hsapiens.v75.v20)
mhg <- MirhostDb.Hsapiens.v75.v20
miRNA.exprs <- c(8, 8.3, 5.6, 9.5, 4.6, 13.1)
names(miRNA.exprs) <- c("hsa-miR-15b-3p", "hsa-miR-15b-5p", "hsa-miR-16-5p",
                        "hsa-miR-16-1-3p", "hsa-miR-223-3p", "hsa-miR-223-5p")

transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
               toNames="pre_mirna_name")
## we can also use chooseOrderedValue as solve function to select for each
## pre-miRNA only the mature miRNA with the highest expression. That way we get
## a unique value for each pre-miRNA.
transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
               toNames="pre_mirna_name", solveFun=chooseOrderedValue)

## A "standard" filter that might be useful to map miRNA data to genes.
df <- list(DatabaseFilter("core"),
           GeneBiotypeFilter("miRNA", condition="!="),
           ArrayFilter("HG-U133_Plus_2"))


}
\keyword{classes}

