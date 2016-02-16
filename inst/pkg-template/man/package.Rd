\name{@PKGNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}

\title{@PKGTITLE@}

\description{
  This package loads an SQL connection to a database containing the
  definition of host genes for all miRNAs from the species
  @SPECIES@. Host genes have been defined using Ensembl databases,
  miRNAs were defined by the miRBase (versions: @PROVIDERVERSION@).
  For examples and help on functions see the help pages from the
  \code{mirhostgenes} package!

}

\note{
  This data package was made using the perl scripts from the
  \code{mirhostgenes} package.

}

\author{@AUTHOR@}

\examples{
## load the library
library(@PKGNAME@)
## list the contents that are loaded into memory
ls('package:@PKGNAME@')
## show the db object that is loaded by calling it's name
@PKGNAME@

## for more examples see the mirhostgenes package.


}

\keyword{package}
\keyword{data}
