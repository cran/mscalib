\name{getdiff.massvector}
\alias{getdiff.massvector}
\title{ Massdifferences}
\description{ Computes mass differences in the object. }
\usage{\method{getdiff}{massvector}(object,range=c(0,100),...)}
\arguments{
\item{...}{}
\item{object}{ massvector}
\item{range}{}
}
\details{ Removes one of the masses contributing to a mass difference given in the list of diffs.  Can be used if a variable modification are present in the massvector but can not be considered by the identification software. }
\value{
\item{ massvector }{ filtered massvector.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\keyword{misc}
