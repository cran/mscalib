\name{mvFilter.massvector}
\alias{mvFilter.massvector}
\title{ Filtering Massvector}
\description{ Filters massvector for masses given in a second massvector. }
\usage{\method{mvFilter}{massvector}(object,fmass,match=FALSE,error=250,ppm=TRUE,uniq=FALSE,...)}
\arguments{
\item{...}{}
\item{error}{ mesurment error.}
\item{fmass}{}
\item{match}{ logical; \code{TRUE} - than returns masses matching to the masses in massvector abundant, \code{FALSE} - returns masses not matching.}
\item{object}{ massvector}
\item{ppm}{ given either in ppm (\code{TRUE}) or as absolut error (F).}
\item{uniq}{ logical; \code{FALSE} - returns all masses in the range given by error. \code{TRUE} - returns only the closest mass.}
}
\value{
\item{ massvector }{ with matching or not matchin masses.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(mv2)
 mvFilter(mv1,mv2,error=250,match=FALSE)
 mvFilter(mv1,mv2,error=250,match=TRUE)
}
\keyword{misc}
