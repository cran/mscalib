\name{mvFilter.massvectorlist}
\alias{mvFilter.massvectorlist}
\title{ Filtering Massvector}
\description{ Filters massvector for masses given in a second massvector. }
\usage{\method{mvFilter}{massvectorlist}(object,fmass,match=FALSE,error=250,ppm=TRUE,uniq=FALSE,...)}
\arguments{
\item{...}{ further parameters.}
\item{error}{ mesurment error.}
\item{fmass}{}
\item{match}{ logical;\code{TRUE} - than returns masses matching to the masses in massvector abundant, \code{FALSE} - returns masses not matching.}
\item{object}{ massvectorlist}
\item{ppm}{ given either in ppm (\code{TRUE}) or as absolut error (\code{FALSE}).}
\item{uniq}{ logical; \code{FALSE} - returns or removes all masses in the range given by error. \code{TRUE} - returns ore removes only the closest mass.}
}
\value{
\item{ massvectorlist }{ with massvectors with matching or not matching masses.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{mvFilter.massvector}},
}
\examples{
 data(mvl)
 mvl<-mvFilter(mvl,mvl[[1]],match=FALSE,error=250)
 length(mvl[[1]])
}
\keyword{misc}
