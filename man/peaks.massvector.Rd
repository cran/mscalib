\name{peaks.massvector}
\alias{peaks.massvector}
\title{ Data Access}
\description{ Access to the mass field of the massvector. }
\usage{\method{peaks}{massvector}(object,masses,...)}
\arguments{
\item{...}{}
\item{masses}{ matrix with masses in first column. If missing the matrix of the massvector is returned.}
\item{object}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{mass.massvector}},
}
\examples{
 mv1<-massvector()
 mv1<- peaks(mv1,cbind(1:10,1:10))
 peaks(mv1)
}
\keyword{misc}
