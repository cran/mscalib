\name{mass.massvector}
\alias{mass.massvector}
\title{ Mass Access.}
\description{ Access to the mass field of a massvector. }
\usage{\method{mass}{massvector}(object,mas,...)}
\arguments{
\item{...}{}
\item{mas}{ A array with masses. If missing function returns masses.}
\item{object}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 mass(mv1)
 mass(mv1,1:10)
}
\keyword{misc}
