\name{as.vector.massvector}
\alias{as.vector.massvector}
\title{ Vector}
\description{ Casts massvector into vector. It does it by calling the summary.massvector method first and unlisting the result. }
\usage{\method{as.vector}{massvector}(x, mode="any")}
\arguments{
\item{mode}{}
\item{x}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 as.vector(mv1)
}
\keyword{misc}
