\name{mget.massvector}
\alias{mget.massvector}
\title{ Field Access}
\description{ Access to the fields in the massvector }
\usage{\method{mget}{massvector}(object,attrn,...)}
\arguments{
\item{...}{ further parameters}
\item{attrn}{ The value of which field to return. If missing the fields of the objects are returned.}
\item{object}{ massvector}
}
\value{
\item{ xxx }{ depends which field in the massvector are accessed.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 mget(mv1,"info")
 mget(mv1,"peaks")
}
\keyword{misc}
