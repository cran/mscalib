\name{mget.mylistobj}
\alias{mget.mylistobj}
\title{ Field Access}
\description{ Access to the fields in the mylistobj }
\usage{\method{mget}{mylistobj}(object,attrn,...)}
\arguments{
\item{...}{ further parameters}
\item{attrn}{ The of the field to access. If missing the fields of the objects are returned.}
\item{object}{ mylistobj}
}
\details{ mylistobj provides basic functionality for objects implemented using list to store the attributes. }
\value{
\item{ xxx }{ depends which field in the massvector are accessed.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 res <- getrecalib(mv1)
 class(res)
 mget(mv1,"lenghtmv")
 mget(mv1,"Coef.Intercept")
}
\keyword{misc}
