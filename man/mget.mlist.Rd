\name{mget.mlist}
\alias{mget.mlist}
\title{ Field Access}
\description{ Acces fields in object of class myobj }
\usage{\method{mget}{mlist}(object,attrn,...)}
\arguments{
\item{...}{}
\item{attrn}{ name of field (Attribute)}
\item{object}{ object of class  mlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 mget(mvl)
 mget(mvl,"info")
 mget(mvl,"tcoor")
}
\keyword{misc}
