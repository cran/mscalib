\name{setParms<-.mylistobj}
\alias{setParms<-.mylistobj}
\title{ Field Access}
\description{ Set attributes in object of class mylistobj }
\usage{}
\arguments{
\item{object}{ object of class myobj}
\item{value}{ a list where list names are attributes names.}
}
\details{ mylistobj provides basic functionality for objects implemented using list to store the attributes. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{setParms<-.myobj}},
\code{\link{mget.mylistobj}},
}
\examples{
 data(mv1)
 res<-getrecalib(mv1)
 setParms(res) <- list(info="test")
 res
}
\keyword{misc}
