\name{setParms<-.myobj}
\alias{setParms<-.myobj}
\title{ Field Access}
\description{ Set attributes in object of class myobj }
\usage{}
\arguments{
\item{object}{ object of class myobj}
\item{value}{ a list where list names are attributes names.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 mv1
 setParms(mv1) <- list(info="test")
 mv1
}
\keyword{misc}
