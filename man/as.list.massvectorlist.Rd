\name{as.list.massvectorlist}
\alias{as.list.massvectorlist}
\title{ List}
\description{ Turns the massvectorlist into a list }
\usage{\method{as.list}{massvectorlist}(x,...)}
\arguments{
\item{...}{ further parameter.}
\item{x}{ massvectorlist.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 res<-as.list(mvl)
 class(res)
}
\keyword{misc}
