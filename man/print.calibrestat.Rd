\name{print.calibrestat}
\alias{print.calibrestat}
\title{ Print Calibrestat Object}
\description{ Prints the fields of the calibrestat object. }
\usage{\method{print}{calibrestat}(x,...)}
\arguments{
\item{...}{}
\item{x}{ object of class calibrestat}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{print}},
}
\examples{
 data(mv1)
 res<-getrecalib(mv1)
 print(res)
}
\keyword{misc}
