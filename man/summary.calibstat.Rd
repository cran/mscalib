\name{summary.calibstat}
\alias{summary.calibstat}
\title{ Calibstat summary}
\description{ Generates Summary. }
\usage{\method{summary}{calibstat}(object,...)}
\arguments{
\item{...}{}
\item{object}{ object of class calibstat.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(cal)
 test<-getintcalib(mv1,cal)
 summary(test)
}
\keyword{misc}
