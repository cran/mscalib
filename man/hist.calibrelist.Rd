\name{hist.calibrelist}
\alias{hist.calibrelist}
\title{ Histogram Plot}
\description{ Computes Histograms of the lengthmv, Coef.Intercept, Coef.Slope, PQM. }
\usage{\method{hist}{calibrelist}(x,...)}
\arguments{
\item{...}{ further graphical parameters.}
\item{x}{ calibrelist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 res2 <- getrecalib(mvl) # get recalibration model for not filtered data
 hist(res2)
}
\keyword{misc}
