\name{hist.calibintlist}
\alias{hist.calibintlist}
\title{ Histogram Plot}
\description{ Computes histograms of the lengthmv, Coef.Intercept, Coef.Slope ... }
\usage{\method{hist}{calibintlist}(x,...)}
\arguments{
\item{...}{ further graphical parameters.}
\item{x}{ calibrelist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 data(cal)
 mvl <- mvl[1:100]
 data(cal)
 ires <- getintcalib(mvl,cal,error=250)
 hist(ires)
}
\keyword{misc}
