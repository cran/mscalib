\name{plot.calibintstat}
\alias{plot.calibintstat}
\title{ Calibinstat Plotting}
\description{ Plots a line with intercept and slope given by the error model. }
\usage{\method{plot}{calibintstat}(x,...)}
\arguments{
\item{...}{ graphical parameters can be given as arguments to plot.}
\item{x}{ object of class calibintstat}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
 plot(res)
}
\keyword{misc}
