\name{plot.calibspline}
\alias{plot.calibspline}
\title{ Calibspline Plotting}
\description{ Function for plotting object of class calibspline. }
\usage{\method{plot}{calibspline}(x,...)}
\arguments{
\item{...}{ graphical parameters can be given as arguments to `plot'.}
\item{x}{ calibspline}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(ppg)
 tmp <- getextcalib(ppg,error=200)
 plot(tmp)
}
\keyword{misc}
