\name{plot.calibintlist}
\alias{plot.calibintlist}
\title{ Calibintlist Plotting}
\description{ A matrix of scatterplots is produced. }
\usage{\method{plot}{calibintlist}(x,...)}
\arguments{
\item{...}{ graphical parameters can be given as arguments to plot.}
\item{x}{ object of class calibintlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 data(cal)
 mvl <- mvl[1:100]
 ires <- getintcalib(mvl,cal,error=250)
 plot(ires)
}
\keyword{misc}
