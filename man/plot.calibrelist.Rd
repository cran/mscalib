\name{plot.calibrelist}
\alias{plot.calibrelist}
\title{ Plot}
\description{ Shows a matrix of scatterplots for different varialbes.. }
\usage{\method{plot}{calibrelist}(x,...)}
\arguments{
\item{...}{ graphical parameters can be given as arguments to plot.}
\item{x}{ calibrelist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 res2 <- getrecalib(mvl) # get recalibration model for not filtered data
 plot(res2)
}
\keyword{misc}
