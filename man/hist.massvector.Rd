\name{hist.massvector}
\alias{hist.massvector}
\title{ Histograms}
\description{ Histograms }
\usage{\method{hist}{massvector}(x,accur = 0.1,abund = 0, main=info(x) ,xlab="m/z",xlim=c(min(mass(x)),max(mass(x))),add=FALSE,col=1,...)}
\arguments{
\item{...}{ further plotting arguments.}
\item{abund}{ draws a horizontal line at the frequency given by abund.}
\item{accur}{ sets the bin width of the histogramm.}
\item{add}{ T-adds the histogram to an existing image.}
\item{col}{ the color of the histogram.}
\item{main}{}
\item{x}{}
\item{xlab}{ sets the xlabels.}
\item{xlim}{ sets the min and max value to be displayed.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{hist}},
}
\examples{
 data(mv1)
 hist(mv1)
}
\keyword{misc}
