\name{hist.massvectorlist}
\alias{hist.massvectorlist}
\title{ Histograms}
\description{ Histogram of mass frequencies in the massvectorlist. }
\usage{\method{hist}{massvectorlist}(x,accur = 0.1, main=info(mvl) ,xlab="m/z",xlim=c(700,4500),add=FALSE,col=1,...)}
\arguments{
\item{...}{ further parameters.}
\item{accur}{ bin width for plotting mass frequencies.}
\item{add}{ logical; If \code{TRUE} add to already existing plot.}
\item{col}{ color of the histogram,}
\item{main}{ title of graph}
\item{x}{ massvectorlist}
\item{xlab}{ a title for the x axis.}
\item{xlim}{ the range to be encompassed by the x axis.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{hist.massvector}},
}
\examples{
 data(mvl)
 hist(mvl)
}
\keyword{misc}
