\name{plot.massvectorlist}
\alias{plot.massvectorlist}
\title{ Massvectorlist Plotting}
\description{ Plots masses (m/z) in massvector against sample in massvectorlist. }
\usage{\method{plot}{massvectorlist}(x, main=info(mvl) ,xlab="m/z",xlim=c(700,4500),add=FALSE,col=1,cex=0.5,...)}
\arguments{
\item{...}{ graphical parameters can be given as arguments to `plot'.}
\item{add}{ \code{TRUE} - masses of a new massvectorlist are added to existing plot.}
\item{cex}{}
\item{col}{ color of the points denoting the masses.}
\item{main}{ an overall title for the plot.}
\item{x}{ massvectorlist}
\item{xlab}{ a title for the x axis.}
\item{xlim}{ the range to be encompassed by the x axis.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{hist.massvectorlist}},
\code{\link{plot.massvector}},
\code{\link{image.mlist}},
\code{\link{plot}},
\code{\link{par}},
}
\examples{
 data(mvl)
 plot(mvl,col= 3 )
 tt<- gamasses(mvl,abund=50)
 plot(mvl,col=1, xlim=c(tt[1,1]-0.4,tt[1,1]+0.4))
}
\keyword{misc}
