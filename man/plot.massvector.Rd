\name{plot.massvector}
\alias{plot.massvector}
\title{ Massvector Plotting}
\description{ Function for plotting massvectors. If one massvector are given it shows a stick masspectrum.  If two massvectors are given their masses are plotted against each other. }
\usage{\method{plot}{massvector}(x,...)}
\arguments{
\item{...}{ a second massvector and graphical parameters can be given as arguments to `plot'.}
\item{x}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(mv2)
 plot(mv1)
 plot(mv1,mv2)
}
\keyword{misc}
