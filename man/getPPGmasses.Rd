\name{getPPGmasses}
\alias{getPPGmasses}
\title{ PPG masses}
\description{ Computes poly-(propylene glycol) masses. }
\usage{getPPGmasses(start=10,end=100)}
\arguments{
\item{end}{ length of longest polymer.}
\item{start}{ length of shortest polymer.}
}
\value{
\item{ massvector }{ massvector of ppg masses}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 plot(getPPGmasses(start=12,end=100))
}
\keyword{misc}
