\name{c.massvector}
\alias{c.massvector}
\title{ Combine Massvectors into one Massvector.}
\description{ Combines Massvectors into one Massvector. }
\usage{\method{c}{massvector}(x,...)}
\arguments{
\item{...}{ massvectors to be concatenated.}
\item{x}{ massvector}
}
\value{
\item{ massvector }{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{rbind}},
}
\examples{
 data(mv1)
 data(mv2)
 par(mfrow=c(2,1))
 plot(mv1)
 plot(mv2,add=TRUE)
 plot(c(mv1,mv2))
}
\keyword{misc}
