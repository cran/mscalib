\name{as.matrix.massvector}
\alias{as.matrix.massvector}
\title{ Matrices}
\description{ Turns the massvector into a matrix }
\usage{\method{as.matrix}{massvector}(x)}
\arguments{
\item{x}{ massvector}
}
\value{
\item{ matrix }{ matrix with masses and intensities.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{peaks.massvector}},
}
\examples{
 data(mv1)
 res<-as.matrix(mv1)
 class(res)
}
\keyword{misc}
