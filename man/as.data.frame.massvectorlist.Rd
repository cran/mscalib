\name{as.data.frame.massvectorlist}
\alias{as.data.frame.massvectorlist}
\title{ Data Frames}
\description{ Turns the massvectorlist into a data.frame }
\usage{\method{as.data.frame}{massvectorlist}(x,...)}
\arguments{
\item{...}{ further arguments}
\item{x}{ object of class massvectorlist}
}
\details{ These functions create a data frame, tightly coupled  collections of variables which share many of the properties of  matrices and of lists, used as the fundamental data structure by  most of R's modeling software. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.matrix.massvectorlist}},
\code{\link{as.matrix}},
}
\examples{
 data(mvl)
 tmp <- as.data.frame(mvl)
 names(mvl)
 plot(tmp$lengthmv,tmp$mass.Min.)
}
\keyword{misc}
