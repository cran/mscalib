\name{as.matrix.massvectorlist}
\alias{as.matrix.massvectorlist}
\title{ Matrices}
\description{ Turns the massvectorlist into a matrix. On each massvectorlist the \code{summary} method is called.  The result of \code{summary} are combined into a matrix. }
\usage{\method{as.matrix}{massvectorlist}(x)}
\arguments{
\item{x}{ massvectorlist}
}
\value{
\item{ matrix }{ matrix with columns }
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.data.frame.massvectorlist}},
\code{\link[base]{as.matrix}},
\code{link[base]{as.data.frame}},
}
\examples{
 data(mvl)
 tmp<-as.matrix(mvl)
 print(colnames(mvl))
 plot(tmp[,"lengthmv"],tmp[,"mass.Min."])
}
\keyword{misc}
