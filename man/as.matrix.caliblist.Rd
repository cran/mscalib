\name{as.matrix.caliblist}
\alias{as.matrix.caliblist}
\title{ Matrices}
\description{ Turns the caliblist into a matrix }
\usage{\method{as.matrix}{caliblist}(x)}
\arguments{
\item{x}{ caliblist}
}
\value{
\item{ matrix }{ matrix with masses and intensities.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.data.frame.caliblist}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 test<-getintcalib(mvl,cal,error=500)
 tmp<-as.matrix(test)
 colnames(tmp)
 dim(tmp)
}
\keyword{misc}
