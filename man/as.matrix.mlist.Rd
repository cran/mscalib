\name{as.matrix.mlist}
\alias{as.matrix.mlist}
\title{ Matrices}
\description{ Turns the caliblist into a matrix }
\usage{\method{as.matrix}{mlist}(x)}
\arguments{
\item{x}{ caliblist}
}
\value{
\item{ matrix }{ matrix}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.data.frame.mlist}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 test<-getintcalib(mvl,cal,error=500)
 tmp<-as.matrix(test)
 colnames(tmp)
 dim(tmp)
 data(mvl)
 tmp<-as.matrix(mvl)
 print(colnames(mvl))
 plot(tmp[,"lengthmv"],tmp[,"mass.Min."])
}
\keyword{misc}
