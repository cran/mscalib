\name{as.data.frame.mlist}
\alias{as.data.frame.mlist}
\title{ Data Frames}
\description{ Turns the mlist object into a data.frame }
\usage{\method{as.data.frame}{mlist}(x,row.names=NULL,optional = FALSE)}
\arguments{
\item{optional}{}
\item{row.names}{}
\item{x}{ object of class massvectorlist}
}
\details{ These functions create a data frame, tightly coupled  collections of variables which share many of the properties of  matrices and of lists, used as the fundamental data structure by  most of R's modeling software. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.matrix.mlist}},
\code{\link{as.matrix}},
}
\examples{
 data(mvl)
 tmp <- as.data.frame(mvl)
 names(mvl)
 plot(tmp$lengthmv,tmp$mass.Min.)
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 test<-getintcalib(mvl,cal,error=500)
 tmp<-as.data.frame(test)
 names(tmp)
}
\keyword{misc}
