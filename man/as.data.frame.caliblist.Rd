\name{as.data.frame.caliblist}
\alias{as.data.frame.caliblist}
\title{ Data Frames}
\description{ Turns the caliblist into a data.frame }
\usage{\method{as.data.frame}{caliblist}(x,row.names=NULL,optional=FALSE)}
\arguments{
\item{optional}{ logical. If \code{TRUE}, setting row names is optional. }
\item{row.names}{ \code{NULL} or a character vector giving the row names for the data frame.  Missing values are not allowed.}
\item{x}{ object of class caliblist}
}
\details{ These functions create a data frame, tightly coupled  collections of variables which share many of the properties of  matrices and of lists, used as the fundamental data structure by  most of R's modeling software. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{as.matrix.caliblist}},
\code{\link[base]{as.data.frame}},
\code{\link{calibrelist}},
\code{\link{calibintlist}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 test<-getintcalib(mvl,cal,error=500)
 tmp<-as.data.frame(test)
 names(tmp)
}
\keyword{misc}
