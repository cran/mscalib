\name{caliblist}
\alias{caliblist}
\alias{calibrelist},
\alias{calibintlist},
\title{ Constructor}
\description{ Returns object of class caliblist. }
\usage{caliblist(class,experiment,data,...)}
\arguments{
\item{...}{}
\item{class}{ string with class name "calibrestat","calibintstat".}
\item{data}{ a list with calibstat objects}
\item{experiment}{ string with experiment name.}
}
\details{ Returns objects of class calibrelist and calibintlist. The name of the class  have to be specified as the first argument in the constructor.  Constructor is used by function getrecalib.massvectorlist & getintcalib.massvectorlist }
\value{
\item{ calibrelist }{ If called with first (\code{class}) argument set to "calibrelist".}
\item{ calibintlist }{ If called with first (\code{class}) argument set to "calibintlist".}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getrecalib.massvectorlist}},
\code{\link{getintcalib.massvectorlist}},
}
\examples{
 #Example calibrelist class:
 data(mvl)
 mvl<-mvl[1:10]
 res <- getrecalib(mvl)
 print(res)
 summary(res)
 image(res,what="Coef.Intercept")
 image(res,what="Coef.Slope")
 plot(res)
 hist(res)
 dres<-as.data.frame(res)
 plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM")
 #greate subset.
 res2<-subset(res,PQM>10)
 length(res2)
 plot(res2)
 test<-applyrecalib(mvl, res2)
}
\keyword{misc}
