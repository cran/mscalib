\name{summary.calibintstat}
\alias{summary.calibintstat}
\title{ Calibintstat Summaries}
\description{ Generates a summary: info, Intercept, Slope etc. for object of class calibintastat }
\usage{\method{summary}{calibintstat}(object,...)}
\arguments{
\item{...}{ further parameters.}
\item{object}{ calibinstat.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{summary}},
}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
 summary(res)
}
\keyword{misc}
