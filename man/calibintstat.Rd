\name{calibintstat}
\alias{calibintstat}
\title{ Constructor}
\description{ Returns an object of class calibintstat. Is used by function \code{getintcalib}. }
\usage{calibintstat(info, lm , ...)}
\arguments{
\item{...}{ further parameters Coef.Intercept, Coef.Slope, lengthmv, nrmatch, error.mean, error.stdv, ppm.}
\item{info}{ identifier of the calibstat object. Links it with the massvector.}
\item{lm}{ object of class 'lm' (linear model)}
}
\value{
\item{ calibintstat }{ object of class calibintstat.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link[base]{lm}},
}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=2,ppm=FALSE)
 class(res)
 plot(res)
 res<-getintcalib(mv1,cal,error=400,ppm=TRUE)
 class(res)
 plot(res)
}
\keyword{misc}
