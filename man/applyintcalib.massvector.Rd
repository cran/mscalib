\name{applyintcalib.massvector}
\alias{applyintcalib.massvector}
\title{ Internal Calibration}
\description{ Corrects the massvector for the error model stored in calibintstat object. }
\usage{\method{applyintcalib}{massvector}(object,cal,...)}
\arguments{
\item{...}{ further parameters.}
\item{cal}{ object of class calibintstat}
\item{object}{ massvector}
}
\details{ \bold{Internal calibration} aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvector }{ calibrated massvector. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applycalib.calibintstat}},
\code{\link{getintcalib.massvector}},
\code{\link{correctinternal.massvector}},
}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=300)
 mv2<- applyintcalib(mv1,res)
 plot(mv1[,1],mv2[,1]-mv1[,1])
}
\keyword{misc}
