\name{correctinternal.massvector}
\alias{correctinternal.massvector}
\title{ Internal Calibration}
\description{ Corrects the masses of the massvector. It first obtains the model of the  measrument error by calling \code{getintcalib.massvector}. It than corrects the masses  by a call to \code{applyintcalib.massvector}. }
\usage{\method{correctinternal}{massvector}(object,calib,error=500,uniq=FALSE,ppm=TRUE,...)}
\arguments{
\item{...}{}
\item{calib}{ massvector with calibration masses}
\item{error}{}
\item{object}{ massvector}
\item{ppm}{}
\item{uniq}{}
}
\details{ \bold{Internal calibration} aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvector }{ calibrated massvector. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getintcalib.massvector}},
\code{\link{calibintstat}},
\code{\link{applycalib.calibintstat}},
}
\examples{
 data(mv1)
 data(cal)
 res <- correctinternal(mv1,cal,error=200)
}
\keyword{misc}
