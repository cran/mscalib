\name{getintcalib.massvectorlist}
\alias{getintcalib.massvectorlist}
\title{ Internal Calbiration}
\description{ Obtains error model using massvector with known masses. }
\usage{\method{getintcalib}{massvectorlist}(object,calib,error=500,ppm=TRUE, ... )}
\arguments{
\item{...}{}
\item{calib}{ massvector with calibration masses}
\item{error}{ assumed measurment error.}
\item{object}{ massvectorlist}
\item{ppm}{ \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.}
}
\details{ Internal calibration aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ calibintstat }{ object of class calibintstat. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyintcalib.massvector}},
\code{\link{getintcalib.massvector}},
\code{\link{correctinternal.massvectorlist}},
\code{\link{calibintstat}},
\code{\link{caliblist}},
}
\examples{
 data(mvl)
 data(cal)
 res<-getintcalib(mvl,cal,error=400,ppm=TRUE)
 plot(res)
}
\keyword{misc}
