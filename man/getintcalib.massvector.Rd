\name{getintcalib.massvector}
\alias{getintcalib.massvector}
\title{ Internal Calbiration.}
\description{ Obtains error model by alingning masses in massvector to known masses (calibration list). }
\usage{\method{getintcalib}{massvector}(object,calib,error=500,uniq=FALSE,ppm=TRUE,...)}
\arguments{
\item{...}{}
\item{calib}{ massvector with calibration masses}
\item{error}{ assumed measurment error.}
\item{object}{ massvector}
\item{ppm}{ \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.}
\item{uniq}{ \code{TRUE}- use only mass closest to calibration mass. \code{FALSE}- use all masses closer to the calibration mass then given error.}
}
\details{ \bold{Internal calibration} aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
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
}
\keyword{misc}
