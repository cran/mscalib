\name{getglobalcalib.massvectorlist}
\alias{getglobalcalib.massvectorlist}
\title{ Set Based Internal Calbiration}
\description{ Obtains error model using massvector with known masses. }
\usage{\method{getglobalcalib}{massvectorlist}(object, calib , error=500,ppm=TRUE, labund=12, abund = length(object)/5 , accur = ifelse(ppm,error/2000,error) ,...)}
\arguments{
\item{...}{ further parameters.}
\item{abund}{ default =  length(object)/5}
\item{accur}{ default =  ifelse(ppm,error,error/2000)}
\item{calib}{ massvector with calibration masses}
\item{error}{ assumed measurment error.}
\item{labund}{ how many abundant masses to use for calibration. default=12.}
\item{object}{ massvectorlist}
\item{ppm}{ logical; \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.}
}
\details{ Set based Calibration copes with the problem of missing calibration  masses. It first extracts about 15 most abundant masses of the  massvectorlist, then they are internally calibrated and used  as new calibration masses. In this fashion more massvectors  can be internally calibrated.  Internal calibration aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
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
 res<-getglobalcalib(mvl,cal,error=500,ppm=TRUE)
 hist(res)
}
\keyword{misc}
