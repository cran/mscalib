\name{globalcalib.massvectorlist}
\alias{globalcalib.massvectorlist}
\title{ Set Based Internal Calibration}
\description{ Determines the error and corrects for it. }
\usage{\method{globalcalib}{massvectorlist}(object, calib , error=500, labund = 12, ppm=TRUE , abund=length(object)/5,accur = ifelse(ppm,error/2000,error) ,...)}
\arguments{
\item{...}{ further parameters.}
\item{abund}{ default  }
\item{accur}{ default }
\item{calib}{ massvector with calibration masses}
\item{error}{ assumed measurment error.}
\item{labund}{ how many abundant masses use for calibration. (8-12 masses are sufficient).}
\item{object}{ massvectorlist}
\item{ppm}{ \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.}
}
\details{ Set based Calibration copes with the problem of missing calibration  masses. It first extracts about 15 most abundant masses of the  massvectorlist, then they are internally calibrated and used  as new calibration masses. In this fashion more massvectors  can be internally calibrated.  Internal calibration aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getintcalib.massvectorlist}},
\code{\link{correctinternal.massvectorlist}},
\code{\link{getintcalib.massvector}},
\code{\link{calibintstat}},
\code{\link{caliblist}},
}
\examples{
 data(mvl)
 data(cal)
 mvl2<-globalcalib(mvl,cal,error=500,labund=12)
}
\keyword{misc}
