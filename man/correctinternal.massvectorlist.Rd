\name{correctinternal.massvectorlist}
\alias{correctinternal.massvectorlist}
\title{ Internal Calibration}
\description{ Determines the measurment error of the masses using \code{getintcalib.massvectorlist}.  It refines the error model and applies it to the massvector by \code{applyintcalib.massvectorlist}. }
\usage{\method{correctinternal}{massvectorlist}(object,calib,error=500,ppm=TRUE,...)}
\arguments{
\item{...}{}
\item{calib}{ massvector with calibration masses}
\item{error}{ assumed measurment error.}
\item{object}{ massvectorlist}
\item{ppm}{ logical;\code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.}
}
\details{ \bold{Internal calibration} aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getintcalib.massvectorlist}},
\code{\link{correctinternal.massvector}},
\code{\link{getintcalib.massvector}},
\code{\link{calibintstat}},
\code{\link{caliblist}},
}
\examples{
 data(mvl)
 data(cal)
 mvl2 <- correctinternal( mvl,cal, error=500 , ppm=TRUE )
}
\keyword{misc}
