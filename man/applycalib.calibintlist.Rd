\name{applycalib.calibintlist}
\alias{applycalib.calibintlist}
\title{ Internal Calibration}
\description{ Corrects the massvectors in the massvectorlist \code{mvl} using the error model  of the \code{calibintstat} objects stored in the \code{calibintlist}. }
\usage{\method{applycalib}{calibintlist}(object,mvl,...)}
\arguments{
\item{...}{ further params}
\item{mvl}{ massvectorlist}
\item{object}{ calibintlist}
}
\details{ \bold{Internal calibration} aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyintcalib.massvectorlist}},
\code{\link{getintcalib.massvectorlist}},
\code{\link{correctinternal.massvectorlist}},
\code{\link{calibintstat}},
\code{\link{caliblist}},
}
\examples{
 data(mvl)
 data(cal)
 res<-getintcalib(mvl,cal,error=300)
 mvl2<-applycalib(res,mvl)
}
\keyword{misc}
