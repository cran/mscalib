\name{applyintcalib.massvectorlist}
\alias{applyintcalib.massvectorlist}
\title{ Internal Calibration}
\description{ Corrects the massvectors in the list for the error model stored in calibintstat object. }
\usage{\method{applyintcalib}{massvectorlist}(object,calc,...)}
\arguments{
\item{...}{ further params}
\item{calc}{ caliblist with objects of class calibintstat}
\item{object}{ massvectorlist}
}
\details{ Internal calibration aligns masses of  peaks to known masses and determines by linear regression a affine  function that describing the relative error. The internal  correction fails when no calibration peaks can be found. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyintcalib.massvector}},
\code{\link{getintcalib.massvectorlist}},
\code{\link{correctinternal.massvectorlist}},
\code{\link{calibintstat}},
\code{\link{caliblist}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 res <- getintcalib(mvl,cal,error=300,ppm=FALSE)
 mvl2<-applyintcalib(mvl,res)
}
\keyword{misc}
