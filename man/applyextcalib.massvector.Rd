\name{applyextcalib.massvector}
\alias{applyextcalib.massvector}
\title{ External Calbiration}
\description{ Applys object of class calibspline to massvector to correct ofr measurment errors.  The error model stored in the calibspline are obtained by the function \code{getextcalib} }
\usage{\method{applyextcalib}{massvector}(object,cS,...)}
\arguments{
\item{...}{}
\item{cS}{ calibspline}
\item{object}{ massvector}
}
\details{ In case of \bold{external calibration} some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ massvector }{ calibrated massvector}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. \emph{Anal Chem.} 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyextcalib.massvectorlist}},
\code{\link{applycalib.calibspline}},
\code{\link{getextcalib.massvector}},
\code{\link{getextcalib.massvectorlist}},
}
\examples{
 data(ppg)
 data(mv1)
 res<-getextcalib(ppg[[1]],getPPGmasses(),error=200)
 mv2<-applyextcalib(mv1,res)
 compare(mv1,mv2)
}
\keyword{misc}
