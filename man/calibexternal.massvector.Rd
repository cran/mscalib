\name{calibexternal.massvector}
\alias{calibexternal.massvector}
\title{ External Calbiration}
\description{ Perfroms external calibration. Obtains the error model by calling \code{getextcalib}  and corrects the masses in the massvector for the errror. }
\usage{\method{calibexternal}{massvector}(object,ppg,calib,...)}
\arguments{
\item{...}{ further parameters.}
\item{calib}{ a massvector with exact masses of the calibration sample (ppg).}
\item{object}{ massvector}
\item{ppg}{ either a massvector or massvectorlist with masses of the calibration sample (e.g. poly-(propylene glycol)).}
}
\details{ In case of \bold{external calibration} some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ massvector}{ calibrated massvector}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS.\emph{Anal Chem.} 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyextcalib.massvectorlist}},
\code{\link{getextcalib.massvector}},
\code{\link{getextcalib.massvectorlist}},
\code{\link{calibexternal.massvector}},
\code{\link{applycalib.calibspline}},
}
\examples{
 data(ppg)
 data(mv1)
 mv2<-calibexternal(mv1,ppg)
 compare(mv1,mv2)
}
\keyword{misc}
