\name{calibexternal.massvectorlist}
\alias{calibexternal.massvectorlist}
\title{ External Calbiration}
\description{ Perfroms external calibration of massvectorlist. Obtains the model of the measurment error by \code{getextcalib} and corrects the masses for this error. }
\usage{\method{calibexternal}{massvectorlist}(object,ppg,calib,error=300,...)}
\arguments{
\item{...}{ further parameters.}
\item{calib}{ a massvector with exact masses for the massvectors of the calibration samples (ppg).}
\item{error}{ relative error of the measurment in (ppm).}
\item{object}{ massvectorlist}
\item{ppg}{ either a massvector or massvectorlist with masses of the calibration sample (e.g poly-(propylene glycol) ppg)}
}
\details{ In case of external calibration some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ massvectorlist}{ calibrated massvector}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{calibexternal.massvector}},
\code{\link{applyextcalib.massvectorlist}},
\code{\link{getextcalib.massvector}},
\code{\link{getextcalib.massvectorlist}},
\code{\link{calibexternal.massvector}},
\code{\link{applycalib.calibspline}},
}
\examples{
 data(mvl)
 data(ppg)
 res <- calibexternal(mvl,ppg,getPPGmasses(),error=150)
}
\keyword{misc}
