\name{applyextcalib.massvectorlist}
\alias{applyextcalib.massvectorlist}
\title{ External Calbiration}
\description{ Corrects the massvectorlist for the measurment error stored calibspline object. }
\usage{\method{applyextcalib}{massvectorlist}(object,cS,...)}
\arguments{
\item{...}{ further parameters}
\item{cS}{ object of class calibspline}
\item{object}{ massvectorlist}
}
\details{ In case of external calibration some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyextcalib.massvector}},
\code{\link{applyextcalib.massvectorlist}},
\code{\link{getextcalib.massvector}},
\code{\link{getextcalib.massvectorlist}},
}
\examples{
 data(ppg)
 data(mvl)
 res <- getextcalib(ppg,getPPGmasses(),error=150)
 mvl2 <- applyextcalib(mvl,res)
}
\keyword{misc}
