\name{getextcalib.massvectorlist}
\alias{getextcalib.massvectorlist}
\title{ External Error Model}
\description{ Returns the error model obtained from the calibration sample. }
\usage{\method{getextcalib}{massvectorlist}(object,calib,error=250,...)}
\arguments{
\item{...}{ further parameters.}
\item{calib}{ massvector with calibration masses (ppg)}
\item{error}{ relative measurment error in ppm. (peaks in this range are assumed as matching)}
\item{object}{ massvectorlist with calibration sample (ppg) masses}
}
\value{
\item{ calibspline }{ can be used to calibrate peaklists}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{calibspline}},
\code{\link{applyextcalib.massvectorlist}},
\code{\link{getextcalib.massvector}},
}
\keyword{misc}
