\name{getextcalib.massvector}
\alias{getextcalib.massvector}
\title{ External Error Model}
\description{ Returns the error model obtained from the calibration sample. }
\usage{\method{getextcalib}{massvector}(object,calib,error=300,...)}
\arguments{
\item{...}{}
\item{calib}{ massvector with calibration masses}
\item{error}{ relative measurment error in ppm.}
\item{object}{ massvector}
}
\details{ In case of \bold{external calibration} some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ calibspline }{ can be used to calibrate peaklists}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. \emph{Anal Chem.} 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{calibspline}},
\code{\link{applyextcalib.massvector}},
}
\examples{
 data(mv1)
 data(ppg)
 res<-getextcalib(ppg[[1]])
 mv2 <- applycalib(res,mv1)
 plot(mv1[,1],mv2[,1]-mv1[,1])
}
\keyword{misc}
