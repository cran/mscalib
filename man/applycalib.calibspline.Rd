\name{applycalib.calibspline}
\alias{applycalib.calibspline}
\title{ External Calbiration}
\description{ Applys object of class calibspline to massvector or massvectorlist to correct for measurment errors. }
\usage{\method{applycalib}{calibspline}(object,mv,...)}
\arguments{
\item{...}{}
\item{mv}{}
\item{object}{ massvector}
}
\details{ In case of \bold{external calibration} some sample spots are only dedicated  to calibration. Calibration samples which produces equidistant  peaks, which exact masses are known, can be used to precisely  estimate the mass dependent error function. }
\value{
\item{ massvector }{ calibrated massvector}
}
\references{\item{-}{ Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyextcalib.massvectorlist}},
\code{\link{getextcalib.massvector}},
\code{\link{getextcalib.massvectorlist}},
}
\examples{
 data(mv1)
 data(ppg)
 res<- getextcalib(ppg,getPPGmasses(),error=150)
 plot(res)
 mv2<-applycalib(res,mv1)
 compare(mv1,mv2,error=300)
 rm(mv1,mv2)
 data(mvl)
 mvl<-mvl[1:100]
 res<-applycalib(res,mvl)
}
\keyword{misc}
