\name{recalibrate.massvector}
\alias{recalibrate.massvector}
\title{ Precalibration}
\description{ Obtains the error and performs the calibration in one step. }
\usage{\method{recalibrate}{massvector}(object,PQM=7,...)}
\arguments{
\item{...}{ further parameters}
\item{PQM}{ Peak Quality Measure. Indicates how well the wavelenght of the massvector was determined.}
\item{object}{ massvector}
}
\details{ Precalibration method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. {\em Proteomics.} 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 mv2<-recalibrate(mv1)
 plot(mv1[,1],mv1[,1]-mv2[,1],type="l")
}
\keyword{misc}
