\name{recalibrate.massvectorlist}
\alias{recalibrate.massvectorlist}
\title{ Precalibration}
\description{ Obtains the error and performs the calibration in one step. }
\usage{\method{recalibrate}{massvectorlist}(object,PQM=7,...)}
\arguments{
\item{...}{ further parameters.}
\item{PQM}{ Peak Quality Measure. Indicates how well the wavelenght of the massvector was determined.}
\item{object}{ massvectorlist.}
}
\details{ Precalibration method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\value{
\item{ massvectorlist }{ recalibrated massvectorlist.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{recalibrate.massvector}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 pp <- recalibrate(mvl,PQM=8)
 compare(pp[[1]],mvl[[1]],error=1,ppm=FALSE)
}
\keyword{misc}
