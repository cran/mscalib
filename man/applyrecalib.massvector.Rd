\name{applyrecalib.massvector}
\alias{applyrecalib.massvector}
\title{ Precalibration}
\description{ Uses the model of the measurment error stored in a calibrestat object  to correct the masses in the massvector. }
\usage{\method{applyrecalib}{massvector}(object,calc,...)}
\arguments{
\item{...}{ further arguments.}
\item{calc}{ object calibrestat class.}
\item{object}{ massvector.}
}
\details{ \bold{Precalibration} method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting.\emph{Proteomics.} 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getrecalib.massvector}},
\code{\link{getrecalib.massvectorlist}},
\code{\link{recalibrate.massvector}},
\code{\link{recalibrate.massvectorlist}},
\code{\link{calibrestat}},
\code{\link{calibrelist}},
}
\examples{
 data(mv1)
 res <- getrecalib(mv1)
 plot(res)
 mv2<-applyrecalib(mv1,res)
 compare(mv1,mv2)
}
\keyword{misc}
