\name{applyrecalib.massvectorlist}
\alias{applyrecalib.massvectorlist}
\title{ Recalibration}
\description{ Corrects the massvector for model of the measurment error stored in the calibrestat object. }
\usage{\method{applyrecalib}{massvectorlist}(object,calc,...)}
\arguments{
\item{...}{ further parameters.}
\item{calc}{ caliblist with objects of class calibretstat.}
\item{object}{ massvectorlist.}
}
\details{ Precalibration method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\value{
\item{ massvectorlist }{ calibrated massvectorlist. }
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{recalibrate.massvectorlist}},
\code{\link{getrecalib.massvectorlist}},
\code{\link{correctinternal.massvectorlist}},
}
\examples{
 data(mvl)
 res<-getrecalib(mvl)
 summary(res)
 mvl2<-applyrecalib(mvl,res)
}
\keyword{misc}
