\name{applycalib.calibrelist}
\alias{applycalib.calibrelist}
\title{ Precalibration}
\description{ Uses the error model obtained by the method \code{getrecalib} and stored in a \code{calibrelist}  to correct masses in the massvectorlist. }
\usage{\method{applycalib}{calibrelist}(object,mvl,...)}
\arguments{
\item{...}{ further parameters}
\item{mvl}{ massvectorlist}
\item{object}{ calibrelist}
}
\details{ \bold{Precalibration} method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the \emph{massesvector} can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
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
 mvl<-mvl[1:100]
 res<-getrecalib(mvl)
 plot(res)
 image(res,what="PQM")
 mvlr<-applycalib(res,mvl)
}
\keyword{misc}
