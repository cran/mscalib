\name{getrecalib.massvector}
\alias{getrecalib.massvector}
\title{ Precalibration}
\description{ Obtains the an error model using the wavelength analysis of the peaklist. }
\usage{\method{getrecalib}{massvector}(object,plot=FALSE,...)}
\arguments{
\item{...}{}
\item{object}{}
\item{plot}{ \code{TRUE} - A graphic showing the \eqn{F(\omega)}{F(omega)} function. default \code{FALSE}}
}
\details{ \bold{Precalibration} method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\value{
\item{ calibrestat }{ object of class calibrestat.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting.\emph{ Proteomics.} 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{applyrecalib.massvector}},
\code{\link{massvector}},
\code{\link{calibrestat}},
\code{\link{wsFilter.massvector}},
}
\examples{
 data(mv1)
 res <- getrecalib(mv1)
 print(res)
 as.vector(res)
 summary(res)
 image(res)
 plot(res)
}
\keyword{misc}
