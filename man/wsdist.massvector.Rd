\name{wsdist.massvector}
\alias{wsdist.massvector}
\title{ Wool Smilanski Distance Matrix}
\description{ This function computes and returns the distance matrix  using the intra massvecotor distance  measure. }
\usage{\method{wsdist}{massvector}(object,...)}
\arguments{
\item{...}{}
\item{object}{ massvector}
}
\value{
\item{ dist }{ an object of class distance.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. \emph{Proteomics.} 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{wsFilter.massvector}},
\code{\link{wsiFilter.massvector}},
}
\examples{
 data(mv1)
 plot(hclust(wsdist(mv1),method="single"))
}
\keyword{misc}
