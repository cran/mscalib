\name{wsiFilter.massvector}
\alias{wsiFilter.massvector}
\title{ Smilanski Filtering}
\description{ The function returns the inidces of masses identified as chemical noise. }
\usage{\method{wsiFilter}{massvector}(object , mdist=0.25 , fraction=0.2 , ... )}
\arguments{
\item{...}{ further arguments.}
\item{fraction}{ maximal fraction (nr masses in branch)/(length of massvector) of branche  to be prune.}
\item{mdist}{ minimal distance to branch to be prune. The unit of this distance are Daltons.}
\item{object}{ massvector.}
}
\details{ Chemical noise can be removed from the peptide mass lists  due to the strong clustering of mono-isotopic peptide  peaks. Following the distance measure and filtering  method proposed by Wool Smilanski we developed an algorithm to  classify masses as peptide and non-peptide. The algorithm is based  on a modified distance measure and hierarchical clustering of all  intra massvector distances. }
\value{
\item{ indices }{  Indices of masses which are identified as being nonpeptide.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. \emph{Proteomics.} 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{wsFilter.massvector}},
\code{\link{wsdist.massvector}},
\code{\link{wsFilter.massvectorlist}},
}
\examples{
 data(mv1)
 data(mv2)
 length(mv1)
 length(wsFilter(mv1))
 length(mv2)
 length(wsFilter(mv2))
}
\keyword{misc}
