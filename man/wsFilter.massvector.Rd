\name{wsFilter.massvector}
\alias{wsFilter.massvector}
\title{ Smilanski Filtering}
\description{ Removes chemical noise from the massvectorlist. }
\usage{\method{wsFilter}{massvector}(object,mdist=0.25,fraction=0.2, peptides=TRUE,...)}
\arguments{
\item{...}{ further parameters.}
\item{fraction}{ maximal size of branche (nr masses in branch)/(length of massvector) to be cut.}
\item{mdist}{ minimal distance of branch to be cut.}
\item{object}{ massvector.}
\item{peptides}{ \code{TRUE} - returns peptides, \code{FALSE} - returns chemical noise.}
}
\details{ Chemical noise can be removed from the peptide mass lists  due to the strong clustering of mono-isotopic peptide  peaks. Following the distance measure and filtering  method proposed by Wool Smilanski we developed an algorithm to  classify masses as peptide and non-peptide. The algorithm is based  on a modified distance measure and hierarchical clustering of all  intra massvector distances. }
\value{
\item{ massvector }{  returns either a massvector of nonpeptide masses or massvector of peptide masses.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{wsiFilter.massvector}},
\code{\link{wsFilter.massvectorlist}},
}
\examples{
 data(mv1)
 tmp <- wsFilter(mv1,peptide=FALSE)
 plot(tmp)
 tmp <- wsFilter(mv1,peptide=TRUE)
 plot(tmp)
}
\keyword{misc}
