\name{wsFilter.massvectorlist}
\alias{wsFilter.massvectorlist}
\title{ Smilanski Filtering}
\description{ Removes chemical noise from massvectors in the massvectorlist (if \code{peptides} argument \code{TRUE}) or returns it. }
\usage{\method{wsFilter}{massvectorlist}(object,mdist=0.25,fraction=0.2, peptides=TRUE,... )}
\arguments{
\item{...}{ further parameters.}
\item{fraction}{ Maximal fraction (nr masses in branch)/(length of massvector) of branche to be prune.}
\item{mdist}{ Minimal distance to branch to prune. The unit of the distance is Dalton.}
\item{object}{ massvectorlist.}
\item{peptides}{ logical; \code{TRUE} - returns peptides, \code{FALSE} - returns chemical noise.}
}
\details{ Chemical noise can be removed from the peptide mass lists  due to the strong clustering of mono-isotopic peptide  peaks. Following the distance measure and filtering  method proposed by Wool Smilanski we developed an algorithm to  classify masses as peptide and non-peptide. The algorithm is based  on a modified distance measure and hierarchical clustering of all  intra massvector distances. }
\value{
\item{ massvectorlist }{  Returns a massvectorlist where the massvectors either contain the peptides or the non-peptides.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{wsiFilter.massvector}},
\code{\link{wsFilter.massvector}},
}
\examples{
 data(mvl)
 res <- wsFilter(mvl,peptides = TRUE)
 plot(res)
 res2 <- wsFilter(mvl,peptides = FALSE)
 plot(res2,col=2,add=TRUE)
 image(res2,what="lengthmv")
 hist(res)
 hist(res2,col=2,add=TRUE)
}
\keyword{misc}
