\name{getdiff.massvectorlist}
\alias{getdiff.massvectorlist}
\title{ Massdifferences}
\description{ Computes mass differences in the mv. }
\usage{\method{getdiff}{massvectorlist}(object,masst="massrea",range=c(0,100),...)}
\arguments{
\item{...}{}
\item{masst}{}
\item{object}{}
\item{range}{}
}
\details{ Removes one of the masses contributing to a mass difference given in the list of diffs.  Can be used if a variable modification are present in the massvector but can not be considered by the identification software.  Abundant intra massvector mass differences indicate the  presence of variable modifications in the data set. This  information can be used to optimize the search strategy. }
\value{
\item{ massvector }{ filtered massvector.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getdiff.massvector}},
}
\examples{
 data(mvl)
 res<-getdiff(mvl,range=c(0,100))
 plot(res)
 tt<-gamasses(res,abund=40)
 plot(tt)
}
\keyword{misc}
