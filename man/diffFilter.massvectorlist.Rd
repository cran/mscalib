\name{diffFilter.massvectorlist}
\alias{diffFilter.massvectorlist}
\title{ Abundant Differences}
\description{ Removes masses from the massvector. }
\usage{\method{diffFilter}{massvectorlist}(object,listofdiffs,higher=TRUE,error=0.05,uniq=TRUE,prune=FALSE,...)}
\arguments{
\item{...}{}
\item{error}{ How much the differences can diviate from the differences given in listofdiffs}
\item{higher}{ \code{TRUE} - remove higher mass, \code{FALSE} = remove lower mass.}
\item{listofdiffs}{}
\item{object}{}
\item{prune}{ logical;default = \code{TRUE} - removes modified mass; \code{FALSE} - return modified mass.}
\item{uniq}{}
}
\details{ Removes one of the masses contributing to a mass difference given in the list of diffs.  Can be used if a variable modification are present in the massvector but can not be considered by the identification software.  Abundant intra massvector mass differences indicate the  presence of variable modifications in the data set. This  information can be used to optimize the search strategy. }
\value{
\item{ massvector }{ filtered massvector.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{diffFilter.massvector}},
\code{\link{getdiff.massvector}},
\code{\link{getdiff.massvectorlist}},
}
\examples{
 data(mvl)
 res<-getdiff(mvl,range=c(0,100))
 hist(res)
 res<-gamasses(res,abund=100)
 test<-diffFilter(mvl,res,higher=TRUE,error=0.1,uniq=TRUE)
 test
}
\keyword{misc}
