\name{diffFilter.massvector}
\alias{diffFilter.massvector}
\title{ Abundant Differences}
\description{ Removes mass differences from the massvector. }
\usage{\method{diffFilter}{massvector}(object,listofdiffs,higher=TRUE,error=0.05,uniq=TRUE,prune=TRUE,...)}
\arguments{
\item{...}{ further parameters.}
\item{error}{ How much the differences can diviate from the differences given in listofdiffs.}
\item{higher}{ logical;\code{TRUE} - remove higher mass, \code{FALSE} = remove lower mass.}
\item{listofdiffs}{}
\item{object}{ massvector}
\item{prune}{ logical;\code{TRUE} - remove mass, \code{FALSE} = return modified mass.}
\item{uniq}{}
}
\details{ Removes one of the masses contributing to a mass difference given in the list of differences.  Can be used if a variable modification are present in the massvector but can not be considered by the identification software. It also can be used to return the modified masses. }
\value{
\item{ massvector }{ filtered massvector.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getdiff.massvector}},
\code{\link{getdiff.massvectorlist}},
\code{\link{diffFilter.massvectorlist}},
}
\examples{
 data(mv1)
 res<-getdiff(mv1,range=c(0,100))
 diffFilter(mv1,res,higher=TRUE)
 diffFilter(mv1,res,higher=FALSE)
}
\keyword{misc}
