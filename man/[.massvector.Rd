\name{[.massvector}
\alias{[.massvector}
\title{ Extract Parts of an Massvector}
\description{ The massvector extends matrix. The masses and intensities can be accessed like in case of a matrix.  Do not use mv[1:10] (not handled properly) }
\usage{}
\arguments{
\item{i}{ row access}
\item{j}{ column access}
\item{peak}{ object from which to extract elements.}
}
\value{
\item{ xxx }{ If rows are selected then a massvector is returned. If columns are accessed than arrays are returned.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{peaks.massvector}},
\code{\link{mass.massvector}},
\code{\link{mget.massvector}},
}
\examples{
 data(mv1)
 ls()
 mv1[1:10,] # returns a massvector of length 10
 mv1[1:10,1] # the first ten masses are returned.
 mv1[,1] # the masses are returned.
 mv1[,2] # the peak area are returned.
}
\keyword{misc}
