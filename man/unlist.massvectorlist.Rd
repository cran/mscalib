\name{unlist.massvectorlist}
\alias{unlist.massvectorlist}
\title{ Flatten massvectorlist}
\description{  Given a list structure `x', `unlist' simplifies it to produce a   massvector which contains all the atomic components which occur in x. }
\usage{\method{unlist}{massvectorlist}(x,...)}
\arguments{
\item{...}{ further arguments.}
\item{x}{ massvectorlist.}
}
\value{
\item{ massvector }{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link[base]{unlist}},
}
\examples{
 data(mvl)
 amv <- unlist(mvl)
 length(amv)
}
\keyword{misc}
