\name{c.massvectorlist}
\alias{c.massvectorlist}
\title{ Combine }
\description{ Combines Massvectorlist into one Massvectorlist. }
\usage{\method{c}{massvectorlist}(mvl,...)}
\arguments{
\item{...}{ massvectorlists to be concatenated.}
\item{mvl}{ massvector}
}
\details{ It does not check if the massvectors in the list are unique. }
\value{
\item{ massvectorlist }{ massvectorlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{c}},
}
\examples{
 data(mvl)
 mvl2<-c(mvl[1:100],mvl[1:50])
 mvl2
}
\keyword{misc}
