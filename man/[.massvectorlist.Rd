\name{[.massvectorlist}
\alias{[.massvectorlist}
\title{ Extract Parts of an Massvectorlist.}
\description{ The massvectorlist extends list. The massvectors in the massvectorlist can therefore be accessed like list elements. }
\usage{}
\arguments{
\item{i}{ indices of massvectors to extract}
\item{mvl}{ massvectorlist}
}
\value{
\item{ massvectorlist }{ massvectorlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{[<-.massvectorlist},
\code{[<-.list},
}
\examples{
 data(mvl)
 mvl<-mvl[1:10] # returns a massvectorlist of length 10
 mvl
 class(mvl)
}
\keyword{misc}
