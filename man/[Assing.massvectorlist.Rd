\name{[Assing.massvectorlist}
\alias{[Assing.massvectorlist}
\alias{[<-.massvectorlist}
\title{ Replace Parts of Massvectorlist}
\description{ The massvectorlist extends list. The massvectors in the list can therefore be accessed like list elements. }
\usage{}
\arguments{
\item{i}{ elements to replace.}
\item{mvl}{ massvectorlist.}
\item{value}{ replace by value.}
}
\value{
\item{ massvectorlist }{ massvectorlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{[.massvector},
\code{[.list},
}
\examples{
 data(mvl)
 mvl2<-mvl
 mvl2[11:20]<-mvl[1:10]
 compare(mvl[[11]],mvl[[1]])
}
\keyword{misc}
