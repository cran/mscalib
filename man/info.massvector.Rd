\name{info.massvector}
\alias{info.massvector}
\title{ Info Acces}
\description{ access to the info field of the massvector. }
\usage{\method{info}{massvector}(object,info,...)}
\arguments{
\item{...}{}
\item{info}{ info character. If missing function returns the current info. If not missing function returns massvector with new info field content.}
\item{object}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 info(mv1)
 mv1<-info(mv1,"testname")
}
\keyword{misc}
