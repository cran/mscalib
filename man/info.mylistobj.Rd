\name{info.mylistobj}
\alias{info.mylistobj}
\title{ Info Acces}
\description{ access to the info field of the mylistobj }
\usage{\method{info}{mylistobj}(object,info , ...)}
\arguments{
\item{...}{}
\item{info}{ info character. If missing function returns the current info. If not missing function returns mylistobj with new info field content.}
\item{object}{ mylistobj}
}
\details{ mylistobj provides basic functionality for objects implemented using list to store the attributes. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 res <- getrecalib(mv1)
 info(res)
 res<-info(res,"testname")
}
\keyword{misc}
