\name{experiment.mlist}
\alias{experiment.mlist}
\alias{info.mlist},
\title{ Info Acces}
\description{ Access to the experiment field of the mlist. Can be used for setting or getting the experiment field. }
\usage{\method{experiment}{mlist}(object,experiment,...)}
\arguments{
\item{...}{}
\item{experiment}{ New experiment name. If not missing function returns massvector with new info field content.}
\item{object}{ object of class mlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 experiment(mvl)
 mvl<-experiment(mvl,"newname")
}
\keyword{misc}
