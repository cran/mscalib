\name{readF.massvector}
\alias{readF.massvector}
\title{ Read Massvector}
\description{ Reads massvector written with the function \code{writeF.massvector} }
\usage{\method{readF}{massvector}(object,path,file=info(object),ext="txt",...)}
\arguments{
\item{...}{}
\item{ext}{ file extension; default }
\item{file}{ file name; default }
\item{object}{ object of class massvector. Use constructor \code{massvector()}}
\item{path}{ path to file.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{}
\examples{
 data(mv1)
 writeF(mv1,".")
 readF(massvector(info(mv1)),".")
 file.remove(paste(info(mv1),".txt",sep=""))
}
\keyword{misc}
