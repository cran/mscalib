\name{writeF.massvector}
\alias{writeF.massvector}
\title{ Write massvector}
\description{ Write massvector to File }
\usage{\method{writeF}{massvector}(object,path,file=info(object),ext="txt",...)}
\arguments{
\item{...}{}
\item{ext}{ file extension.}
\item{file}{ file name. defualt info(object)}
\item{object}{ massvector}
\item{path}{ path to folder.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{readBruker.massvector}},
\code{\link{readBruker.massvectorlist}},
\code{\link{writeF.massvectorlist}},
}
\examples{
 data(mv1)
 writeF(mv1,".") # writes the file in the home directory.
 readF(massvector(info(mv1)),".")
 file.remove(paste(massvector(info(mv1)),"txt",sep=""))
}
\keyword{misc}
