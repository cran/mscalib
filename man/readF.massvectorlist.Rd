\name{readF.massvectorlist}
\alias{readF.massvectorlist}
\title{ Reads massvectorlist}
\description{ written to disk with \code{writeF.massvectorlist} }
\usage{\method{readF}{massvectorlist}(object,path,file=experiment(object),ext="txt",...)}
\arguments{
\item{...}{}
\item{ext}{}
\item{file}{ file to read. default =experiment(object)}
\item{object}{massvectorlist}
\item{path}{ path to directory.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{readF.massvector}},
\code{\link{writeF.massvectorlist}},
\code{\link{readBruker.massvectorlist}},
}
\examples{
 data( mvl )
 mvl
 writeF( mvl , "." )
 test <- readF(massvectorlist(info(mvl)),".")
 test
 unlink(paste(info(test),".txt",sep=""))
}
\keyword{misc}
