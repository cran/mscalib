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
\details{ The read and write functions for all the different peak-list formats are not provided by the package. This is because  there are oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package  which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{readF.massvector}},
\code{\link{writeF.massvectorlist}},
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
