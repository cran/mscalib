\name{writeF.massvectorlist}
\alias{writeF.massvectorlist}
\title{ Write massvectorlist}
\description{ Write massvectorlist to File }
\usage{\method{writeF}{massvectorlist}(object,path,file=experiment(object),ext="txt",...)}
\arguments{
\item{...}{}
\item{ext}{ file extension. default txt.}
\item{file}{ file name. default experiment(object)}
\item{object}{ massvectorlist.}
\item{path}{ path to directory.}
}
\details{ The read and write functions for all the different peak-list formats are not provided by the package. This is because  there are oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package  which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{readF.massvectorlist}},
\code{\link{readF.massvector}},
}
\keyword{misc}
