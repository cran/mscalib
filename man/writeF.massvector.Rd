\name{writeF.massvector}
\alias{writeF.massvector}
\title{ Write massvector}
\description{ Write massvector to File. }
\usage{\method{writeF}{massvector}(object,path,file=info(object),ext="txt",...)}
\arguments{
\item{...}{}
\item{ext}{ file extension.}
\item{file}{ file name. defualt info(object)}
\item{object}{ massvector}
\item{path}{ path to folder.}
}
\details{ The read and write functions for different peak-list formats are not provided by the package. This is because  there are a oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package  which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{writeF.massvectorlist}},
\code{\link{readF.massvector}},
}
\examples{
 data(mv1)
 writeF(mv1,".") # writes the file in the home directory.
 readF(massvector(info(mv1)),".")
 file.remove(paste(info(mv1),".txt",sep=""))
}
\keyword{misc}
