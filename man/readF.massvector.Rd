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
\details{ The read and write functions for all the different peak-list formats are not provided by the package. This is because  there are oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package  which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}.  The file format of the file to be read must be:  \code{>0_A1_1SRef:1,1}  \code{842.4257236	650.66}  \code{987.3931319	180.3}  \code{...  ...}  The char between > and : is read into the info field and used as a key in the massvectorlist. So it must be unique.  After the double colon the coordinates have to written.  The first column are the masses. The second column are the intensities. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{writeF.massvector}},
}
\examples{
 data(mv1)
 writeF(mv1,".")
 readF(massvector(info(mv1)),".")
 file.remove(paste(info(mv1),".txt",sep=""))
}
\keyword{misc}
