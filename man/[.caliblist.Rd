\name{[.caliblist}
\alias{[.caliblist}
\title{ Extract Parts of a Caliblist.}
\description{ The caliblist extends list. The calibstat objects in the list can therefore be accessed like list elements. }
\usage{}
\arguments{
\item{i}{ indices of calibstat object to extract}
\item{x}{ caliblist}
}
\value{
\item{ caliblist }{ caliblist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{[<-.caliblist},
\code{\link[base]{[.list}},
}
\examples{
 data(mvl)
 data(cal)
 print(cal)
 res<-getintcalib(mvl,cal,error=300)
 res<-res[1:10]
 class(res)
 plot(res)
}
\keyword{misc}
