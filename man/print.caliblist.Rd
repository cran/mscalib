\name{print.caliblist}
\alias{print.caliblist}
\title{ Print caliblist}
\description{ `print' prints its argument and returns it invisibly (via `invisible(x)') }
\usage{\method{print}{caliblist}(x,...)}
\arguments{
\item{...}{}
\item{x}{ caliblist}
}
\value{
\item{ list }{ list}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 data(cal)
 test <- getintcalib(mvl,cal,error=500)
 print(test)
}
\keyword{misc}
