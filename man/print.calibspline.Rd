\name{print.calibspline}
\alias{print.calibspline}
\title{ Print Calibspline Object}
\description{ `print' prints its argument and returns it invisibly (via `invisible(x)') }
\usage{\method{print}{calibspline}(x,...)}
\arguments{
\item{...}{}
\item{x}{ calibspline}
}
\value{
\item{ info }{ the id of the calibspline}
\item{ smooth.spline }{ smooth.spline}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(ppg)
 tmp <- getextcalib(ppg,error=200)
 print(tmp)
}
\keyword{misc}
