\name{print.calibintstat}
\alias{print.calibintstat}
\title{ Print calibintstat}
\description{ `print' prints its argument and returns it invisibly (via `invisible(x)') }
\usage{\method{print}{calibintstat}(x,...)}
\arguments{
\item{...}{}
\item{x}{ calibintstat object}
}
\value{
\item{ info }{ the id of the massvector}
\item{ type }{ type of error}
\item{ Intercept }{ Intercept of the mass dependent error}
\item{ Slope }{ The slope of the mass dependent error}
\item{ Lenght pl }{ Lenght of the peaklist}
\item{ mean }{ mean of the error}
\item{ stdv }{ stdv of the error}
\item{ nrmatch }{ nr matches.}
\item{ Xcoor }{ x coordinate on sample support}
\item{ Ycoor }{ y coordinate on sample support}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
 print(res)
}
\keyword{misc}
