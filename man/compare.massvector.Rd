\name{compare.massvector}
\alias{compare.massvector}
\title{ Compares massvectors}
\description{ Compares the masses in the massvectors. Returns basic statistics about the matching peaks.  Plots the relative or absolute error of matchin peaks. }
\usage{\method{compare}{massvector}(object,mv2,plot=TRUE,error=1000,ppm=TRUE,uniq=FALSE,...)}
\arguments{
\item{...}{ further parameters.}
\item{error}{ size of the measurment error (default 150 ppm)}
\item{mv2}{ massvector}
\item{object}{ massvector}
\item{plot}{ True - plot the relatvetor absolute error. \code{FALSE} - no plotting.}
\item{ppm}{ \code{TRUE} - relative error in parts per million, \code{FALSE} - absolute error.}
\item{uniq}{ logical }
}
\value{
\item{ FMSTAT }{ Fowlkes & Mallows statistik (nr matching)/sqrt(length(object)*length(mv2))}
\item{ min }{ smallest error}
\item{ ... }{ 1st qu. , mean, median, 3rd qu., max and stdv of error.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{resid}},
}
\examples{
 data(mv1)
 data(mv2)
 compare(mv1,mv2,error=5000,ppm=TRUE,uniq=TRUE)
 compare(mv2,mv1,error=1,ppm=FALSE,uniq=TRUE)
}
\keyword{misc}
