\name{getaccC}
\alias{getaccC}
\title{ Find Matching Masses}
\description{ Returns the indices of masses matching to each other. }
\usage{getaccC(pl,cal,error=500,ppm=TRUE,uniq=FALSE)}
\arguments{
\item{cal}{ vector with masses}
\item{error}{ error in Da or ppm.}
\item{pl}{ vector with masses}
\item{ppm}{ \code{TRUE} - then error has to be given as relative error in ppm. \code{FALSE} - error are absolute error in dalton.}
\item{uniq}{ \code{TRUE} - return only closest match to mass in cal. \code{FALSE} - return all matches in the error range.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 getaccC(1001:1010,1001:1010,error=300,ppm=TRUE,uniq=TRUE)
}
\keyword{misc}
