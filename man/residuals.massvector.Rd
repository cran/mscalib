\name{residuals.massvector}
\alias{residuals.massvector}
\title{ Residues of 2 massvectors.}
\description{ Computes Residues of 2 massvectors. }
\usage{\method{residuals}{massvector}(object,fmass,error=250,ppm=TRUE,uniq=TRUE,...)}
\arguments{
\item{...}{}
\item{error}{ the assumed difference between the aligned peaks either in ppm or in Da.}
\item{fmass}{ object of class massvector.}
\item{object}{ object of class massvector. Use constructor \code{massvector()}}
\item{ppm}{ logical; if TRUE then the error must be specified in ppm.}
\item{uniq}{}
}
\details{ Computes the mass differences between aligned masses of two massvectorlists. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getaccC}},
\code{\link{compare}},
}
\examples{
 data(mv1)
 data(mv2)
 plot(resid(mv1,mv2,error=250,uniq=TRUE))
}
\keyword{misc}
