\name{calibrestat}
\alias{calibrestat}
\title{ Constructor}
\description{ Returns object of class calibrestat. Used by function \code{getrecalib.massvector}.  A calibrestat object is returned by the \code{getrecalib.massvector} method. }
\usage{calibrestat(info,...)}
\arguments{
\item{...}{ can be Coeff.Intercept,Coeff.Slope,lengthvm,PQM,tcoor.}
\item{info}{ unique identifier.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getrecalib.massvector}},
}
\examples{
 data(mv1)
 res <- getrecalib(mv1)
 print(res)
 as.vector(res)
 summary(res)
 image(res)
 plot(res)
}
\keyword{misc}
