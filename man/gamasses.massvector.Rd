\name{gamasses.massvector}
\alias{gamasses.massvector}
\title{ Abundant masses}
\description{ Determines abundant masses in a massvector. }
\usage{\method{gamasses}{massvector}(object,accur = 0.1,abund = 50,...)}
\arguments{
\item{...}{ further parameters.}
\item{abund}{ how many times a mass have to occur in a mass bin to be considered as an abundant mass.}
\item{accur}{ measurment accuracy in m/z}
\item{object}{ massvector}
}
\details{ Abundant masses are masses that occur in a large  fraction of the massvectors. Typical abundant mass  are derived from tryptic autoproteolysis products.  Abundant masses can also often be assigned to keratin  isoforms (human hair- and skin proteins).  Many of the abundant masses cannot be assigned to any protein.  Abundant masses can be used for calibration.  Removing them may increase the identification specificity. }
\value{
\item{ massvector }{ massvector with abundant masses.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{gamasses.massvectorlist}},
}
\examples{
 data(mvl)
 mv<-unlist(mvl) 
 res<-gamasses(mv,abund=30,accur=0.3)
 plot(res)
}
\keyword{misc}
