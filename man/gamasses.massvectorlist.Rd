\name{gamasses.massvectorlist}
\alias{gamasses.massvectorlist}
\title{ Abundant masses}
\description{ Determines abundant masses in a massvector. }
\usage{\method{gamasses}{massvectorlist}(object,accur=0.1,abund=50,...)}
\arguments{
\item{...}{}
\item{abund}{ how many times a mass have to occur to be an abundant mass.}
\item{accur}{ measurment accuracy}
\item{object}{ massvectorlist}
}
\details{ Abundant masses are masses that occur in a large  fraction of the massvectors. Typical abundant mass  are derived from tryptic autoproteolysis products.  Abundant masses can also often be assigned to keratin  isoforms (human hair- and skin proteins).  Many of the abundant masses cannot be assigned to any protein.  Abundant masses can be used for calibration.  Removing them may increase the identification specificity. }
\value{
\item{ massvector }{ massvector with abundant masses.}
}
\references{\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{gamasses.massvector}},
}
\examples{
 data(mvl)
 #Filtering for abundant masses.
 res<-gamasses(mvl,abund=50)
 plot(res)
 mvFilter(mvl[[1]],res)
 res2<-mvFilter(mvl,res,abundant=TRUE)
 image(res2,what="lengthmv")
 image(mvl,what="lengthmv")
 image(image(res2,what="lengthmv")/image(mvl,what="lengthmv"))
 
 hist(mvl,accur=0.3)
 hist(res2,add=TRUE,col=2,accur=0.3)
}
\keyword{misc}
