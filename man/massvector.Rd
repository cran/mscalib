\name{massvector}
\alias{massvector}
\title{ Constructor}
\description{ massvector extends matrix }
\usage{massvector(info,masses,tcoor)}
\arguments{
\item{info}{ a unique identifier for the list.}
\item{masses}{ a matrix or a array with masses. First column must contain masses.}
\item{tcoor}{ a array of length two if integer sample support coordinates.}
}
\details{ The class massvector keeps the masses and their intensities in a matrix.  It stores additional attributes like the coordinates of the mass spectrometric sample on the support.  It also stores a identifier of the massvector. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 massvector("hello march",NULL)
 massvector("hello march",1:100)
 massvector("hello march",cbind(1:10,10:1))
 tmp<-cbind(1:10,1:10)
 colnames(tmp)<-c("mass","test")
 rr<-massvector("hello march",tmp)
 rr2<-massvector("hello bart",cbind(1:12,1:12))
 # plot functions for massvector
 plot(rr)
 hist(rr)
 summary(rr)
 image(rr)
 info(rr)
 # setting new masses
 mass(rr)
 mass(rr,1:10)
 peaks(rr,cbind(1:10,11:20))
 # plotting with masses
 data(mv1)
 data(mv2)
 plot(mv1,mv2)
 image(mv1  ,mv2)
 summary(mv1)
 print(mv1)
 hist(mv1)
 plot(mv1)
 image(mv1,mv2,error=199,ppm=FALSE)
}
\keyword{misc}
