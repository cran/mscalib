\name{massvectorlist}
\alias{massvectorlist}
\title{ Constructor}
\description{ constructor of object massvectorlist (extends list). }
\usage{massvectorlist(experiment,data,project,...)}
\arguments{
\item{...}{}
\item{data}{ list of massvectors}
\item{experiment}{ name of experiment}
\item{project}{ name of project.}
}
\value{
\item{ massvectorlist }{ object of class massvector.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 # testing constructor.
 massvectorlist("my1experiment")
 data(mvl)
 massvectorlist("my2experiment",mvl,"hello project")
 plot(mvl)
 summary(mvl)
 hist(mvl)
 hist(mvl)
 image(mvl,what="lengthmv")
 mvl2<-mvl[1:100]
 plot(mvl2)
 summary(mvl2)
 hist(mvl2)
 image(mvl2,what="lengthmv")
 #testing assingments
 mvl2[[11]]<-mvl2[[1]]
 plot(mvl2[[11]],mvl2[[1]])
 image(mvl2[[11]],mvl2[[1]])
 #make one massvector out of the peaklist
 tt<-unlist(mvl)
 plot(tt)
}
\keyword{misc}
