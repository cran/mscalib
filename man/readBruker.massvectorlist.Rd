\name{readBruker.massvectorlist}
\alias{readBruker.massvectorlist}
\title{ Read Bruker XML File}
\description{ Read Bruker XML peaklist.xml file format. }
\usage{\method{readBruker}{massvectorlist}(object,folder,expname,project,filename="peaklist.xml",...)}
\arguments{
\item{...}{ further parameters}
\item{expname}{ name of the experiment. If missing the name are extracted from the folder path.}
\item{filename}{ default = "peaklist.xml"}
\item{folder}{ path to folder with Bruker file structure.}
\item{object}{ object of class massvectorlist}
\item{project}{ name of the project.}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\keyword{misc}
