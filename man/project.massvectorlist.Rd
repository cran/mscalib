\name{project.massvectorlist}
\alias{project.massvectorlist}
\title{ Acces}
\description{ Access to the project field of the massvectorlist. }
\usage{\method{project}{massvectorlist}(object,project,...)}
\arguments{
\item{...}{ further arguments}
\item{object}{ massvectorlist}
\item{project}{ info character. If missing function returns the current info. If not missing function returns massvector with new project field content.}
}
\details{ Can be used for setting or getting the project field. }
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{info.mlist}},
\code{\link{experiment.mlist}},
}
\examples{
 data(mvl)
 project(mvl)
 mvl <- project(mvl,"newprojectname")
}
\keyword{misc}
