\name{summary.calibspline}
\alias{summary.calibspline}
\title{ Calibspline Summaries}
\description{ Generates a summary, min, max etc. }
\usage{\method{summary}{calibspline}(object,...)}
\arguments{
\item{...}{}
\item{object}{ massvector}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{summary}},
}
\examples{
 data(ppg)
 tmp <- getextcalib(ppg,error=200)
 summary(tmp)
}
\keyword{misc}
