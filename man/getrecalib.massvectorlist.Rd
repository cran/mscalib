\name{getrecalib.massvectorlist}
\alias{getrecalib.massvectorlist}
\title{ Precalibration}
\description{ Obtain the error model using the wavelength analysis of the peaklist. }
\usage{\method{getrecalib}{massvectorlist}(object,plot=FALSE,...)}
\arguments{
\item{...}{ further parameters.}
\item{object}{ massvectorlist}
\item{plot}{ logical; \code{TRUE} - A graphic showing the \eqn{F(\omega)}{F(omega)}. default \code{FALSE}}
}
\details{ Precalibration method utilizes the knowledge that masses  of peptides are in equidistant spaced clusters. The wavelength of  the massesvector can be determined as described by  Wool. The comparision of the experimental wavelength with  the theoretical one, makes possible to find an affine function  that corrects the masses. Chemical noise in the spectra may hamper  the determination of mass list frequency. The package provides a  function to filter chemical noise. }
\value{
\item{ caliblist }{ caliblist with objects of class calibrestat.}
}
\references{\item{-}{ Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.}
\item{-}{Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\seealso{\code{\link{getrecalib.massvector}},
\code{\link{applyrecalib.massvector}},
\code{\link{massvector}},
\code{\link{calibrestat}},
\code{\link{wsFilter.massvectorlist}},
}
\examples{
 data(mvl)
 mvl<-mvl[1:100]
 res <- getrecalib(mvl)
 print(res)
 summary(res)
 image(res,what="Coef.Intercept")
 image(res,what="Coef.Slope")
 plot(res)
 hist(res)
 dres<-as.data.frame(res)
 plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM")
 #create subset.
 res2<-subset(res,PQM>10)
 length(res2)
 plot(res2)
 test<-applyrecalib(mvl, res2)
 plot(test)
}
\keyword{misc}
