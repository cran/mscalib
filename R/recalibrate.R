##calibrestat

##calibrelist

calibrestat<-function(info,...)
  {
    ##t Constructor
    ##- Returns object of class calibrestat. Used by function \code{getrecalib.massvector}.
    ##- A calibrestat object is returned by the \code{getrecalib.massvector} method.
    ##+ info : unique identifier.
    ##+ ... : can be Coeff.Intercept,Coeff.Slope,lengthvm,PQM,tcoor.
    ##sa getrecalib.massvector
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e print(res)
    ##e as.vector(res)
    ##e summary(res)
    ##e image(res)
    ##e plot(res)

    
    res<-list()
    tmp<-list(...)
    allow<-c("info","Coeff.Intercept","Coeff.Slope","lengthmv","PQM","tcoor")
    res$allow<-allow
    if(!missing(info))
      res$info<-info
    for(x in names(tmp))
      {
        if(x %in% allow)
          {
            res[[x]]<-tmp[[x]]
          }
      }
    class(res)<-c("calibrestat","calibstat","mylistobj")
    res
  }


print.calibrestat<-function(x,...)
  {
    ##t Print Calibrestat Object
    ##- Prints the fields of the calibrestat object.
    ##+ x : object of class calibrestat
    ##sa print
    ##e data(mv1)
    ##e res<-getrecalib(mv1)
    ##e print(res)
    
    tmp<-x$allow
    cat("class           :", class(x),"\n")
    cat("info            :",mget(x,"info"),"\n")
    cat("lengthmv        :",mget(x,"lengthmv"),"\n")
    cat("Coeff.Intercept :",mget(x,"Coeff.Intercept"),"\n")
    cat("Coeff.Slope     :",mget(x,"Coeff.Slope"),"\n")
    cat("PQM             :",mget(x,"PQM"),"\n")
    cat("tcoor           :",mget(x,"tcoor"),"\n")
    invisible(x)
  }


as.vector.calibrestat<-function(x,...)
  {
    ##t Coerces to vector
    ##- coerces the calibrestat object into a vector.
    ##+ x : calibrestat
    ##e data(mv1)
    ##e res<-getrecalib(mv1)
    ##e as.vector(res)
    res<-c(
           ifelse(is.null(x$lengthmv),NA,x$lengthmv),
           ifelse(is.null(x$Coeff.Intercept),NA,x$Coeff.Intercept),
           ifelse(is.null(x$Coeff.Slope),NA,x$Coeff.Slope),
           ifelse(is.null(x$PQM),NA,x$PQM),
           if(is.null(x$tcoor)){c(NA,NA)}else{x$tcoor}
           )
    names(res)<-c("lengthmv","Coef.Intercept","Coef.Slope","PQM","Xcoor","Ycoor")
    res
  }



applycalib.calibrestat<-function(object,mv,...)
  {
    ##t Precalibration
    ##d \bold{Precalibration} method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the \emph{massesvector} can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##- Uses the error model obtained by the method \code{getrecalib} to correct masses in the massvector.
    ##+ object : calibrestat
    ##+ mv : massvector
    ##+ ...: further arguments
    ##w massvector : recalibrated massvector.
    ##sa applycalib.calibintstat, getrecalib.massvector, getrecalib.massvectorlist, recalibrate.massvector, recalibrate.massvectorlist, calibrestat, calibrelist
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. \bold{Proteomics.} 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e plot(res)
    ##e mv2<-applycalib(res,mv1)
    ##e plot(mv1[,1],mv2[,1]-mv1[,1])
    
    peak<- (mv[,1]-mget(object,"Coeff.Intercept"))*(mget(object,"Coeff.Slope")/1e6+1)
    mv[,1] <- peak
    return(mv)
  }

###############################################
## calibrelist
#


applycalib.calibrelist<-function(object,mvl,...)
{
  ##t Precalibration
  ##- Uses the error model obtained by the method \code{getrecalib} and stored in a \code{calibrelist}
  ##- to correct masses in the massvectorlist.
  ##d \bold{Precalibration} method utilizes the knowledge that masses
  ##d of peptides are in equidistant spaced clusters. The wavelength of
  ##d the \emph{massesvector} can be determined as described by
  ##d Wool. The comparision of the experimental wavelength with
  ##d the theoretical one, makes possible to find an affine function
  ##d that corrects the masses. Chemical noise in the spectra may hamper
  ##d the determination of mass list frequency. The package provides a
  ##d function to filter chemical noise.
  ##+ object : calibrelist
  ##+ mvl : massvectorlist
  ##+ ... : further parameters
  ##v massvectorlist : calibrated massvectorlist. 
  ##sa recalibrate.massvectorlist, getrecalib.massvectorlist, correctinternal.massvectorlist
  ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
  ##e data(mvl)
  ##e mvl<-mvl[1:100]
  ##e res<-getrecalib(mvl)
  ##e plot(res)
  ##e image(res,what="PQM")
  ##e mvlr<-applycalib(res,mvl)
  
   if(!inherits(mvl,"massvectorlist"))
      stop(as.character(substitute(mvl)),"have to be a object of class massvectorlist!\n")
    for(x in 1:length(object))
      {
        nami<-names(object)[x]
        tmp <- mvl[[nami]]
        tmp<-applyrecalib(tmp,object[[x]])
        mvl[[nami]]<-tmp
        if(x%%10==0)
          cat(formatC(x,width=3)," ",sep="")
        if(x%%100==0)
          cat("\n")
      }
    cat("\n")
    mvl
}

plot.calibrelist<-function(x,...)
  {
    ##t Plot
    ##- Shows a matrix of scatterplots for different varialbes..
    ##+ x : calibrelist
    ##+ ... : graphical parameters can be given as arguments to plot.
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e res2 <- getrecalib(mvl) # get recalibration model for not filtered data
    ##e plot(res2)

    dat<-as.matrix(x)[,1:4]
    colnames(dat)<-names(as.vector(x[[1]]))[1:4]
    plot(data.frame(dat),pch="*",...)
    par(mfrow=c(1,1))
  }
  
hist.calibrelist<-function(x,...)
  {
    ##t Histogram Plot
    ##- Computes Histograms of the lengthmv, Coef.Intercept, Coef.Slope, PQM.
    ##+ x : calibrelist
    ##+ ... : further graphical parameters.
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e res2 <- getrecalib(mvl) # get recalibration model for not filtered data
    ##e hist(res2)
    
    dat<-as.matrix(x)
    colnames(dat)<-names(as.vector(x[[1]]))
    par(mfrow=c(2,2))
    par(cex.main=0.6)
    hist(dat[,1],main="Recalibration",xlab="lengthmv",...)
    hist(dat[,2],main="Recalibration",xlab="Coef.Intercept",...)
    hist(dat[,3],main="Recalibration",xlab="Coef.Slope",...)
    hist(dat[,4],main="Recalibration",xlab="PQM",...)
    par(cex.main=1)
    par(mfrow=c(1,1))
  }
        


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
