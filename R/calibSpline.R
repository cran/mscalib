##not neccessary - its the same for each peaklist
##



calibspline <- function(ispl,mv,error,theo)
  {
    ##t Constructor
    ##- Returns object of class calibspline. It is used by the function \code{getextcalib.massvector}.
    ##n Used by by method getextcalib.
    ##+ ispl : object of class "smooth.spline"
    ##+ mv : massvector
    ##+ error : vector with errors.
    ##+ theo : vector with theoretical masses
    ##sa getextcalib.massvector
    ##ex data(ppg)
    ##ex getextcalib(ppg)
    ##ex if(!inherits(ppg,"calibspline")){stop("its not a calibspline")}
    ##ex plot(ppg)
    
    if(!inherits(ispl,"smooth.spline") || !(inherits(mv,"massvectorlist")||inherits(mv,"massvector") ))
      {
        print(class(ispl))
        print(class(mv))
        stop("param ispl not of class smooth.spline or mv not of class massvector!\n")
      }
    ispl$allow <- c("info","error","theo")
    ispl$info <- info(mv)
    if(length(error) != length(theo))
      stop("arg error and arg theo must have the same length\n")
    
    ispl$error <- error
    ispl$theo <- theo
    class(ispl)<-c("calibspline","smooth.spline","mylistobj")
    return(ispl)
  }



print.calibspline <- function(x,...)
  {
    ##t Print Calibspline Object
    ##- `print' prints its argument and returns it invisibly (via `invisible(x)')
    ##+ x : calibspline
    ##v info : the id of the calibspline
    ##v smooth.spline : smooth.spline
    ##e data(ppg)
    ##e tmp <- getextcalib(ppg,error=200)
    ##e print(tmp)
    
    cat("info :",x$info,"\n")
    res <- NextMethod("print")
    invisible(x)
  }

summary.calibspline <- function(object,...)
  {
    ##t Calibspline Summaries
    ##- Generates a summary, min, max etc.
    ##+ object : massvector
    ##sa summary
    ##e data(ppg)
    ##e tmp <- getextcalib(ppg,error=200)
    ##e summary(tmp)
    NextMethod(object)
  }


##MASSVECTOR

plot.calibspline <- function(x,...)
  {
    ##t Calibspline Plotting
    ##- Function for plotting object of class calibspline.
    ##+ x : calibspline
    ##+ ... : graphical parameters can be given as arguments to `plot'.
    ##e data(ppg)
    ##e tmp <- getextcalib(ppg,error=200)
    ##e plot(tmp)
    cS<-x
    x<-NULL
    plot(cS$theo,cS$error,pch="*",xlab="thoretical mass [m/z]",ylab="error [ppm]",main=mget(cS,"info"),...)
    lines(predict(cS,cS$x),col=2)
  }



applycalib.calibspline<-function(object,mv,...)
  {
    ##t External Calbiration
    ##- Applys object of class calibspline to massvector or massvectorlist to correct for measurment errors.
    ##d In case of \bold{external calibration} some sample spots are only dedicated
    ##d to calibration. Calibration samples which produces equidistant
    ##d peaks, which exact masses are known, can be used to precisely
    ##d estimate the mass dependent error function.
    ##+ object : massvector
    ##+ cS : calibspline
    ##v massvector : calibrated massvector
    ##sa applyextcalib.massvectorlist, getextcalib.massvector, getextcalib.massvectorlist
    ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.
    ##r Wolski http://www.molgen.mpg.de/~wolski/mscalib
    ##e data(mv1)
    ##e data(ppg)
    ##e res<- getextcalib(ppg,getPPGmasses(),error=150)
    ##e plot(res)
    ##e mv2<-applycalib(res,mv1)
    ##e compare(mv1,mv2,error=300)
    ##e rm(mv1,mv2)
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e res<-applycalib(res,mvl)
    
    if(!(inherits(mv,"massvectorlist")||inherits(mv,"massvector")))
      {
        stop(as.character(substitute(mv))," should be an object of class massvector or massvectorlist but is of class : ",class(mv),"\n")
      }
                                        #peaklist- array of peaks masses.
                                        #spline - spline to predict the error.
    if(inherits(mv,"massvector"))
      {
        error <- predict(object,mass(mv))
        error <- error$y
        masspred <- mv[,1]/(1+error/1e6)
        mv[,1] <- masspred
        return(mv)
      }
    if(inherits(mv,"massvectorlist"))
      {
        res<-lapply(mv,applyextcalib,object)
        mv <- massvectorlist(experiment(mv),res,project(mv))
        rm(res)
        return(mv)
      }
  }



#
#caliblist
#
#
#calibextlist  function(experiment,calib,...)
#  {
#    
#    res<-list()
#    tmp<-list(...)
#    allow<-c("experiment","project","data","calib")
#    attr(res,"allow") <- allow
#    if(!missing(experiment))
#      attr(res,"experiment") <- experiment
#    if(!missing(calib))
#      attr(res,"calib") <- calib
#    for(x in names(tmp))
#      {
#        if(x %in% allow)
#          {
#            attr(res,x)<-tmp[[x]]
#          }
#      }
#    class(res) <- c("calibextlist" , "caliblist" , "myobj" , "list")
#    res
#  }

