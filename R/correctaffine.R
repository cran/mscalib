#Copyright 2001, W. Wolski, all rights reserved.
##calibintstat - object with calibration statistics
##calibintlist - collection with calibration data.
#this is an calibration object.

calibintstat<-function(info , lm , ...)
  {
    ##t Constructor
    ##- Returns an object of class calibintstat. Is used by function \code{getintcalib}.
    ##+ info : identifier of the calibstat object. Links it with the massvector.
    ##+ lm : object of class 'lm' (linear model)
    ##+ ... : further parameters Coef.Intercept, Coef.Slope, lengthmv, nrmatch, error.mean, error.stdv, ppm.
    ##sa \link[base]{lm}
    ##v calibintstat : object of class calibintstat.
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=2,ppm=FALSE)
    ##e class(res)
    ##e plot(res)
    ##e res<-getintcalib(mv1,cal,error=400,ppm=TRUE)
    ##e class(res)
    ##e plot(res)
    allow<-c("info","Coeff.Intercept", "Coeff.Slope","info","lengthmv","nrmatch","ppm","tcoor")

    if(missing(lm))
      {
        res<-list()
      }
    else
      {
        res<-lm
        res$Coeff.Intercept<-coef(lm)[1]
        res$Coeff.Slope<-coef(lm)[2]
      }
    tmp<-list(...)
    res$allow<- allow 
    if(!missing(info))
      res$info<-info
    for(x in names(tmp))
      {
        if(x %in% allow)
          {
            res[[x]]<-tmp[[x]]
          }
      }
    class(res)<-c("calibintstat","calibstat","lm","mylistobj")
    res
  }


as.lm.calibintstat<-function(object,...)
  {
    ##t Fitting Linear Models
    ##- attempts to coerce its argument into a linear model.
    ##+ object of class calibintstat
    ##v object of class "lm".
    al <- object$allow
    for(i in al)
      {
        object[[i]]<-NULL
      }
    object$allow<-NULL
    class(object)<-"lm"
    return(object)
  }

print.calibintstat<-function(x,...)
  {
    ##t Print calibintstat
    ##- `print' prints its argument and returns it invisibly (via `invisible(x)')
    ##+ x : calibintstat object
    ##v info : the id of the massvector
    ##v type : type of error
    ##v Intercept : Intercept of the mass dependent error
    ##v Slope : The slope of the mass dependent error
    ##v Lenght pl : Lenght of the peaklist
    ##v mean : mean of the error
    ##v stdv : stdv of the error
    ##v nrmatch : nr matches.
    ##v Xcoor : x coordinate on sample support
    ##v Ycoor : y coordinate on sample support
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
    ##e print(res)

    inf <- mget(x,"info")
    cat("info         :",inf,"\n")
    typ <- ifelse(mget(x,"ppm"),"ppm","abs")
    cat("type         :",typ,"\n")
    int <- mget(x,"Coeff.Intercept")
    cat("Intercept    :",int,"\n")
    slop <- mget(x,"Coeff.Slope")
    cat("Slope        :",slop,"\n")
    ll <- mget(x,"lengthmv")
    cat("Length pl    :",ll,"\n")
    nrmat<- mget(x,"nrmatch")
    cat("nr match     :",nrmat,"\n")
    tcoor<-mget(x,"tcoor")
    cat("tcoor        :",tcoor,"\n")
    invisible(x)
  }

summary.calibintstat<-function(object,...)
  {
    ##t Calibintstat Summaries
    ##- Generates a summary: info, Intercept, Slope etc. for object of class calibintastat
    ##+ object : calibinstat.
    ##+ ... : further parameters.
    ##sa summary
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
    ##e summary(res)

    res<-c(info(object),as.vector(object))
    names(res)[1]<-"info"
    class(res) <- "table"
    res
  }

as.vector.calibintstat<-function(x,mode="any")
  {
    ##t Coerces to vectors
    ##- 'as.vector', a generic, attempts to coerce its argument into a
    ##-  vector of mode 'mode' (the default is to coerce to whichever mode
    ##- is most convenient).  The attributes of 'x' are removed.
    ##+ x : object of class calibintstat
    ##+ mode : any
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
    ##e as.vector(res)
    temp<-summary(as.lm(x))
    res <- c(ifelse(is.null(x$lengthmv),NA, x$lengthmv ),
             ifelse(is.null(x$Coeff.Intercept),NA , x$Coeff.Intercept ),
             ifelse(is.null(x$Coeff.Slope),NA , x$Coeff.Slope ),
             ifelse(is.null(x$nrmatch),NA,x$nrmatch),
             if(is.null(x$tcoor)){c(NA,NA)}else{x$tcoor},
             temp$r.squared,
             temp$adj.r.squared,
             if(!is.null(temp$fstatistic[1]))
             {
               pf(temp$fstatistic[1],temp$fstatistic[2],temp$fstatistic[3],lower.tail=FALSE)
             }
             else
             {
               NA
             }
             )
    names(res) <- c("lengthmv","Coef.Intercept","Coef.Slope","nrmatch","Xcoor","Ycoor","R.Squared","Adjusted.R.squared","p.value")
    res
  }

applycalib.calibintstat<-function(object,mv,...)
  {
    ##t Internal Calibration
    ##- Corrects the massvector for the error model stored in calibintstat object.
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object :calibintstat
    ##+ mv : massvector
    ##+ ... : further parameters
    ##v massvector : calibrated massvector. 
    ##sa applyintcalib.massvector, getintcalib.massvector, correctinternal.massvector
    ##r Wolski
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=300)
    ##e mv2<- applycalib(res,mv1)
    ##e plot(mv1[,1],mv2[,1]-mv1[,1])
    if(!inherits(mv,"massvector"))
      {
        print(class(mv))
        stop(as.character(substitute(mv)),"should be of class massvector!\n")
      }
    if(inherits(object,"lm"))
      {
        errp <- predict(object,data.frame(masstheo=mass(mv) ) )
        if(mget(object,"ppm"))
          mv[,1] <- mv[,1]/(1 - errp/1e6)
        else
          mv[,1] <- mv[,1] + errp
      }
    mv
  }


plot.calibintstat<-function(x,...)
  {
    y<-x
    rm(x)
    ##t Calibinstat Plotting
    ##- Plots a line with intercept and slope given by the error model.
    ##+ x : object of class calibintstat
    ##+ ... : graphical parameters can be given as arguments to plot.
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
    ##e plot(res)
    ylab <- ifelse(mget(y,"ppm"),"error ppm","error abs")
    curve(y$Coeff.Intercept + y$Coeff.Slope*x,xlim=c(800,4000),xlab="m/z",ylab=ylab,main="mass dependent error functon",...)
  }


#################################################
## calibintlist
#

applycalib.calibintlist<-function(object,mvl,...)
  {
    ##t Internal Calibration
    ##- Corrects the massvectors in the massvectorlist \code{mvl} using the error model
    ##- of the \code{calibintstat} objects stored in the \code{calibintlist}.
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object: calibintlist
    ##+ mvl : massvectorlist
    ##+ ... : further params
    ##v massvectorlist : calibrated massvectorlist. 
    ##sa applyintcalib.massvectorlist, getintcalib.massvectorlist, correctinternal.massvectorlist, calibintstat, caliblist
    ##r Wolski
    ##e data(mvl)
    ##e data(cal)
    ##e res<-getintcalib(mvl,cal,error=300)
    ##e mvl2<-applycalib(res,mvl)
    

    if(!inherits(mvl,"massvectorlist"))
      stop(as.character(substitute(mvl)),"have to be a object of class massvectorlist!!!\n")
                                       #calconstants -  output of function recalibrate
    for(x in 1:length(object))
      {
        tmp <- mvl[[names(object)[x]]]
        mvl[[names(object)[x]]] <- applyintcalib(tmp,object[[x]])
        if(x%%10==0)
          cat(formatC(x,width=3)," ",sep="")
        if(x%%100==0)
          cat("\n")
      }
    cat("\n")
    mvl
  }

plot.calibintlist<-function(x,...)
  {
    ##t Calibintlist Plotting
    ##- A matrix of scatterplots is produced.
    ##+ x : object of class calibintlist
    ##+ ... : graphical parameters can be given as arguments to plot.
    ##e data(mvl)
    ##e data(cal)
    ##e mvl <- mvl[1:100]
    ##e ires <- getintcalib(mvl,cal,error=250)
    ##e plot(ires)

    dat<-as.matrix(x)
    colnames(dat)<-names(as.vector(x[[1]]))
    plot(data.frame(dat),pch="*",...)
    invisible(dat)
  }


hist.calibintlist<-function(x,...)
  {
    ##t Histogram Plot
    ##- Computes histograms of the lengthmv, Coef.Intercept, Coef.Slope ...
    ##+ x : calibrelist
    ##+ ... : further graphical parameters.
    ##e data(mvl)
    ##e data(cal)
    ##e mvl <- mvl[1:100]
    ##e data(cal)
    ##e ires <- getintcalib(mvl,cal,error=250)
    ##e hist(ires)

    dat<-as.matrix(x)
    colnames(dat)<-names(as.vector(x[[1]]))
    par(mfrow=c(3,2))
    hist(dat[,1],main=colnames(dat)[1])
    hist(dat[,2],main=colnames(dat)[2])
    hist(dat[,3],main=colnames(dat)[3])
    hist(dat[,4],main=colnames(dat)[4])
    hist(dat[,5],main=colnames(dat)[5])
    hist(dat[,6],main=colnames(dat)[6])
    par(mfrow=c(1,1))
    invisible(dat)
  }




