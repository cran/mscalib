#Copyright 2004, W. Wolski, all rights reserved.
############################################
##calibstat
#
summary.calibstat<-function(object,...)
  {
    ##t Calibstat summary
    ##- Generates Summary.
    ##+ object : object of class calibstat.
    ##e data(mv1)
    ##e data(cal)
    ##e test<-getintcalib(mv1,cal)
    ##e summary(test)
    cat("info :",info(object),"\n")
    print(as.vector(object))
    invisible(c(info(object),as.vector(object)))
  }

hist.calibstat <- function(x,...)
  {
    ##t Histogram
    ##- Histogram
    ##+ x : object of class calibstat
    warning("Not implemented!!\n")
  }
print.calibstat <- function(x,...)
  {
    ##t print
    ##- print
    ##+ x : object of class calibstat
    warning("Not implemented!!\n")
  }
plot.calibstat <- function(x,...)
  {
    ##t plot
    ##- plot
    ##+ x : object of class calibstat
    warning("Not implemented!!\n")
  }


################################################
#caliblist
#


plot.caliblist <- function(x,...)
  {
    ##t plot
    ##- plot
    ##+ x : object of class caliblist
    warning("not implemented yet!!\n")
  }





print.caliblist <- function(x,...)
  {
    ##t Print caliblist
    ##- `print' prints its argument and returns it invisibly (via `invisible(x)')
    ##+ x : caliblist
    ##v list : list
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e data(cal)
    ##e test <- getintcalib(mvl,cal,error=500)
    ##e print(test)
    cl<-class(x)[1]
    cat("Class  :", cl ,"\n")
    exp <- mget(x,"experiment")
    cat("Experiment :", exp ,"\n")
    pro <- mget(x,"project")
    cat("Project    :", pro ,"\n")
    ll<- length(x)
    cat("Calibstat List lenght : ",ll ,"\n")
     if(length(x)>0)
      {
        cmv<- class(x[[1]])[1]
        cat("Class calibstat object: ", cmv ,"\n")
        allow<- x[[1]][["allow"]]
        cat("Fields in calibstat objects: ", join(allow,sep=" ") ,"\n")
      }
    else{
      cmv<-NULL
      allow<-NULL
    }
    invisible(list(class=cl,experiment=exp,project=pro,length=ll,classcalibstat=cmv,allow=allow))
  }

"[[<-.caliblist" <- function(x,i,value)
  {
    ##t Replace calibstat object in the caliblist.
    ##- Replace a calibstat object in the caliblist with a different one.
    ##+ x : caliblist
    ##+ i : index or name (info) of calibstat object to replace
    ##+ value : object of class calibstat (e.g. calibrestat, calibextstat).
    ##e data(mvl)
    ##e data(cal)
    ##e res<-getintcalib(mvl,cal,error=300)
    ##e res[[1]]<-res[[10]]
    
    if(!inherits(value,"calibstat"))
      stop("only calibstats allowed to assing")
    x <- NextMethod("[[<-")
    if(is.numeric(i))
      {
        names(x)[i]<- mget(value,"info")
      }
    x
  }

"[.caliblist"<-function(x,i)
  {
    ##t Extract Parts of a Caliblist.
    ##- The caliblist extends list. The calibstat objects in the list can therefore be accessed like list elements.
    ##+ x : caliblist
    ##+ i : indices of calibstat object to extract
    ##v caliblist : caliblist
    ##sa [<-.caliblist, \link[base]{[.list}
    ##e data(mvl)
    ##e data(cal)
    ##e print(cal)
    ##e res<-getintcalib(mvl,cal,error=300)
    ##e res<-res[1:10]
    ##e class(res)
    ##e plot(res)
    
    tmp<-NextMethod("[")
    res<- caliblist(class(x)[1],mget(x,"experiment"),tmp,project=mget(x,"project"))
    res
  }


caliblist <- function(class,experiment,data,...)
  {
    ##t Constructor
    ##- Returns object of class caliblist.
    ##a calibrelist, calibintlist
    ##d Returns objects of class calibrelist and calibintlist. The name of the class
    ##d have to be specified as the first argument in the constructor.
    ##d Constructor is used by function getrecalib.massvectorlist & getintcalib.massvectorlist
    ##+ class : string with class name "calibrestat","calibintstat".
    ##+ experiment : string with experiment name.
    ##+ data : a list with calibstat objects
    ##v calibrelist : If called with first (\code{class}) argument set to "calibrelist".
    ##v calibintlist : If called with first (\code{class}) argument set to "calibintlist".
    ##sa getrecalib.massvectorlist, getintcalib.massvectorlist
    ##e #Example calibrelist class:
    ##e data(mvl)
    ##e mvl<-mvl[1:10]
    ##e res <- getrecalib(mvl)
    ##e print(res)
    ##e summary(res)
    ##e image(res,what="Coef.Intercept")
    ##e image(res,what="Coef.Slope")
    ##e plot(res)
    ##e hist(res)
    ##e dres<-as.data.frame(res)
    ##e plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM")
    ##e #greate subset.
    ##e res2<-subset(res,PQM>10)
    ##e length(res2)
    ##e plot(res2)
    ##e test<-applyrecalib(mvl, res2)

    if(missing(data))
      {
        res<-list()
      }
    else
      {
        if(!inherits(data[[1]],"calibstat"))
          stop(as.character(substitute(data))," must be a list of calibstat objects")
        res<-data
      }
    allow <- c("experiment","project","data")
    tmp<-list(...)
    for(x in names(tmp))
      {
        if(x %in% allow)
          {
            attr(res,x) <- tmp[[x]]
          }
      }
    attr(res,"allow") <- allow
    if(!missing(experiment))
      attr(res,"experiment") <- experiment

    if(missing(class))
      {
        class(res) <- c("caliblist","mlist","list","myobj")
      }
    else
      {
        class(res) <- c(class,"caliblist","mlist","list","myobj")
      }
    res
  }
