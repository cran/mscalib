#Copyright 2004, W. Wolski, all rights reserved.
                                        #massvector list.
massvectorlist <- function(experiment,data,project,...)
  {
    ##t Constructor
    ##- constructor of object massvectorlist (extends list).
    ##+ experiment : name of experiment
    ##+ data : list of massvectors
    ##+ project : name of project.
    ##v massvectorlist : object of class massvector.
    ##e # testing constructor.
    ##e massvectorlist("my1experiment")
    ##e data(mvl)
    ##e massvectorlist("my2experiment",mvl,"hello project")
    ##e plot(mvl)
    ##e summary(mvl)
    ##e hist(mvl)
    ##e hist(mvl)
    ##e image(mvl,what="lengthmv")
    ##e mvl2<-mvl[1:100]
    ##e plot(mvl2)
    ##e summary(mvl2)
    ##e hist(mvl2)
    ##e image(mvl2,what="lengthmv")
    ##e #testing assingments
    ##e mvl2[[11]]<-mvl2[[1]]
    ##e plot(mvl2[[11]],mvl2[[1]])
    ##e image(mvl2[[11]],mvl2[[1]])
    ##e #make one massvector out of the peaklist
    ##e tt<-unlist(mvl)
    ##e plot(tt)
    
    if(!missing(data))
      {
        res <- data
      }
    else
      res <- list()
    attr(res,"allow") <- c("data","experiment","project")
    if(!missing(project))
      attr(res,"project") <- project #project.
    if(!missing(experiment))
      attr(res,"experiment") <- experiment # the big project may consists of several experiments
    class(res) <- c("massvectorlist","mlist","list","myobj")
    return(res)
  }

c.massvectorlist <- function(mvl,...)
  {
    ##t Combine 
    ##- Combines Massvectorlist into one Massvectorlist.
    ##d It does not check if the massvectors in the list are unique.
    ##+ mvl : massvector
    ##+ ... : massvectorlists to be concatenated.
    ##sa c
    ##v massvectorlist : massvectorlist
    ##e data(mvl)
    ##e mvl2<-c(mvl[1:100],mvl[1:50])
    ##e mvl2
    
    tmp <- list(...)
    for(x in tmp)
      {
        if(!inherits(x,"massvectorlist"))
          stop("All Arguments have to be a massvectorlist!!!")
      }
    res <- NextMethod("c")
    res<-massvectorlist(info(mvl),res,mget(mvl,"project"))
    res
  }

  
print.massvectorlist<-function(x,...)
{
  ##t Print massvectorlist
  ##- `print' prints its argument and returns it invisibly (via `invisible(x)')
  ##+ x : massvectorlist
  ##sa \link[base]{print}
  ##e data(mvl)
  ##e print(mvl)
  exp <- mget(x,"experiment")
  print(paste("Experiment : ", ifelse(is.null(exp),"()",paste("(",exp,")",sep="")),sep=""),quote=FALSE,...)
  proj<-mget(x,"project")
  print(paste("Project    : ", ifelse(is.null(proj),"()",paste("(",proj,")",sep="")),sep=""),quote=FALSE,...)
  print(paste("Nr (mv)    : ", paste("(",length(x),")",sep=""),sep=""),quote=FALSE,...)
  invisible(x)
}





as.list.massvectorlist<-function(x,...)
  {
    ##t List
    ##- Turns the massvectorlist into a list
    ##+ x : massvectorlist.
    ##+ ... : further parameter.
    ##e data(mvl)
    ##e res<-as.list(mvl)
    ##e class(res)
    al <- attr(x,"allow")
    print(al)
    for(u in al)
      {
        attr(x,u)<-NULL
      }
    class(x)<-"list"
    x
  }



project.massvectorlist <- function(object,project,...)
  {
    ##t Acces
    ##- Access to the project field of the massvectorlist.
    ##d Can be used for setting or getting the project field.
    ##+ object : massvectorlist
    ##+ project : info character. If missing function returns the current info. If not missing function returns massvector with new project field content.
    ##+ ... : further arguments
    ##sa info.mlist, experiment.mlist
    ##e data(mvl)
    ##e project(mvl)
    ##e mvl <- project(mvl,"newprojectname")

    if(missing(project))
      return(mget(object,"project"))
    else
      {
        setParms(object)<-list(project=project)
      }
    object
  }

histLengths <- function(x,main=info(x),xlab="length of massvectors",...)
  {
    ##t Histograms
    ##- Histogram of the massvector lengths.
    ##+ x : massvectorlist
    ##+ main : info(x)
    ##+ xlab : length of massvectors.
    ##+ ... : further parameters.
    ##sa \link[base]{hist}

    ll <- lapply(mvl,length)
    res <- hist.default(unlist(ll),xlab="massvector length",main=main,border=1,...)
  }


hist.massvectorlist<-function(x,accur = 0.1, main=info(mvl) ,xlab="m/z",xlim=c(700,4500),add=FALSE,col=1,...)
  {
    ##t Histograms
    ##- Histogram of mass frequencies in the massvectorlist.
    ##+ x : massvectorlist
    ##+ accur : bin width for plotting mass frequencies.
    ##+ col : color of the histogram,
    ##+ main : title of graph
    ##+ xlim : the range to be encompassed by the x axis.
    ##+ xlab : a title for the x axis.
    ##+ add : logical; If \code{TRUE} add to already existing plot.
    ##+ ... : further parameters.
    ##sa hist.massvector
    ##e data(mvl)
    ##e hist(mvl)
    mvl<-x
    rm(x)
    dat<-unlist(mvl)
                                        #attributes(dat)<-NULL
    mhist<-list(NULL)
                                        #assign indices of bins with high peak abundance
    mhist[[1]] <- hist(mass(dat),breaks=seq(min(mass(dat))-accur/2,max(mass(dat))+1.5*accur,accur),plot=TRUE,main=main,xlab=xlab,xlim=xlim,add=add,col=col,border=col,...)
    mhist[[2]] <- hist(mass(dat),breaks=seq(min(mass(dat))-accur,max(mass(dat))+accur,accur),add=TRUE,col=col,border=col)
  }

plot.massvectorlist <-function(x, main=info(mvl) ,xlab="m/z",xlim=c(700,4500),add=FALSE,col=1,cex=0.5,...)
{
  ##t Massvectorlist Plotting
  ##- Plots masses (m/z) in massvector against sample in massvectorlist.
  ##+ x : massvectorlist
  ##+ main : an overall title for the plot.
  ##+ xlab : a title for the x axis.
  ##+ xlim : the range to be encompassed by the x axis.
  ##+ add : \code{TRUE} - masses of a new massvectorlist are added to existing plot.
  ##+ col : color of the points denoting the masses.
  ##+ ... : graphical parameters can be given as arguments to `plot'.
  ##sa hist.massvectorlist, plot.massvector, image.mlist, plot, par
  ##e data(mvl)
  ##e plot(mvl,col= 3 )
  ##e tt<- gamasses(mvl,abund=50)
  ##e plot(mvl,col=1, xlim=c(tt[1,1]-0.4,tt[1,1]+0.4))
  mvl<-x
  rm(x)
  par(bg="lightgray")
  if(add)
    {
      for( x in 1:length(mvl))
        {
          points(mass(mvl[[x]]), rep(x,length(mass(mvl[[x]]))),pch=15,cex=cex,col=col)
        }
    }
  else
    {
      
      allm  <- unlist(mvl)
      mmin  <- min(allm)[1]
      mmax <- max(allm)[1]
      lod<-mvl
      if(is.null(xlim))
        {

          plot.default(mass(lod[[1]]), rep(1,length(mass(lod[[1]]))) , xlim = c(mmin,mmax) , ylim = c(1,length(lod)),pch=15,cex=cex,xlab=xlab,ylab="sample",main=main,col=col,...)
        }
      else
        {
          plot.default(mass(lod[[1]]), rep(1,length(mass(lod[[1]]))) , xlim = xlim , ylim = c(1,length(lod)) , pch=15, cex=cex , xlab=xlab , ylab="sample" , main=main,col=col,...)
        }
      for(x in 2:length(lod))
        {
          points(mass(lod[[x]]), rep(x,length(mass(lod[[x]]))),pch=15,cex=cex)
        }
      abline( h = seq(0,length(lod),10) , lty=3)
    }
}


"[.massvectorlist"<-function(mvl,i)
  {
    ##t Extract Parts of an Massvectorlist.
    ##- The massvectorlist extends list. The massvectors in the massvectorlist can therefore be accessed like list elements.
    ##+ mvl : massvectorlist
    ##+ i : indices of massvectors to extract
    ##v massvectorlist : massvectorlist
    ##sa [<-.massvectorlist, [<-.list
    ##e data(mvl)
    ##e mvl<-mvl[1:10] # returns a massvectorlist of length 10
    ##e mvl
    ##e class(mvl)

 
    tmp<-NextMethod("[");
    res<- massvectorlist(experiment(mvl),tmp,project(mvl))
    res
  }


"[<-.massvectorlist"<-function(mvl,i,value)
  {
    ##t Replace Parts of Massvectorlist
    ##- The massvectorlist extends list. The massvectors in the list can therefore be accessed like list elements.
    ##+ mvl : massvectorlist.
    ##+ i : elements to replace.
    ##+ value : replace by value.
    ##sa [.massvector,[.list
    ##v massvectorlist : massvectorlist
    ##e data(mvl)
    ##e mvl2<-mvl
    ##e mvl2[11:20]<-mvl[1:10]
    ##e compare(mvl[[11]],mvl[[1]])
    
    if(!inherits(value,"massvectorlist"))
      {
        stop(as.character(substitute(value)),"not a massvectorlist")
      }
    NextMethod("[<-")
  }

mvFilter.massvectorlist <-  function(object,fmass,match=FALSE,error=250,ppm=TRUE,uniq=FALSE,...)
{
  ##t Filtering Massvector
  ##- Filters massvector for masses given in a second massvector.
  ##+ object : massvectorlist
  ##+ abund : massvector
  ##+ error : mesurment error.
  ##+ ppm : given either in ppm (\code{TRUE}) or as absolut error (\code{FALSE}).
  ##+ match : logical;\code{TRUE} - than returns masses matching to the masses in massvector abundant, \code{FALSE} - returns masses not matching.
  ##+ uniq : logical; \code{FALSE} - returns or removes all masses in the range given by error. \code{TRUE} - returns ore removes only the closest mass.
  ##+ ... : further parameters.
  ##v massvectorlist : with massvectors with matching or not matching masses.
  ##sa mvFilter.massvector
  ##e data(mvl)
  ##e mvl<-mvFilter(mvl,mvl[[1]],match=FALSE,error=250)
  ##e length(mvl[[1]])
   res<-lapply(object,mvFilter,fmass,match=match,error=error,ppm=ppm,uniq=uniq)
   res<-  massvectorlist(experiment(object),res,project(object))
   res
}
  




unlist.massvectorlist<-function(x,...)
  {
    ##t Flatten massvectorlist
    ##-  Given a list structure `x', `unlist' simplifies it to produce a
    ##-  massvector which contains all the atomic components which occur in x.
    ##+ x : massvectorlist.
    ##+ ... : further arguments.
    ##sa \link[base]{unlist}
    ##v massvector : massvector
    ##e data(mvl)
    ##e amv <- unlist(mvl)
    ##e length(amv)
                                        #     tmp<-NextMethod("unlist")
    tmp<-NULL
    for(u in x)
      {
        if(length(u)!=0)
          tmp<-rbind(tmp,u)
      }
    tmp<-massvector(info(x),tmp)
    tmp
  }



"[[<-.massvectorlist"<-function(x,i,value)
  {
    ##t Replace Parts of an Object
    ##- Replace a massvector in the massvectorlist with a different one.
    ##+ x : massvectorlist
    ##+ i : index or name (info) of massvector to replace
    ##+ value : massvector
    ##e data(mvl)
    ##e data(mv1)
    ##e mvl[[10]]<-mv1
    
    
    if(!inherits(value,"massvector"))
      stop("only massvectors allowed to assing")
    x <- NextMethod("[[<-")
    if(is.numeric(i))
      names(x)[i]<-mget(value,"info")
    x
  }

  
wsFilter.massvectorlist<-function(object,mdist=0.25,fraction=0.2, peptides=TRUE,... )
  {
    ##t Smilanski Filtering
    ##- Removes chemical noise from massvectors in the massvectorlist (if \code{peptides} argument \code{TRUE}) or returns it.
    ##d Chemical noise can be removed from the peptide mass lists
    ##d due to the strong clustering of mono-isotopic peptide
    ##d peaks. Following the distance measure and filtering
    ##d method proposed by Wool Smilanski we developed an algorithm to
    ##d classify masses as peptide and non-peptide. The algorithm is based
    ##d on a modified distance measure and hierarchical clustering of all
    ##d intra massvector distances.
    ##+ object : massvectorlist.
    ##+ mdist : Minimal distance to branch to prune. The unit of the distance is Dalton.
    ##+ fraction : Maximal fraction (nr masses in branch)/(length of massvector) of branche to be prune.
    ##+ peptides : logical; \code{TRUE} - returns peptides, \code{FALSE} - returns chemical noise.
    ##+ ... : further parameters.
    ##v massvectorlist :  Returns a massvectorlist where the massvectors either contain the peptides or the non-peptides.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.
    ##sa wsiFilter.massvector, wsFilter.massvector
    ##e data(mvl)
    ##e res <- wsFilter(mvl,peptides = TRUE)
    ##e plot(res)
    ##e res2 <- wsFilter(mvl,peptides = FALSE)
    ##e plot(res2,col=2,add=TRUE)
    ##e image(res2,what="lengthmv")
    ##e hist(res)
    ##e hist(res2,col=2,add=TRUE)
    res<-lapply(object,wsFilter,mdist=mdist,fraction=fraction,peptides=peptides)
    res<- massvectorlist(experiment(object),res,project(object))
    res
  }


#----------------------------------------------
##
#calculates the spline function from the return value of getListforCalib
#it returns an spline object this spline object can be used for calibration of ppg spectra by the function
#returns a list with the spline prediction used for calibration the theoretical time and the averaged time
##

getextcalib.massvectorlist <- function(object,calib,error=250,...)
{
  ##t External Error Model
  ##- Returns the error model obtained from the calibration sample.
  ##+ object : massvectorlist with calibration sample (ppg) masses
  ##+ calib : massvector with calibration masses (ppg)
  ##+ error : relative measurment error in ppm. (peaks in this range are assumed as matching)
  ##+ ... : further parameters.
  ##v calibspline : can be used to calibrate peaklists
  ##sa calibspline, applyextcalib.massvectorlist ,getextcalib.massvector
  ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.
  ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}

  if(missing(calib))
    {
      calib <- getPPGmasses()
    }
                                        #ppg<-getPPGmasses()
  ppg<-calib
  theo<-NULL
  expd<-NULL
  for(x in 1:length(object))
     {

      exp <- object[[x]]      
      exp <- correctinternal(exp,calib,error=error,uniq=TRUE)
      match <- getaccC(mass(exp),mass(calib),error=error,uniq=TRUE)
      expd <- c(expd,exp[match$plind,1])
      theo <- c(theo,ppg[match$calind,1])
    }
  mord <- order(theo)
  expd <- expd[mord]
  theo <- theo[mord]
  require(modreg)
                                        #calculate mass dependent error
  error <- (expd-theo)*1e6/theo
  ispl <- smooth.spline(theo,error)
  ispl <- calibspline(ispl,object,error,theo)
  return(ispl)
}



applyextcalib.massvectorlist <- function(object,cS,...)
  {
    ##t External Calbiration
    ##- Corrects the massvectorlist for the measurment error stored calibspline object.
    ##d In case of external calibration some sample spots are only dedicated
    ##d to calibration. Calibration samples which produces equidistant
    ##d peaks, which exact masses are known, can be used to precisely
    ##d estimate the mass dependent error function.
    ##+ object : massvectorlist
    ##+ cS : object of class calibspline
    ##+ ... : further parameters
    ##v massvectorlist : calibrated massvectorlist
    ##sa applyextcalib.massvector, applyextcalib.massvectorlist, getextcalib.massvector, getextcalib.massvectorlist
    ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(ppg)
    ##e data(mvl)
    ##e res <- getextcalib(ppg,getPPGmasses(),error=150)
    ##e mvl2 <- applyextcalib(mvl,res)
    
    if(!inherits(cS,"calibspline"))
      {
        print("second param should be a calibspline")
        return()
      }
    res<-lapply(object,applyextcalib,cS)
    res<- massvectorlist(experiment(object),res,project(object))
    res
  }


calibexternal.massvectorlist <- function(object,ppg,calib,error=300,...)
{
  ##t External Calbiration
  ##- Perfroms external calibration of massvectorlist. Obtains the model of the measurment error by \code{getextcalib} and corrects the masses for this error.
  ##d In case of external calibration some sample spots are only dedicated
  ##d to calibration. Calibration samples which produces equidistant
  ##d peaks, which exact masses are known, can be used to precisely
  ##d estimate the mass dependent error function.
  ##+ object : massvectorlist
  ##+ ppg : either a massvector or massvectorlist with masses of the calibration sample (e.g poly-(propylene glycol) ppg)
  ##+ calib : a massvector with exact masses for the massvectors of the calibration samples (ppg).
  ##+ error : relative error of the measurment in (ppm).
  ##+ ... : further parameters.
  ##v massvectorlist: calibrated massvector
  ##sa calibexternal.massvector, applyextcalib.massvectorlist, getextcalib.massvector, getextcalib.massvectorlist, calibexternal.massvector, applycalib.calibspline
  ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. Anal Chem. 74(15):3915-23.
  ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
  ##e data(mvl)
  ##e data(ppg)
  ##e res <- calibexternal(mvl,ppg,getPPGmasses(),error=150)
  cS <- getextcalib(ppg,calib,error=error)
  object <- applyextcalib(object,cS)
  object
}

##########################################################
##Recalibration

recalibrate.massvectorlist <- function(object,PQM=7,...)
  {
    ##t Precalibration
    ##d Precalibration method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the massesvector can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##sa recalibrate.massvector
    ##- Obtains the error and performs the calibration in one step.
    ##+ object : massvectorlist.
    ##+ PQM : Peak Quality Measure. Indicates how well the wavelenght of the massvector was determined.
    ##+ ... : further parameters.
    ##v massvectorlist : recalibrated massvectorlist.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.
    ##r Wolski  \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e pp <- recalibrate(mvl,PQM=8)
    ##e compare(pp[[1]],mvl[[1]],error=1,ppm=FALSE)

    res<-lapply(object,recalibrate,PQM=PQM)
    res<- massvectorlist(experiment(object),res,project(object))
    res
  }


getrecalib.massvectorlist<-function(object,plot=FALSE,...)
  {
    ##t Precalibration
    ##d Precalibration method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the massesvector can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##- Obtain the error model using the wavelength analysis of the peaklist.
    ##+ object : massvectorlist
    ##+ plot: logical; \code{TRUE} - A graphic showing the \eqn{F(\omega)}{F(omega)}. default \code{FALSE}
    ##+ ... : further parameters.
    ##v caliblist : caliblist with objects of class calibrestat.
    ##sa getrecalib.massvector,  applyrecalib.massvector, massvector, calibrestat, wsFilter.massvectorlist
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e res <- getrecalib(mvl)
    ##e print(res)
    ##e summary(res)
    ##e image(res,what="Coef.Intercept")
    ##e image(res,what="Coef.Slope")
    ##e plot(res)
    ##e hist(res)
    ##e dres<-as.data.frame(res)
    ##e plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM")
    ##e #create subset.
    ##e res2<-subset(res,PQM>10)
    ##e length(res2)
    ##e plot(res2)
    ##e test<-applyrecalib(mvl, res2)
    ##e plot(test)
    res<-lapply(object,getrecalib,plot=FALSE)
    res<-caliblist("calibrelist",experiment(object),res,project=project(object))
    #res
  }



applyrecalib.massvectorlist <- function(object,calc,...)
  {
    ##t Recalibration
    ##- Corrects the massvector for model of the measurment error stored in the calibrestat object.
    ##d Precalibration method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the massesvector can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##+ object : massvectorlist.
    ##+ calc : caliblist with objects of class calibretstat.
    ##+ ... : further parameters.
    ##v massvectorlist : calibrated massvectorlist. 
    ##sa recalibrate.massvectorlist, getrecalib.massvectorlist, correctinternal.massvectorlist
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mvl)
    ##e res<-getrecalib(mvl)
    ##e summary(res)
    ##e mvl2<-applyrecalib(mvl,res)


    
    if(!inherits(calc,"calibrelist"))
      stop("Second argument have to be a calibrelist object!!! (use getrecalib)\n")
    for(x in 1:length(calc))
      {
        nami<-names(calc)[x]
        tmp <- object[[nami]]
        tmp<-applyrecalib(tmp,calc[[x]])
        object[[nami]]<-tmp
        #if(x%%10==0)
        #  cat(formatC(x,width=3)," ",sep="")
        #if(x%%100==0)
        #  cat("\n")
      }
    cat("\n")
    object
  }


#########################################
#affine calibration
#

applyintcalib.massvectorlist <- function(object,calc,...)
  {
    ##t Internal Calibration
    ##- Corrects the massvectors in the list for the error model stored in calibintstat object.
    ##d Internal calibration aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvectorlist
    ##+ calc : caliblist with objects of class calibintstat
    ##+ ... : further params
    ##v massvectorlist : calibrated massvectorlist. 
    ##sa applyintcalib.massvector,getintcalib.massvectorlist, correctinternal.massvectorlist, calibintstat, caliblist
    ##r Wolski
    ##e data(mvl)
    ##e mvl<-mvl[1:100]
    ##e data(cal)
    ##e res <- getintcalib(mvl,cal,error=300,ppm=FALSE)
    ##e mvl2<-applyintcalib(mvl,res)
    if(!inherits(calc,"calibintlist"))
      stop("Second argument have to be a calibintlist object!!! (use getintcalib)\n")
                                       #calconstants -  output of function recalibrate
    for(x in 1:length(calc))
      {
        tmp <- object[[names(calc)[x]]]
        tmp <- applyintcalib(tmp,calc[[x]])
        object[[names(calc)[x]]] <- tmp
                                        #        if(x%%10==0)
                                        #          cat(formatC(x,width=3)," ",sep="")
                                        #        if(x%%100==0)
                                        #          cat("\n")
      }
    cat("\n")
    object
  }


correctinternal.massvectorlist <- function(object,calib,error=500,ppm=TRUE,...)
  {
    ##t Internal Calibration
    ##- Determines the measurment error of the masses using \code{getintcalib.massvectorlist}.
    ##- It refines the error model and applies it to the massvector by \code{applyintcalib.massvectorlist}.
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvectorlist
    ##+ calib : massvector with calibration masses
    ##+ error : assumed measurment error.
    ##+ uniq : logical;\code{TRUE}- use only mass closest to calibration mass. \code{FALSE}- use all masses closer to the calibration mass then given error.
    ##+ ppm : logical;\code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.
    ##v massvectorlist : calibrated massvectorlist. 
    ##sa getintcalib.massvectorlist, correctinternal.massvector, getintcalib.massvector,calibintstat,caliblist
    ##r wolski
    ##e data(mvl)
    ##e data(cal)
    ##e mvl2 <- correctinternal( mvl,cal, error=500 , ppm=TRUE )
    tmp<-getintcalib(object,calib,error=error,ppm=ppm)
    object<-applyintcalib(object,tmp)
    object
  }

getintcalib.massvectorlist <- function(object,calib,error=500,ppm=TRUE, ... )
  {
    ##t Internal Calbiration
    ##- Obtains error model using massvector with known masses.
    ##d Internal calibration aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvectorlist
    ##+ calib : massvector with calibration masses
    ##+ error : assumed measurment error.
    ##+ uniq : \code{TRUE}- use only mass closest to calibration mass. \code{FALSE}- use all masses closer to the calibration mass then given error.
    ##+ ppm : \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.
    ##v calibintstat : object of class calibintstat. 
    ##sa applyintcalib.massvector, getintcalib.massvector, correctinternal.massvectorlist, calibintstat, caliblist
    ##r Wolski
    ##e data(mvl)
    ##e data(cal)
    ##e res<-getintcalib(mvl,cal,error=400,ppm=TRUE)
    ##e plot(res)
    
    res <- lapply(object,getintcalib,calib,error=error,ppm=ppm)
    res <- caliblist("calibintlist",experiment(object) , res, project=project(object))
    res
  }



##################################################
#Global calibration
#
getglobalcalib.massvectorlist<-function(object , calib , error=500,ppm=TRUE, labund=12, abund = length(object)/5 , accur = ifelse(ppm,error/2000,error) ,...)
  {
    ##t Set Based Internal Calbiration
    ##- Obtains error model using massvector with known masses.
    ##d Set based Calibration copes with the problem of missing calibration
    ##d masses. It first extracts about 15 most abundant masses of the
    ##d massvectorlist, then they are internally calibrated and used
    ##d as new calibration masses. In this fashion more massvectors
    ##d can be internally calibrated.
    ##d Internal calibration aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvectorlist
    ##+ calib : massvector with calibration masses
    ##+ error : assumed measurment error.
    ##+ uniq : \code{TRUE}- use only mass closest to calibration mass. \code{FALSE}- use all masses closer to the calibration mass then given error.
    ##+ ppm : logical; \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.
    ##+ labund : how many abundant masses to use for calibration. default=12.
    ##+ abund : default =  length(object)/5
    ##+ accur : default =  ifelse(ppm,error,error/2000)
    ##+ ... : further parameters.
    ##v calibintstat : object of class calibintstat. 
    ##sa applyintcalib.massvector, getintcalib.massvector, correctinternal.massvectorlist, calibintstat, caliblist
    ##r Wolski
    ##e data(mvl)
    ##e data(cal)
    ##e res<-getglobalcalib(mvl,cal,error=500,ppm=TRUE)
    ##e hist(res)
    if(accur<0.01) stop("error to small")
    mabund <- gamasses(object , accur=accur , abund=abund , ...)
    while(length(mabund)<15)
      {
        print("finding abundant masses")
        abund <- abund - 5
        print(abund)
        mabund <- gamasses(object , accur=accur , abund=abund , ... )
      }
    if(length(mabund)>labund)
      {
        ord<-order(mabund[,2])
        mabund <- mabund[ord,]
        mabund<- mabund[c((length(mabund)-(labund-1)):length(mabund)),]
        ord<-order(mabund[,1])
        mabund <- mabund[ord,]
      }
    abund <- correctinternal(mabund,calib,error=error,ppm=ppm)
    res <- getintcalib(object,abund,error=error,ppm=ppm)
    res
  }



globalcalib.massvectorlist<-function(object , calib , error=500, labund = 12, ppm=TRUE , abund=length(object)/5,accur = ifelse(ppm,error/2000,error) ,...)
  {
    ##t Set Based Internal Calibration
    ##- Determines the error and corrects for it.
    ##d Set based Calibration copes with the problem of missing calibration
    ##d masses. It first extracts about 15 most abundant masses of the
    ##d massvectorlist, then they are internally calibrated and used
    ##d as new calibration masses. In this fashion more massvectors
    ##d can be internally calibrated.
    ##d Internal calibration aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvectorlist
    ##+ calib : massvector with calibration masses
    ##+ error : assumed measurment error.
    ##+ labund : how many abundant masses use for calibration. (8-12 masses are sufficient).
    ##+ accur : default := ifelse(ppm,error/2000,error), used to determine abundant masses.
    ##+ abund : default  := length(object)/5
    ##+ ppm : \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.
    ##+ ... : further parameters.
    ##v massvectorlist : calibrated massvectorlist. 
    ##sa getintcalib.massvectorlist, correctinternal.massvectorlist, getintcalib.massvector,calibintstat,caliblist
    ##r wolski
    ##e data(mvl)
    ##e data(cal)
    ##e mvl2<-globalcalib(mvl,cal,error=500,labund=12)

    mabund <- gamasses(object , accur=accur , abund=abund , ...)
    print(mabund)
    while(length(mabund)<10)
      {
        abund <- abund - 10
        mabund <- gamasses(object , accur=accur , abund=abund , ... )
      }
    if(length(mabund)>labund)
      {
        ord<-order(mabund[,2])
        mabund <- mabund[ord,]
        mabund<- mabund[c((length(mabund)- (labund-1)):length(mabund)),]
      }
    print(mabund)
    abund <- correctinternal(mabund,calib,error=error)
    object <- correctinternal(object,mabund,error=error)
    object
  }

#############################################
#getting abundant masses
#

gamasses.massvectorlist <- function(object,accur=0.1,abund=50,...)
  {
    ##t Abundant masses
    ##- Determines abundant masses in a massvector.
    ##d Abundant masses are masses that occur in a large
    ##d fraction of the massvectors. Typical abundant mass
    ##d are derived from tryptic autoproteolysis products.
    ##d Abundant masses can also often be assigned to keratin
    ##d isoforms (human hair- and skin proteins).
    ##d Many of the abundant masses cannot be assigned to any protein.
    ##d Abundant masses can be used for calibration.
    ##d Removing them may increase the identification specificity.
    ##sa gamasses.massvector
    ##+ object : massvectorlist
    ##+ accur : measurment accuracy
    ##+ abund : how many times a mass have to occur to be an abundant mass.
    ##v massvector : massvector with abundant masses.
    ##r Wolski
    ##e data(mvl)
    ##e #Filtering for abundant masses.
    ##e res<-gamasses(mvl,abund=50)
    ##e plot(res)
    ##e mvFilter(mvl[[1]],res)
    ##e res2<-mvFilter(mvl,res,abundant=TRUE)
    ##e image(res2,what="lengthmv")
    ##e image(mvl,what="lengthmv")
    ##e image(image(res2,what="lengthmv")/image(mvl,what="lengthmv"))
    ##e 
    ##e hist(mvl,accur=0.3)
    ##e hist(res2,add=TRUE,col=2,accur=0.3)
    tmp<-unlist(object)
    res <- gamasses(tmp,accur=accur,abund=abund,main=main,xlab=xlab,xlim=xlim,...)
    res
  }




diffFilter.massvectorlist <- function(object,listofdiffs,higher=TRUE,error=0.05,uniq=TRUE,...)
{
  ##t Abundant Differences
  ##- Removes masses from the massvector.
  ##d Removes one of the masses contributing to a mass difference given in the list of diffs.
  ##d Can be used if a variable modification are present in the massvector but can not be considered by the identification software.
  ##d Abundant intra massvector mass differences indicate the
  ##d presence of variable modifications in the data set. This
  ##d information can be used to optimize the search strategy.
  ##sa diffFilter.massvector,getdiff.massvector, getdiff.massvectorlist
  ##+ mv : massvector
  ##+ listoffdiffs : massvector with mass differences
  ##+ higher : \code{TRUE} - remove higher mass, \code{FALSE} = remove lower mass.
  ##+ error : How much the differences can diviate from the differences given in listofdiffs
  ##v massvector : filtered massvector.
  ##r Wolski
  ##e data(mvl)
  ##e res<-getdiff(mvl,range=c(0,100))
  ##e hist(res)
  ##e res<-gamasses(res,abund=400)
  ##e test<-diffFilter(mvl,res,higher=TRUE,error=0.1,uniq=TRUE)
  ##e test
  res <- lapply(object,diffFilter, listofdiffs, higher=higher, error=error,uniq=uniq)
  res <- massvectorlist(experiment(object),res,project=project(object))
  res
}

getdiff.massvectorlist<-function(object,masst="massrea",range=c(0,100),...)
{
  ##t Massdifferences
  ##- Computes mass differences in the mv.
  ##d Removes one of the masses contributing to a mass difference given in the list of diffs.
  ##d Can be used if a variable modification are present in the massvector but can not be considered by the identification software.
  ##d Abundant intra massvector mass differences indicate the
  ##d presence of variable modifications in the data set. This
  ##d information can be used to optimize the search strategy.
  ##sa getdiff.massvector
  ##+ mv : massvector
  ##+ listoffdiffs : massvector with mass differences
  ##+ higher : \code{TRUE} - remove higher mass, \code{FALSE} - remove lower mass.
  ##+ error : How much the differences can diviate from the differences given in listofdiffs
  ##v massvector : filtered massvector.
  ##r Wolski
  ##e data(mvl)
  ##e res<-getdiff(mvl,range=c(0,100))
  ##e plot(res)
  ##e tt<-gamasses(res,abund=40)
  ##e plot(tt)
  res<-massvector(NULL)
  res<-lapply(mvl,getdiff,range=c(0,100))
  res<-massvectorlist(experiment(object),res,project=project(object))
  res<-unlist(res)
  res<-info(res,paste(info(object),"diffs",sep="_"))
  return(res)
}


writeF.massvectorlist<-function(object,path,file=experiment(object),ext="txt",...)
  {
    ##t Write massvectorlist
    ##- Write massvectorlist to File
    ##d The read and write functions for all the different peak-list formats are not provided by the package. This is because
    ##d there are oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package
    ##d which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}.
    ##+ object : massvectorlist.
    ##+ path : path to directory.
    ##+ file : file name. default experiment(object)
    ##+ ext : file extension. default txt.
    ##sa readF.massvectorlist, readF.massvector

    
    filep <- file.path(path,paste(file,".",ext,sep=""),fsep = .Platform$file.sep)
    con <- file(filep, "w")  # open an output file connection

    intern <-function(lobject,con)
      {
        firstline<-paste(">",info(lobject),":",join(mget(lobject,"tcoor"),sep=","),sep="")
        writeLines(firstline,con=con)
        write.table(lobject,file = con, append = TRUE, quote = FALSE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = FALSE, qmethod = c("escape", "double"))
      }
    res<-lapply(object,intern,con)
    close(con)
  }

readF.massvectorlist<-function(object,path,file=experiment(object),ext="txt",...)
  {
    ##t Reads massvectorlist
    ##- written to disk with \code{writeF.massvectorlist}
    ##d The read and write functions for all the different peak-list formats are not provided by the package. This is because
    ##d there are oodles of different formats. I will try to collect read-write functions for as many as possible peak-list format's in an add on package
    ##d which you can find at \url{http://www.molgen.mpg.de/~wolski/mscalib/IO/}.
    ##+ object :massvectorlist
    ##+ path : path to directory.
    ##+ file : file to read. default =experiment(object)
    ##sa readF.massvector,writeF.massvectorlist
    ##e data( mvl )
    ##e mvl
    ##e writeF( mvl , "." )
    ##e test <- readF(massvectorlist(info(mvl)),".")
    ##e test
    ##e unlink(paste(info(test),".txt",sep=""))
    
    filep <- file.path(path,paste(file,".",ext,sep="") , fsep = .Platform$file.sep)
    con<-file(filep,"r")
    mm <- readLines(con=filep,n=-1)
    starts <- grep(">",mm)
    starts <- c(starts,length(mm))
    res<-massvectorlist(file)
    setParms(res)<-list(project=project(object))
    intern<-function(xx)
      {
        object<-massvector()
        res1 <- unlist(strsplit(xx[1],":"))
        object<-info(object,sub(">","",res1[1]))
        setParms(object) <- list(tcoor= unlist(strsplit(res1[2],",")))
        intern2<-function(x){as.numeric(unlist(strsplit(x,"\t")))}
        tt<-xx[2:length(xx)]
        rxx <- t(sapply(tt,intern2))
        rxx<-rxx[order(rxx[,1]),]
        rownames(rxx)<-1:length(rxx[,1])
        object <- peaks(object,rxx)
        return(object)
      }

    for(x in 1:(length(starts)-1))
      {
        res[[x]]<-intern(mm[starts[x]:(starts[x+1]-1)])
      }
    close(con)
    return(res)
  }
