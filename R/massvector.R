

                                        #constructor
                                        #massvector
                                        #infofield
                                       #coordinates on target
massvector <- function(info,masses,tcoor)
{
  ##t Constructor
  ##- massvector extends matrix
  ##d The class massvector keeps the masses and their intensities in a matrix.
  ##d It stores additional attributes like the coordinates of the mass spectrometric sample on the support.
  ##d It also stores a identifier of the massvector.
  ##+ masses : a matrix or a array with masses. First column must contain masses.
  ##+ info : a unique identifier for the list.
  ##+ tcoor : a array of length two if integer sample support coordinates.
  ##e massvector("hello march",NULL)
  ##e massvector("hello march",1:100)
  ##e massvector("hello march",cbind(1:10,10:1))
  ##e tmp<-cbind(1:10,1:10)
  ##e colnames(tmp)<-c("mass","test")
  ##e rr<-massvector("hello march",tmp)
  ##e rr2<-massvector("hello bart",cbind(1:12,1:12))
  ##e # plot functions for massvector
  ##e plot(rr)
  ##e hist(rr)
  ##e summary(rr)
  ##e image(rr)
  ##e info(rr)
  ##e # setting new masses
  ##e mass(rr)
  ##e mass(rr,1:10)
  ##e peaks(rr,cbind(1:10,11:20))
  ##e # plotting with masses
  ##e data(mv1)
  ##e data(mv2)
  ##e plot(mv1,mv2)
  ##e image(mv1  ,mv2)
  ##e summary(mv1)
  ##e print(mv1)
  ##e hist(mv1)
  ##e plot(mv1)
  ##e image(mv1,mv2,error=199,ppm=FALSE)


  
  res<-matrix()
  if(!missing(masses))
    {
      if(is.null(masses))
        {
          masses <- matrix(nrow=0,ncol=2)
          colnames(masses)<-c("mass","C1")
        }
      else
        {
          #if not a matrix.
          if(!inherits(masses,"matrix"))
            {
              masses<-cbind(masses,rep(1,length(masses)))
              colnames(masses)<-c("mass",paste("C",2:length(masses[1,]),sep=""))
            }
          if(dim(masses)[2]<2)
            masses<-cbind(masses,rep(1,length(masses)))
                                        #order the matrix
          or<-order(masses[,1])
          for(i in 1:dim(masses)[2])
            {
              masses[,i]<-masses[,i][or]
            }
          if(length(masses[,1])>0)
            rownames(masses) <- 1:length(masses[,1])

          if(is.null(colnames(masses)))
            {
              colnames(masses) <- c("mass",paste("C",2:length(masses[1,]),sep=""))
            }
          else
            colnames(masses)[1] <- "mass"
        }
    }
  else
    {
      masses<-matrix(nrow=0,ncol=2)
      colnames(masses)<-c("mass","C1")
    }
  if(!missing(info))
    {
      attr(masses,"info")<-info
    }
  if(!missing(tcoor))
    {
      if(length(tcor)!=2){stop("There must be two coordinates!\n")}
      attr(masses,"tcoor") <- tcoor
    }
  attr(masses,"allow")<-c("info","tcoor")
  class(masses)<-c("massvector","myobj","matrix")
  return(masses)
}



summary.massvector<-function(object,...)
  {
    ##t Massvector Summaries
    ##- Generates a summary: min, max etc.
    ##+ object : massvector
    ##sa summary
    tmp <- dim(object)[2]
    res <- list( lengthmv=length(object) )
    for(i in 1:tmp)
      {
        xx<-object[,i]
        res <- c(res,list(summary(ifelse(length(xx)>0,xx,0))))
        names(res)[i+1]<-colnames(object)[i]
      }
    res
  }

as.vector.massvector<-function(x, mode="any")
  {
    ##t Vector
    ##- Casts massvector into vector. It does it by calling the summary.massvector method first and unlisting the result.
    ##+ x : massvector
    ##e data(mv1)
    ##e as.vector(mv1)
    unlist(summary(x))
  }

mass.massvector<-function(object,mas,...)
  {
    ##t Mass Access.
    ##- Access to the mass field of a massvector.
    ##+ object : massvector
    ##+ mas : A array with masses. If missing function returns masses.
    ##e data(mv1)
    ##e mass(mv1)
    ##e mass(mv1,1:10)
    
    if(missing(mas))
      return(object[,1])
    else
      {
        if(inherits(mas,"matrix"))
           stop(substitute(object) ," : has to be an array, to set a matrix use peaks instead!")
        masses<-cbind(sort(mas),rep(1,length(mas)))
        return(masses)
        rownames(masses)<-1:length(mas)
        colnames(masses)<-c("mass","C1")
        tmp <- attributes(object)
        at <- names(attributes(masses))
        for(l in 1:length(tmp))
          {
            if(! names(tmp)[l] %in% at)
              {
                attr(masses,names(tmp)[l])<-tmp[[l]]
              }
          }
        return(masses)
      }
  }


hist.massvector<-function(x,accur = 0.1,abund = 0, main=info(x) ,xlab="m/z",xlim=c(min(mass(x)),max(mass(x))),add=FALSE,col=1,...)
  {
    ##t Histograms
    ##- It either computes a normal histgram of the masses and intensities, or computes a fine graded histogram of the masses to show abundant masses.
    ##+ accur : sets the bin width of the histogramm.
    ##+ abund : draws a horizontal line at the frequency given by abund.
    ##+ xlab : sets the xlabels.
    ##+ xlim : sets the min and max value to be displayed.
    ##+ add : T-adds the histogram to an existing image.
    ##+ col : the color of the histogram.
    ##+ ... : further plotting arguments.
    ##sa hist
    ##e data(mv1)
    ##e hist(mv1,normal=TRUE)
    ##e hist(mv1,normal=FALSE)
    
                                        #attributes(dat)<-NULL
    mhist<-list(NULL)
    
                                        #assign indices of bins with high peak abundance
    mhist[[1]] <- hist(mass(x),breaks=seq(min(mass(x))-accur/2,max(mass(x))+1.5*accur,accur),plot=TRUE,main=main,xlab=xlab,xlim=xlim,add=add,col=col,border=col,...)
    mhist[[2]] <- hist(mass(x),breaks=seq(min(mass(x))-accur,max(mass(x))+accur,accur),add=TRUE,col=col,border=col)
    abline(h=abund,lty=2)
 
  }




#setting and getting masses

peaks.massvector<-function(object,masses,...)
  {
    ##t Data Access
    ##- Access to the mass field of the massvector.
    ##+ object : massvector
    ##+ masses : matrix with masses in first column. If missing the matrix of the massvector is returned.
    ##sa mass.massvector
    ##e mv1<-massvector()
    ##e mv1<- peaks(mv1,cbind(1:10,1:10))
    ##e peaks(mv1)
    if(missing(masses))
      {
        for(y in attr(object,"allow"))
          {
            attr(object,y)<-NULL
          }
        attr(object,"allow")<-NULL
        class(object)<-"matrix"
        return(object)
      }
    else
      {

        if(!inherits(masses,"matrix"))
          {
            stop("The second arg (masses) must be a matrix\n")
          }
        or<-order(masses[,1])
        for(i in 1:length(dim(masses)[2]))
          {
            masses[,i]<-masses[,i][or]
          }
        if(is.null(rownames(masses)))
          rownames(masses)<-1:length(masses[,1])
        if(is.null(colnames(masses)))
          colnames(masses)<-c("mass",paste("C",2:length(masses[1,]),sep=""))
        else
          colnames(masses)[1]<-c("mass")

        tmp <- attributes(object)
        at <- names(attributes(masses))

        for(l in 1:length(tmp))
          {
            if(! names(tmp)[l] %in% at)
              attr(masses,names(tmp)[l])<-tmp[[l]]
          }
        return(masses)
      }
  }

image.massvector<-function(x,mv2,error=NULL,ppm=FALSE,col=topo.colors(100),...)
  {
    ##t Display a Color Image
    ##- Creates a grid of colored or gray-scale rectangles with colors
    ##- corresponding to the mass differences within the peaklist
    ##- or within two peaklists.
    ##+ x : massvector
    ##+ mv2 : massvector
    ##+ error : up to which mass difference display the differences.
    ##+ col : a list of colors such as that generated by `rainbow',`heat.colors', `topo.colors', `terrain.colors' or similar functions.
    ##... : graphical parameters for `plot' may also be passed as arguments to this function.
    ##sa plot.massvector, hist.massvector
    ##e data(mv2)
    ##e data(mv1)
    ##e image(mv1,mv2)
    ##e image(mv1,mv2,error=500)
    
    if(!missing(mv2))
      {
        if(!inherits(mv2,"massvector"))
          stop("Second arg should be a massvector too!\n")
      }
    else
      {
        mv2<-x
      }
    res<-NULL
    for(u in 1:length(x[,1]))
      {
        if(ppm)
          tmp<-abs(mv2[,1]-x[u,1])/mv2[,1]*10e6
        else
          tmp<-abs(mv2[,1]-x[u,1])
        tmp[tmp>error]<-NA
        res<-rbind(res,tmp)
      }
    par(bg="gray")
    nf <- layout(matrix(c(1,2),1,2),widths=c(4,1), TRUE)
    tmar<-par()$mar
    par(mar=c(5,5,1,1))
    image(1:length(x[,1]),1:length(mv2[,1]),res,col=col,xlab=info(x),ylab=info(mv2))
    tres<-na.omit(c(res))
    if(min(tres)!=max(tres))
      {
        scale<-seq(min(tres),max(tres),(max(tres)-min(tres))/9)
      }
    else
      {
        scale<-rep(min(tres),10)
      }
    scale<-matrix(scale,nrow=1)
    par(mar=c(5,1,1,1))
    image(1,1:10,scale,axes=FALSE,xlab="",ylab="",col=col)
    scalet<-format(scale,digits=1)
    for(u in 1:length(scalet))
      {
        text(25,u,scalet[u])
      }
    layout(matrix(1))
    par(mar=tmar)
    invisible(res)
  }

min.massvector<-function(mv,...)
  {
    ##t Maxima and Minima
    ##- Returns the maxima and minima of the masses and intensities
    ##sa max.massvector
    ##+ mv : massvector
    ##v named array with minima of the columns of mv.
    return(apply(mv,2,min))
  }

max.massvector<-function(mv,...)
  {
    ##t Maxima and Minima
    ##- Returns the maxima and minima of the masses and intensities
    ##sa min.massvector
    ##+ mv : massvector
    ##v named array with maxima of the columns of mv.
    return(apply(mv,2,max))
  }

length.massvector<-function(x,...)
  {
    ##t Length of Massvector
    ##- get the Length of a Massvector
    ##+ x : massvector
    ##e data(mv1)
    ##e length(mv1)
    
    return(length(x[,1]))
  }


plot.massvector<-function(x,...)
  {
    ##t Massvector Plotting
    ##- Function for plotting massvectors. If one massvector are given it shows a stick masspectrum.
    ##- If two massvectors are given their masses are plotted against each other.
    ##+ x : massvector
    ##+ ... : a second massvector and graphical parameters can be given as arguments to `plot'.
    ##e data(mv1)
    ##e data(mv2)
    ##e plot(mv1)
    ##e plot(mv1,mv2)
    mv<-x
    rm(x)
    if(length(mv)==0)
      return(paste(as.character(substitute(x))," has length 0!",sep=""))
    pars<-list(...)
    tmp<-FALSE
    if(length(pars)>0)
      {
        tmp<- inherits(pars[[1]],"massvector")
      }
    if(tmp)
      {
        plot(1,1
             ,xlim = if(is.null(pars[["xlim"]])){c(min(mv[,1]),max(mv[,1]))}else{pars[["xlim"]]}
             ,ylim = if(is.null(pars[["ylim"]])){c(min(pars[[1]][,1]),max(pars[[1]][,1]))}else{pars[["ylim"]]}
             ,ylab = info(pars[[1]]),xlab=info(mv))
        abline(v=mv[,1])
        abline(h=pars[[1]][,1])
        if(!is.null(pars[["error"]]))
          {
            abline(coef=c(0,1),col=2)
            abline(coef=c(pars[["error"]],1),col=3)
            abline(coef=c(-pars[["error"]],1),col=3)
          }
        else
          {
            abline(coef=c(0,1),col=2)
            abline(coef=c(1,1),col=3)
            abline(coef=c(-1,1),col=3)
          }
      }
    else
      {
        parsn<-names(pars)
        if("main" %in% parsn)
          {
            m <- pars[["main"]]
            pars <- subset(pars,select=-main)
          }
        else
          m <- info(mv)
        x<-mv[,1]
        y<-mv[,2]
        ylim<-c(0,max(y))
        xlim<-c(min(x),max(x))
                                        #        xlim<-c(min(x),max(x))
        sw<-FALSE
        if(!is.null(pars[["add"]]))
          sw<-pars[["add"]]
        
        if("add" %in% parsn)
          {
            pars <- as.list(subset(as.data.frame(pars),TRUE,select=-add))
          }
        parj<-join(paste(names(pars),pars,sep="="),sep=",")
        if(is.na(parj))
          {
            parj<-""
          }
        if(sw)
          {
            test <- parse(text=paste("points(x,y,type=\"h\",",parj,")",sep="\n"))
          }
        else
          {
            test <- parse(text=paste("plot.default(x,y,type=\"h\",main=\"",m,"\",axes=TRUE,xlab=\"m/z\",ylab=\"",colnames(mv)[2],"\",ylim=\c(",join(ylim,sep=","),"),xlim=\c(",join(xlim,sep=","),"),",parj ,")",sep="") )
          }
        dataf <- eval(test)
        abline(h=0)
      }
  } 


                                        #filtering
mvFilter.massvector<-function(object,fmass,match=FALSE,error=250,ppm=TRUE,uniq=FALSE,...)
  {
    ##t Filtering Massvector
    ##- Filters massvector for masses given in a second massvector.
    ##+ object : massvector
    ##+ abundm : massvector
    ##+ error : mesurment error.
    ##+ ppm : given either in ppm (\code{TRUE}) or as absolut error (F).
    ##+ match: logical; \code{TRUE} - than returns masses matching to the masses in massvector abundant, \code{FALSE} - returns masses not matching.
    ##+ uniq : logical; \code{FALSE} - returns all masses in the range given by error. \code{TRUE} - returns only the closest mass.
    ##v massvector : with matching or not matchin masses.
    ##e data(mv1)
    ##e data(mv2)
    ##e mvFilter(mv1,mv2,error=250,match=FALSE)
    ##e mvFilter(mv1,mv2,error=250,match=TRUE)

    mmatch<-getaccC(mass(object),mass(fmass),error=error,ppm=ppm,uniq=uniq)
    if(length(mmatch$plind)>0)
      {
        if(!match)
          {
            return(object[-mmatch$plind,])
          }
        else
          {
           return(object[mmatch$plind,])
          }
      }
    else
      {
        if(!match)
          return(object)
        else
          return(object[NULL,])
      }
  }

print.massvector <- function(x,quote=FALSE,...)
  {
    ##t Print massvector
    ##- `print' prints its argument and returns it invisibly (via `invisible(x)')
    ##+ x : massvector
    ##sa \link[base]{print}
    ##e data(mv1)
    ##e print(mv1)
        
    print(paste("info   : (",attr(x,"info"),")",sep=""),quote=quote,...)
    inf<-info(x)
    mcor<-attr(x,"tcoor")
    ncor<-names(mcor)
    if(!is.null(ncor))
      print(paste("coor   :  ",ncor[1],"=",mcor[1]," ; ",ncor[2],"=",mcor[2] ,sep=""),quote=quote,...)
    else
      print(paste("coor   :  ",mcor[1]," ; ",mcor[2] ,sep=""),quote=quote,...)
    attr(x,"allow") <- NULL
    attr(x,"info") <- NULL
    attr(x,"tcoor") <- NULL
    class(x) <- NULL
    print(x,...)
    invisible(x)
  }


compare.massvector <- function(object,mv2,plot=TRUE,error=1000,ppm=TRUE,...)
  {
    ##t Compares massvectors
    ##- Compares the masses in the massvectors. Returns basic statistics about the matching peaks.
    ##- Plots the relative or absolute error of matchin peaks.
    ##+ object : massvector
    ##+ mv2 : massvector
    ##+ plot : True - plot the relatvetor absolute error. \code{FALSE} - no plotting.
    ##+ error : size of the measurment error (default 150 ppm)
    ##+ ppm : \code{TRUE} - relative error in parts per million, \code{FALSE} - absolute error.
    ##+ ... : further parameters.
    ##v FMSTAT : Fowlkes & Mallows statistik (nr matching)/sqrt(length(object)*length(mv2))
    ##v min : smallest error
    ##v ... : 1st qu. , mean, median, 3rd qu., max and stdv of error.
    ##e data(mv1)
    ##e data(mv2)
    ##e compare(mv1,mv2,error=5000,ppm=TRUE,uniq=TRUE)
    ##e compare(mv2,mv1,error=1,ppm=FALSE,uniq=TRUE)
    
    if(!inherits(mv2,"massvector"))
      stop("Second arg are not a massvector!\n")
    match<-getaccC(mass(object) , mass(mv2) , error=error , ppm=ppm , uniq=TRUE)
    print(match)
    if(ppm)
      {
        if(plot)
          plot(object[match$plind,1],(object[match$plind,1]-mv2[match$calind,1])*1e6/mv2[match$calind,1],main=paste(info(object),info(mv2)),xlab="mass",ylab="error[ppm]")
        res <- (object[match$plind,1]-mv2[match$calind,1])*1e6/mv2[match$calind,1]
      }
    else
      {
        if(plot)
          plot(object[match$plind,1],(object[match$plind,1]-mv2[match$calind,1]),main=paste(info(object),info(mv2)),xlab="mass",ylab="Da")
        res <- (object[match$plind,1]-mv2[match$calind,1])
        
      }
    if(length(res)>0)
      {
        res<-c(length(match$plind)/sqrt(length(object)*length(mv2)),min(res),quantile(res,0.25),mean(res),median(res),quantile(res,0.75),max(res),sqrt(var(res)))
      }
    else
      {
         res<-c(length(match$plind)/sqrt(length(object)*length(mv2)),NA,NA,NA,NA,NA,NA,NA)
      }
    names(res)<-c("FMSTAT","min","1st qu.","mean","median","3rd qu.","max","stdv")
    invisible(res)
  }

mget.massvector<-function(object,attrn,...)
  {
    ##t Field Access
    ##- Access to the fields in the massvector
    ##+ object : massvector
    ##+ attrn : The value of which field to return. If missing the fields of the objects are returned.
    ##+ ... : further parameters
    ##v xxx : depends which field in the massvector are accessed.
    ##e data(mv1)
    ##e mget(mv1,"info")
    ##e mget(mv1,"peaks")
    
    if(missing(attrn))
      {
        return(attr(object,"allow"))
      }
    else
      {
        if(attrn %in% "peaks")
          {
            al<-attr( object , "allow" )
            for(x in al)
              attr(object,x) <- NULL
            attr(object,"allow")<- NULL
            class(object)<-"matrix"
            return(object)
          }
        else
          NextMethod("mget")
      }
  }


as.matrix.massvector <- function(x)
  {
    ##t Matrices
    ##- Turns the massvector into a matrix
    ##+ x : massvector
    ##v matrix : matrix with masses and intensities.
    ##sa peaks.massvector
    ##e data(mv1)
    ##e res<-as.matrix(mv1)
    ##e class(res)
    return(mget(x,"peaks"))
  }

"[.massvector"<-function(peak,i,j)
  {
    ##t Extract Parts of an Massvector
    ##- The massvector extends matrix. The masses and intensities can be accessed like in case of a matrix.
    ##- Do not use mv[1:10] (not handled properly)
    ##+ peak : object from which to extract elements.
    ##+ i : row access
    ##+ j : column access
    ##v xxx : If rows are selected then a massvector is returned. If columns are accessed than arrays are returned.
    ##sa peaks.massvector, mass.massvector, mget.massvector
    ##e data(mv1)
    ##e ls()
    ##e mv1[1:10,] # returns a massvector of length 10
    ##e mv1[1:10,1] # the first ten masses are returned.
    ##e mv1[,1] # the masses are returned.
    ##e mv1[,2] # the peak area are returned.
 
    
    res<-NextMethod("[")
    if(missing(j) & !missing(i))
                                        #     if(inherits(res,"matrix"))
      {
        if(!inherits(res,"matrix"))
          {
            tt<-matrix(res,ncol=dim(peak)[2])
            colnames(tt)<-colnames(peak)
            rownames(tt)<-1:length(tt[,1])
            res<-tt
          }
        at<-names(attributes(res))
        atr<-attributes(peak)
        for(x in 1:length(atr))
          {
            cur <- names(atr)[x]
            if(!cur %in% at)
              attr(res,cur) <- atr[[x]]
          }
      }
    res
  }

c.massvector<-function(x,...)
  {
    ##t Combine Massvectors into one Massvector.
    ##- Combines Massvectors into one Massvector.
    ##+ x : massvector
    ##+ ... : massvectors to be concatenated.
    ##sa rbind
    ##v massvector : massvector
    ##e data(mv1)
    ##e data(mv2)
    ##e par(mfrow=c(2,1))
    ##e plot(mv1)
    ##e plot(mv2,add=TRUE)
    ##e plot(c(mv1,mv2))
    
    tmp<-list(...)
    for(y in tmp)
      {
        if(inherits(y,"massvector"))
          x <- massvector(paste(info(x),info(y),sep="_"),rbind(x,y))
      }
  }
                                        #------------------------------------------------





getPPGmasses <- function(start=10,end=100)
{
  ##t PPG masses
  ##- Computes poly-(propylene glycol) masses.
  ##+ start : length of shortest polymer.
  ##+ end : length of longest polymer.
  ##v massvector : massvector of ppg masses
  ##e plot(getPPGmasses(start=12,end=100))
  mC  <-  12
  mH <- 1.007825
  mO <- 15.994915
  mNa <- 22.98977
  me <- 0.00054858
  massPPG <- mC*3 + mH*6 + mO
  massOH <- mO+mH

  n <- start:end
  n <- massPPG * n + massOH + mH + mNa - me
  names(n) <- start:end
  tt <-cbind(n,start:end)
  colnames(tt)<-c("mass","n")
  massvector("theoretical ppg masses",tt)
}

##################################################
#
# Wool Smilanski filtering
#

wsiFilter.massvector<-function(object , mdist=0.25 , fraction=0.2 , ... )
  {
    ##t Smilanski Filtering
    ##- The function returns the inidces of masses identified as chemical noise.
    ##d Chemical noise can be removed from the peptide mass lists
    ##d due to the strong clustering of mono-isotopic peptide
    ##d peaks. Following the distance measure and filtering
    ##d method proposed by Wool Smilanski we developed an algorithm to
    ##d classify masses as peptide and non-peptide. The algorithm is based
    ##d on a modified distance measure and hierarchical clustering of all
    ##d intra massvector distances.
    ##+ object : massvector.
    ##+ mdist : minimal distance to branch to be prune. The unit of this distance are Daltons.
    ##+ fraction : maximal fraction (nr masses in branch)/(length of massvector) of branche  to be prune.
    ##+ ... : further arguments.
    ##v indices :  Indices of masses which are identified as being nonpeptide.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. \emph{Proteomics.} 2(10):1365-73.
    ##r Wolski
    ##sa wsFilter.massvector,wsdist.massvector, wsFilter.massvectorlist
    ##e data(mv1)
    ##e data(mv2)
    ##e length(mv1)
    ##e length(wsFilter(mv1))
    ##e length(mv2)
    ##e length(wsFilter(mv2))
    
    if(length(object)<2)
      return(NULL)
                                        #die function gibt die indices der contaminanten zurück.
                                        #der aufrufenden function bleibt überlassen was sie mit ihnen macht.
                                        #pl - peaklist
                                        #mdist - distance of the contaminant cluster.
                                        #fraction - the contaminant cluster shouldnt be greater than fraction.
    if(is.null(rownames(object)))
      rownames(object)<-1:length(object)
    hst <- hclust(wsdist(object),method="single")
    tmp <- which(hst$height > mdist)
    if(length(tmp)>0)
      {
         tmp<-cutree(hst,2)
                                        #count how big are the branches.
         c2<-sum(tmp==2)
         c1<-sum(tmp==1)
         if(min(c(c2,c1))/max(c(c2,c1))< fraction)
           {
             if(c2>c1)
               {
                 res<- as.numeric(rownames(object)[tmp==1])
                 res<-c( res , wsiFilter(object[tmp==2,] , mdist=mdist , fraction=fraction))
                 return(as.numeric(res))
               }else{
                 res<-as.numeric(rownames(object)[tmp==2])
                 res<-c(res , wsiFilter(object[tmp==1,] , mdist=mdist , fraction=fraction))
                 return(res)
               }
           }
      }
                                        # if no contaminations where found.
                                        #wenn keine contaminationen gefunden wruden.
    return(NULL)
  }

wsFilter.massvector <- function(object,mdist=0.25,fraction=0.2, peptides=TRUE,...)
  {
    ##t Smilanski Filtering
    ##- Removes chemical noise from the massvectorlist.
    ##d Chemical noise can be removed from the peptide mass lists
    ##d due to the strong clustering of mono-isotopic peptide
    ##d peaks. Following the distance measure and filtering
    ##d method proposed by Wool Smilanski we developed an algorithm to
    ##d classify masses as peptide and non-peptide. The algorithm is based
    ##d on a modified distance measure and hierarchical clustering of all
    ##d intra massvector distances.
    ##+ object : massvector.
    ##+ mdist : minimal distance of branch to be cut.
    ##+ fraction : maximal size of branche (nr masses in branch)/(length of massvector) to be cut.
    ##+ peptides : \code{TRUE} - returns peptides, \code{FALSE} - returns chemical noise.
    ##+ ... : further parameters.
    ##v massvector :  returns either a massvector of nonpeptide masses or massvector of peptide masses.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. Proteomics. 2(10):1365-73.
    ##sa wsiFilter.massvector,wsFilter.massvectorlist
    ##e data(mv1)
    ##e tmp <- wsFilter(mv1,peptide=FALSE)
    ##e plot(tmp)
    ##e tmp <- wsFilter(mv1,peptide=TRUE)
    ##e plot(tmp)
    
    
    tmp<-wsiFilter(object,mdist=mdist,fraction=fraction)
    if(peptides)
      {
        if(length(tmp)==0)
          return(object)
        else
          return(object[-tmp,])
      }
    else
      {
        return(object[tmp,])
      }
  }




wsdist.massvector <- function(object,...)
  {
    ##t Wool Smilanski Distance Matrix
    ##- This function computes and returns the distance matrix
    ##- using the intra massvecotor distance  measure.
    ##+ object : massvector
    ##v dist : an object of class distance.
    ##sa wsFilter.massvector, wsiFilter.massvector,
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. \emph{Proteomics.} 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e plot(hclust(wsdist(mv1),method="single"))
    ret<-NULL
    pl<-object[,1]
    for(x in pl)
    {
      ret<-rbind(ret,distance(x,pl))
    }
    rownames(ret) <- formatC(pl,digits=2,format="f")
    colnames(ret) <- format(pl,digits=2,format="f")
    ret <- as.dist(ret)
    return(ret)
  }


###########################################
# External calibration


  
getextcalib.massvector <- function(object,calib,error=300,...)
{
  ##t External Error Model
  ##- Returns the error model obtained from the calibration sample.
  ##d In case of \bold{external calibration} some sample spots are only dedicated
  ##d to calibration. Calibration samples which produces equidistant
  ##d peaks, which exact masses are known, can be used to precisely
  ##d estimate the mass dependent error function.
  ##+ object : massvector
  ##+ calib : massvector with calibration masses
  ##+ error : relative measurment error in ppm.
  ##v calibspline : can be used to calibrate peaklists
  ##sa calibspline, applyextcalib.massvector
  ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. \emph{Anal Chem.} 74(15):3915-23.
  ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
  ##e data(mv1)
  ##e data(ppg)
  ##e res<-getextcalib(ppg[[1]])
  ##e mv2 <- applycalib(res,mv1)
  ##e plot(mv1[,1],mv2[,1]-mv1[,1])

  if(missing(calib))
    {
      calib<-getPPGmasses()
    }
  object <- correctinternal(object,calib,error=error)
  theo<-NULL
  expd<-NULL
  exp<-object
  match<-getaccC(mass(exp),mass(calib),error=250,uniq=TRUE)
  expd <- c(expd,exp[match$plind])
  theo <- c(theo,calib[match$calind])

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


applyextcalib.massvector <- function(object,cS,...)
{
  ##t External Calbiration
  ##- Applys object of class calibspline to massvector to correct ofr measurment errors.
  ##- The error model stored in the calibspline are obtained by the function \code{getextcalib}
  ##d In case of \bold{external calibration} some sample spots are only dedicated
  ##d to calibration. Calibration samples which produces equidistant
  ##d peaks, which exact masses are known, can be used to precisely
  ##d estimate the mass dependent error function.
  ##+ object : massvector
  ##+ cS : calibspline
  ##v massvector : calibrated massvector
  ##sa applyextcalib.massvectorlist, applycalib.calibspline,getextcalib.massvector, getextcalib.massvectorlist
  ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS. \emph{Anal Chem.} 74(15):3915-23.
  ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
  ##e data(ppg)
  ##e data(mv1)
  ##e res<-getextcalib(ppg[[1]],getPPGmasses(),error=200)
  ##e mv2<-applyextcalib(mv1,res)
  ##e compare(mv1,mv2)

  if(!inherits(cS,"calibspline"))
    {
      stop(as.character(substitute(cS)), " should be a calibspline but is of class: ", class(cS),"\n" )
    }
  #peaklist- array of peaks masses.
  #spline - spline to predict the error.
  error <- predict(cS,mass(object))
  error <- error$y

  masspred <- object[,1]/(1+error/1e6)
  object[,1] <- masspred
  return(object)
}


calibexternal.massvector <- function(object,ppg,calib,...)
  {
    ##t External Calbiration
    ##- Perfroms external calibration. Obtains the error model by calling \code{getextcalib}
    ##- and corrects the masses in the massvector for the errror.
    ##d In case of \bold{external calibration} some sample spots are only dedicated
    ##d to calibration. Calibration samples which produces equidistant
    ##d peaks, which exact masses are known, can be used to precisely
    ##d estimate the mass dependent error function.
    ##+ object : massvector
    ##+ ppg : either a massvector or massvectorlist with masses of the calibration sample (e.g. poly-(propylene glycol)).
    ##+ calib : a massvector with exact masses of the calibration sample (ppg).
    ##+ ... : further parameters.
    ##v massvector: calibrated massvector
    ##sa applyextcalib.massvectorlist, getextcalib.massvector, getextcalib.massvectorlist, calibexternal.massvector, applycalib.calibspline
    ##r Gobom J, Mueller M, Egelhofer V, Theiss D, Lehrach H, Nordhoff E, 2002. A calibration method that simplifies and improves accurate determination of peptide molecular masses by MALDI-TOF MS.\emph{Anal Chem.} 74(15):3915-23.
    ##r Wolski
    ##e data(ppg)
    ##e data(mv1)
    ##e mv2<-calibexternal(mv1,ppg)
    ##e compare(mv1,mv2)
    
    if(inherits(ppg,"massvector") | inherits(ppg,"massvectorlist"))
      {}
    else
      {
        warning("The second param should be a massvector or massvectorlist!\n")
        return(FALSE)
      }
    if(missing(calib))
      {
        calib<-getPPGmasses()
      }
    cS <- getextcalib(ppg,calib)
    res <- applyextcalib(object,cS)
    res
  }

#recalibration



recalibrate.massvector <- function(object,PQM=7,...)
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
    ##- Obtains the error and performs the calibration in one step.
    ##+ object : massvector
    ##+ PQM : Peak Quality Measure. Indicates how well the wavelenght of the massvector was determined.
    ##+ ... : further parameters
    ##w massvector : recalibrated massvector.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. {\em Proteomics.} 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e mv2<-recalibrate(mv1)
    ##e plot(mv1[,1],mv1[,1]-mv2[,1],type="l")

    cal<-getrecalib(object)
    if(mget(cal,"PQM")>PQM)
      object<-applyrecalib(object,cal)
    return(object)
  }


applyrecalib.massvector<-function(object,calc,...)
  {
    ##t Precalibration
    ##d \bold{Precalibration} method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the massesvector can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##- Uses the model of the measurment error stored in a calibrestat object
    ##- to correct the masses in the massvector.
    ##+ object : massvector.
    ##+ calc: object calibrestat class.
    ##+ ...: further arguments.
    ##w massvector : recalibrated massvector.
    ##sa getrecalib.massvector, getrecalib.massvectorlist, recalibrate.massvector, recalibrate.massvectorlist, calibrestat, calibrelist
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting.\emph{Proteomics.} 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e plot(res)
    ##e mv2<-applyrecalib(mv1,res)
    ##e compare(mv1,mv2)

    
    peak<- object[,1]*(mget(calc,"Coeff.Slope")/1e6+1) + mget(calc,"Coeff.Intercept")
    object[,1] <- peak
    return(object)
  }

getrecalib.massvector<-function(object,plot=FALSE,...)
  {
    ##t Precalibration
    ##d \bold{Precalibration} method utilizes the knowledge that masses
    ##d of peptides are in equidistant spaced clusters. The wavelength of
    ##d the massesvector can be determined as described by
    ##d Wool. The comparision of the experimental wavelength with
    ##d the theoretical one, makes possible to find an affine function
    ##d that corrects the masses. Chemical noise in the spectra may hamper
    ##d the determination of mass list frequency. The package provides a
    ##d function to filter chemical noise.
    ##- Obtains the an error model using the wavelength analysis of the peaklist.
    ##+ mv : massvector
    ##+ plot: \code{TRUE} - A graphic showing the \eqn{F(\omega)}{F(omega)} function. default \code{FALSE}
    ##v calibrestat : object of class calibrestat.
    ##sa applyrecalib.massvector, massvector, calibrestat, wsFilter.massvector
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting.\emph{ Proteomics.} 2(10):1365-73.
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e print(res)
    ##e as.vector(res)
    ##e summary(res)
    ##e image(res)
    ##e plot(res)


    mv<-mass(object)
                                        #returns calibration constants obtained by the wools smilanski method.
    if(length(object)<3)
      {
        res<-calibrestat(info(object))
        setParms(res)<-list(Coeff.Intercept=0)
        setParms(res)<-list(Coeff.Slope=0)
        setParms(res)<-list(lengthmv = length(object))
        setParms(res)<-list(tcoor=mget(object,"tcoor"))
        setParms(res) <- list(PQM = 0)
        return(res)
      }
                                        #   lambda<-seq(1,1.001,0.000001)
    lambda<-seq(0.999495,1.001495,0.000005)
    omega<-2*pi/lambda
    test<-mv%*%t(omega)
    testsin<-sin(test)
    testcos<-cos(test)
    colsumsin <-apply(testsin,2,sum)
    colsumcos <-apply(testcos,2,sum)
    sumcol<-sqrt(colsumsin^2+colsumcos^2)
    sumcol <- cbind(lambda,sumcol)
                                        #determine the maximum and get its index
    mmax <- max(sumcol[,2])
    index<-sumcol[,2] == mmax
    if(FALSE)
      {
        return(sumcol)
      }
                                        #temporary only for testing purposes.
    if(plot){
      len<-length(mv)
      main <- paste("Massvector length : ", len, sep="")
      plot(sumcol,type="l",ylab="Amplitude",xlab="Wavelength",xlim=c(0.9995,1.0015),las=2, main = main)
      points(sumcol[,1][index],sumcol[,2][index],col=2,pch="*")
      stat<-FullWidthatHalfMaximum(sumcol)
      abline(h=stat$hm)
      abline(v=stat$width)
                                        #end temporary
    }
                                        #by which wavelength the maximum occure
    lambdamax<-sumcol[,1][index]
                                        #calculate the phase shift.
    phimax<-atan(colsumsin[index]/colsumcos[index])
                                        #we now see that the peak centers lie on the line
                                        # M= lambdamax*N + bmaxx
    bmax<-lambdamax*phimax/(2*pi)
                                        #to coorect the peaklist apply
                                        #peak<-(peaklist-bzero)/alpha
    alpha <- 1.000495/lambdamax
    res<-calibrestat(info(object))
    tmp<-c(bmax,alpha)
    names(tmp)<- c("Intercept","Slope")
    setParms(res)<-list(Coeff.Intercept=- (tmp[1]*tmp[2]))
    setParms(res)<-list(Coeff.Slope=(tmp[2]*1e6-1e6))
    setParms(res)<-list(lengthmv = length(object))
    setParms(res)<-list(tcoor=mget(object,"tcoor"))
    tmp <- FullWidthatHalfMaximum(sumcol)
    setParms(res) <- list(PQM = tmp$PQM)
                                        #    setParms(res)<-list(stat=tmp)
    return(res)
  }


##############################################
## affine calibration
#


applyintcalib.massvector <- function(object,cal,...)
  {
    ##t Internal Calibration
    ##- Corrects the massvector for the error model stored in calibintstat object.
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvector
    ##+ cal : object of class calibintstat
    ##+ ... : further parameters.
    ##v massvector : calibrated massvector. 
    ##sa applycalib.calibintstat,getintcalib.massvector, correctinternal.massvector
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e data(cal)
    ##e res<-getintcalib(mv1,cal,error=300)
    ##e mv2<- applyintcalib(mv1,res)
    ##e plot(mv1[,1],mv2[,1]-mv1[,1])
    if(length(object)==0)
      {
        return(object)
      }
    if(inherits(cal,"lm"))
      {
        errp <- predict(cal,data.frame(masstheo=mass(object) ) )
        if(mget(cal,"ppm"))
          {
            object[,1] <- object[,1]/(1 - errp/1e6)
          }
        else
          {
            object[,1] <- object[,1] + errp
          }
      }
    object
  }


correctinternal.massvector<-function(object,calib,error=500,uniq=FALSE,ppm=TRUE,...)
  {
    ##t Internal Calibration
    ##- Corrects the masses of the massvector. It first obtains the model of the
    ##- measrument error by calling \code{getintcalib.massvector}. It than corrects the masses
    ##- by a call to \code{applyintcalib.massvector}.
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvector
    ##+ calib : massvector with calibration masses
    ##v massvector : calibrated massvector. 
    ##sa getintcalib.massvector, calibintstat, applycalib.calibintstat
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
    ##e data(mv1)
    ##e data(cal)
    ##e res <- correctinternal(mv1,cal,error=200)
    
    tmp <- getintcalib(object,calib,error=error,uniq=uniq,ppm=ppm)
    tmp <- applyintcalib(object,tmp,ppm=ppm)
    tmp
  }


getintcalib.massvector <- function(object,calib,error=500,uniq=FALSE,ppm=TRUE,...)
  {
    ##t Internal Calbiration.
    ##- Obtains error model by alingning masses in massvector to known masses (calibration list).
    ##d \bold{Internal calibration} aligns masses of
    ##d peaks to known masses and determines by linear regression a affine
    ##d function that describing the relative error. The internal
    ##d correction fails when no calibration peaks can be found.
    ##+ object : massvector
    ##+ calib : massvector with calibration masses
    ##+ error : assumed measurment error.
    ##+ uniq : \code{TRUE}- use only mass closest to calibration mass. \code{FALSE}- use all masses closer to the calibration mass then given error.
    ##+ ppm : \code{TRUE}- describe the error as relative error. \code{FALSE}- describe the error as absolute error.
    ##v calibintstat : object of class calibintstat. 
    ##sa applyintcalib.massvector, getintcalib.massvector, correctinternal.massvectorlist, calibintstat
    ##r Wolski \url{http://www.molgen.mpg.de/~wolski/mscalib}
                                        #peaklist - peaklist.
                                        #calib - list with calibrants.
    massmv<-mass(object)                    #use mass(.. because you still dont shure about implementation.
    masscalib<-mass(calib)
    match<-NULL
    match<-getaccC(massmv,masscalib,error=error,ppm=ppm,uniq=uniq)
    if(length(match$plind)>1)
      {
                                        #calculate mass dependent error function.
                                        #get peaklist matching
        smallerrpl<-massmv[match$plind]
                                        #get calibrants matching
        masstheo<-masscalib[match$calind]
        if(ppm)
          {
            err  <-  (masstheo - smallerrpl) * 1e6 /masstheo
          }
        else
          {
            err  <-  masstheo - smallerrpl
          }
                                        #finally calculate error
                                        #distinguishing 2 cases. if peaks quite close to each other only correct for offset.
        if(abs(diff(range(masstheo)))<200)
          {
            err<-mean(err)
            masstheo <- mean(masstheo)
            mymod <- lm(err~masstheo)
          }else{
            mymod <- lm(err~masstheo)
          }
        coraf <- calibintstat(info(object),mymod)
        stat <- c(mean(abs(err)),sqrt(var(err)))
        names(stat) <- c("mean","stdv")
        setParms(coraf) <- list(Coeff.Intercept= mymod$coefficients[1],
                                Coeff.Slope =  ifelse(ppm,mymod$coefficients[2]*1e4,mymod$coefficients[2]*1e6) ,
                                lengthmv = length(object),
                                nrmatch = length(match$plind),
                                ppm=ppm,
                                tcoor = mget(object,"tcoor")
                                )
      }
    else
      {
        err <- rep(0,10)
        masstheo <- 1:10
        mymod <- lm(err~masstheo)
        coraf <- calibintstat(info(object),mymod) 
        setParms(coraf) <- list( Coeff.Intercept = mymod$coefficients[1] ,
                                Coeff.Slope =  ifelse(ppm,mymod$coefficients[2]*1e4,mymod$coefficients[2]*1e6),
                                lengthmv = length(object),
                                nrmatch = 0 ,
                                ppm=ppm,
                                tcoor = mget(object,"tcoor")
                               )
      }
    return(coraf)
  }


################################################
# gamasses
#


gamasses.massvector <- function(object,accur = 0.1,abund = 50,...)
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
    ##sa gamasses.massvectorlist
    ##+ object : massvector
    ##+ accur : measurment accuracy in m/z
    ##+ abund : how many times a mass have to occur in a mass bin to be considered as an abundant mass.
    ##+ ... : further parameters.
    ##v massvector : massvector with abundant masses.
    ##r Wolski 
    ##e data(mvl)
    ##e mv<-unlist(mvl) 
    ##e res<-gamasses(mv,abund=30,accur=0.3)
    ##e plot(res)
    
    dat<-object[,1]
    attributes(dat)<-NULL
    mhist<-list(NULL)
    mhist[[1]] <-hist(dat,breaks=seq(min(dat)-accur/2,max(dat)+1.5*accur,accur),plot=FALSE)#,main=main,xlab=xlab,xlim=xlim)
    mhist[[2]]<-hist(dat,breaks=seq(min(dat)-accur,max(dat)+accur,accur),plot=FALSE)
    wh1<-which(mhist[[1]]$counts > abund)
    wh2<-which(mhist[[2]]$counts > abund)
    if(length(wh1)==0 & length(wh2)==0)
      {
        res<-matrix(ncol=3)
        colnames(res)<-c("mass","number","wm")
        return(massvector(paste(info(object),"abundant",sep="_"),res))
      }
                                        #    cat("wh1 : ", wh1 ," \nwh2 : ",wh2,"\n")
                                        #bestimme die maxima.
    z1<- rep(0,length(mhist[[1]]$mids))
    z2<- rep(0,length(mhist[[2]]$mids))
                                        #    cat("z1 " , z1 , " z2 " , z2 , "\n" )
                                        #übernehme nur die die häufig sind.
    z1[wh1]<- mhist[[1]]$counts[wh1]
    z2[wh2]<- mhist[[2]]$counts[wh2]
                                        #find the peaks.returns the indexes    
    p1 <- peaks(z1,max=TRUE)
    p2 <- peaks(z2,max=TRUE)
    names(p1)<-rep(1,length(p1))        #name the indexes
    names(p2)<-rep(2,length(p2))        #name the indexes
                                        #Determin which indexes are close to each other.
    p12 <- sort(c(p1,p2))
                                        #if difference between peaks are 0 or 1 than they are close to each other.
    dp <- diff(p12)
                                        #find the indexes of this cases
    ipl <- which(dp==1|dp==0)
    ipu <- ipl+1
                                        #now i have the indices of the peaks in array p12.
                                        #Array p12 by themselves gives me the indices of the peaks.
    p12close <- cbind(p12[ipl],p12[ipu])
    p12closeN <- cbind(as.numeric(names(p12)[ipl]),as.numeric(names(p12)[ipu]))
                                        #now in each row are the indices of the peaks that are close to each other.
                                        #using this indices i can retrieve the mids and the counts
                                        #for this bins out of the histogram and calculate a weighted
                                        #average.
    wm <- NULL

    if(length(p12close)>0)
      {
        for(x in 1:length(p12close[,1]))
          {
            w1 <- mhist[[ p12closeN[x,1] ]]$counts[p12close[x,1]]
            m1 <- mhist[[ p12closeN[x,1] ]]$mids[p12close[x,1]]
            w2 <- mhist[[ p12closeN[x,2] ]]$counts[p12close[x,2]]
            m2 <- mhist[[ p12closeN[x,2] ]]$mids[p12close[x,2]]
            wm <- c( wm , (w1*m1+w2*m2)/(w1+w2) )
          }
      }

    if(length(c(ipl,ipu)>0))
      {
        uniquepeaks <- p12[-c(ipl,ipu)]
      }
    else
      {
        uniquepeaks <- p12
      }
    uN <- as.numeric(names(uniquepeaks))

    if(length(uniquepeaks)>0)
      {
        for(x in 1:length(uniquepeaks))
          {
            wm<-c(wm,mhist[[ uN[x] ]]$mids[uniquepeaks[x]]) 
          }
      }
    rm(mhist)
    res <- NULL
    num<-NULL
                                        #calculating exact mass
    for(x in wm)
      {
        #print(x)
        res <- c(res, mean( dat[(x+accur*0.8) > dat & (x-accur*0.8) < dat ]))
        num<-c(num,length(dat[(x+accur*0.8) > dat & (x-accur*0.8) < dat ]))
      }
                                        #    print(paste(info(object),"abundant",sep="_"))

    res<-cbind(res,num,wm)
    colnames(res)<-c("mass","number","wm")
    return(massvector(paste(info(object),"abundant",sep="_"),res))
  }


diffFilter.massvector<-function(object,listofdiffs,higher=TRUE,error=0.05,uniq=TRUE,...)
  {
    ##t Abundant Differences
    ##- Removes mass differences from the massvector.
    ##d Removes one of the masses contributing to a mass difference given in the list of differences.
    ##d Can be used if a variable modification are present in the massvector but can not be considered by the identification software.
    ##+ object : massvector
    ##+ listoffdiffs : massvector with mass differences
    ##+ higher : logical;\code{TRUE} - remove higher mass, \code{FALSE} = remove lower mass.
    ##+ error : How much the differences can diviate from the differences given in listofdiffs.
    ##+ ... : further parameters.
    ##v massvector : filtered massvector.
    ##sa getdiff.massvector, getdiff.massvectorlist, diffFilter.massvectorlist
    ##r Wolski
    ##e data(mv1)
    ##e res<-getdiff(mv1,range=c(0,100))
    ##e diffFilter(mv1,res,higher=TRUE)
    ##e diffFilter(mv1,res,higher=FALSE)
    
    
    if(length(object)<=1)
      {
        return(object)
      }
    pl <- object[,1]
    pl <- sort(pl)
    res <- NULL
    ldiff<-listofdiffs[,1] # get the difference masses.
    mind <- min(ldiff)-1
    maxd <- max(ldiff)+1
    
    for(x in 1:(length(pl)-1))
      {
        diffpl <- sort(diff(pl,x))
        
        tmp <- getaccC(diffpl,ldiff,error=error,ppm=FALSE,uniq=uniq)
        ind <- as.numeric(names(diffpl))[tmp$plind]
        if(higher)
          {
            res <- c(res,(ind))
          }
        else
          {
            res<-c(res , ind-x )
          }
      }
    if(length(res)>0)
      return(object[-res,])
    else
      return(object)
  }

getdiff.massvector<-function(object,range=c(0,100),...)
  {
    ##t Massdifferences
    ##- Computes mass differences in the object.
    ##d Removes one of the masses contributing to a mass difference given in the list of diffs.
    ##d Can be used if a variable modification are present in the massvector but can not be considered by the identification software.
    ##+ object : massvector
    ##+ listoffdiffs : massvector with mass differences
    ##+ higher : \code{TRUE} - remove higher mass, \code{FALSE} = remove lower mass.
    ##+ error : How much the differences can diviate from the differences given in listofdiffs
    ##v massvector : filtered massvector.
    ##r Wolski

    rrr<-range
    res <- NULL
    rarea <- NULL
    if(!length(object)>1)
      {
        res<-massvector(info(object),NULL)
        setParms(res)<-list(tcoor=mget(object,"tcoor"))
        return(res)
      }
    else
      {
        m <- mass(object)
        a<-object[,2]
        for(x in 1:(length(m)-1))
          {
            diffpl<-m - m[x]
            ratarea<-a / a[x]
            names(diffpl)<-NULL
            names(ratarea)<-NULL
            res <- c(res,diffpl[diffpl>rrr[1] & diffpl<rrr[2]])
            rarea <- c(rarea,ratarea[diffpl>rrr[1] & diffpl<rrr[2]])
          }

      }
    if(length(res)==0)
      {
        res<-massvector(info(object),NULL)
        setParms(res)<-list(tcoor=mget(object,"tcoor"))
        return(res)
      }
    res<-cbind(res,rarea)
    colnames(res)<-c("massd","arear")
    res<-massvector(info(object),res)
    setParms(res)<-list(tcoor=mget(object,"tcoor"))
    return(res)
  }

writeF.massvector<-function(object,path,file=info(object),ext="txt",...)
  {
    ##t Write massvector
    ##- Write massvector to File
    ##+ object : massvector
    ##+ path : path to folder.
    ##+ file : file name. defualt info(object)
    ##+ ext : file extension.
    ##sa readBruker.massvector,readBruker.massvectorlist,writeF.massvectorlist
    ##e data(mv1)
    ##e writeF(mv1,".") # writes the file in the home directory.
    ##e readF(massvector(info(mv1)),".")
    ##e file.remove(paste(info(mv1),".txt",sep=""))
    
    if(missing(path))
      {
        path<-"."
      }
    
    filep <- file.path(path,paste(file,".",ext,sep=""),fsep = .Platform$file.sep)
    con <- file(filep, "w")  # open an output file connection
    firstline<-paste(">",file,":",join(mget(object,"tcoor"),sep=","),sep="")
    writeLines(firstline,con=con)
    close(con)
    write.table(object,file = filep, append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))
  }


readF.massvector<-function(object,path,file=info(object),ext="txt",...)
  {
    ##t Read Massvector
    ##- Reads massvector written with the function \code{writeF.massvector}
    ##+ object : object of class massvector. Use constructor \code{massvector()}
    ##+ path : path to file.
    ##+ file : file name; default : info(object).
    ##+ ext : file extension; default : txt.
    ##sa
    ##e data(mv1)
    ##e writeF(mv1,".")
    ##e readF(massvector(info(mv1)),".")
    ##e file.remove(paste(info(mv1),".txt",sep=""))
    filep <- file.path(path,paste(file,".",ext,sep="") , fsep = .Platform$file.sep)
    con<-file(filep,"r")
    res <- readLines(con=filep,n=-1)
    close(con)
    #return(res)
    
    res1 <- unlist(strsplit(res[1],":"))
    object<-info(object,sub(">","",res1[1]))
    
    setParms(object) <- list(tcoor= unlist(strsplit(res1[2],",")))
   # return(object)
    intern<-function(x){as.numeric(unlist(strsplit(x,"\t")))}
    res <- t(sapply(res[2:length(res)],intern))
    res<-res[order(res[,1]),]
    rownames(res)<-1:length(res[,1])
    object <- peaks(object,res)
    return(object)
  }

