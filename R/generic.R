#Copyright 2004, W. Wolski, all rights reserved.
                                        #generic methods provided by the package
                                        #function used at least by more than one object.
                                        #plot
                                        #summary
                                        #hist
                                        #print

#comparing two massvectors.
compare<-function(object,...)
  {
    UseMethod("compare")
  }

as.lm <- function(object,...)
  {
    UseMethod("as.lm")
  }

#for setting and getting masses
mass <- function(object,...)
  {
    UseMethod("mass")
  }

#each object has a unique human readable field
info<-function(object,...)
  {
    UseMethod("info")
  }
#is the unique identifier either readable or not readable
#id<-function(object,...)
#  {
#    UseMethod("id")
#  }
#is either the id, or the info
#key<-function(object,...)
#  {
#    UseMethod("key")
#  }


##as.data.frame

##(replaces correct affine)
#internalCalib<-function(object,...)
#  {
#    UseMethod("internalCalib")
#  }


##searches the mascot server with the list.#
#maskotSearch<-function(object,...)
#  {
#    UseMethod("maskotSearch")
#  }

###################################################
##
## the Wool and Smilanski filtering method
##

wsFilter<-function(object,...)
  {
    UseMethod("wsFilter")
  }


wsiFilter<-function(object,...)
  UseMethod("wsiFilter")

wsdist<-function(object,...)
    UseMethod("wsdist")


#p1 -have to be a double
#p2 -can be a vector
distance <- function(p1,p2)
  {
    ##t Intra massvector mass distance
    ##- Distance between two peaks o a massvector defined as deviation from the peptide rule.
    ##+ p1 : mass
    ##+ p2 : mass
    ##v distance : distance of the two masses, as deviation from the peptide rule.
    ##r Wool A, Smilansky Z 2002. Precalibration of matrix-assisted laser desorption/ionization-time of flight spectra for peptide mass fingerprinting. {\em Proteomics.} 2(10):1365-73.
    ##r Wolski http://www.molgen.mpg.de/~wolski/mscalib
    d12<-abs(p1-p2)
    m<-1.000495
    return(mmod(d12,m))
  }

mmod <- function(x,m)
  {
    mn <- mod(x,m)
    mnt<-mn
    mns<-mnt[mnt > m/2]
    mn[mnt > m/2] <- (m - mns)
    return(mn*2)
  }

mod <-function(x,m)
  {
    t1<-floor(x/m)
    return(x-t1*m)
  }

#
#calibEternal functions
#

getextcalib <- function(object,...)
  {
    UseMethod("getextcalib")
  }


##calibrates a peaklist by a spline
applyextcalib<-function(object,...)
  {
    UseMethod("applyextcalib")
  }

calibexternal<-function(object,...)
  {
    UseMethod("calibexternal")
  }
###################################################
#
#Recalibration functions
#

recalibrate<-function(object,...)
  {
    UseMethod("recalibrate")
  }
##returns recalibobject.
getrecalib<-function(object,...)
  {
    UseMethod("getrecalib")
  }

applyrecalib<-function(object,...)
  {
    UseMethod("applyrecalib")
  }

applycalib<-function(object,...)
  {
    UseMethod("applycalib")
  }


FullWidthatHalfMaximum<-function(sumcol)
{
   minmax<-range(sumcol[,2])
                                        #   compute the half maximum
   abl1<-diff(sumcol[,2])
   ta1a<-c(abl1,1)
   ta1b<-c(1,abl1)
   minima <- which(ta1a>0 & ta1b<0)
   mmax <- which(sumcol[,2]==max(sumcol[,2]))
   if(mmax==1 | mmax==length(sumcol[,2]))
     {
          return(list(PQM = 0,hm=c(0,0,0),width=c(0,0)))
     }
   else
     {
       ll <- which(minima<mmax)
       rr <- which(minima>mmax)
       if(length(ll)>0)
         {
           lmin<- max(minima[ll])
         }
       else
         {
           lmin<-1
         }
       
       if(length(rr)>0)
         {
           rmin<- min(minima[rr])
         }
       else
         {
           rmin<-length(ta1a)
         }
       
       int<-sumcol[,2]
                                        #peak hight.
       pmax <- int[mmax] # peak maximum
       pmin <- max(c(int[lmin],int[rmin])) #peak minimum
       ph <- pmax - pmin
       if(FALSE)
         {
           plot(sumcol[,2],type="l")
           
           abline(h=pmax)
           abline(h=pmin)
           abline(v=mmax,col=2)
           abline(v=rmin,col=3)
                                        #   maxima <- which(tt<0 & tt2>0)
           abline(v=lmin,col=4)
         }
                                        #peak width
                                        #halfmax
       hm<-ph/2
       sumcc<-sumcol[lmin:rmin,]
       
       
                                        #hm <- (minmax[1] + minmax[2])/2
       dumm<-0.05
       while(TRUE){
         hi<-(pmin + hm + hm*dumm)
         lo<-(pmin + hm - hm*dumm)
         
         width <-sumcc[,1][sumcc[,2]<hi &sumcc[,2]>lo]
         if(length(width)>2) break;
         dumm<-dumm+0.05
       }
       if(FALSE)
         {
           plot(sumcc,type="l")
           abline(h=hi)
           abline(h=lo)
         }

       ind <- which(diff(width)==max(diff(width)))
       width <- c(mean(width[1:ind]),mean(width[(ind+1):length(width)]))
       
       hm<-c(pmin + hm,pmax,pmin)
       names(hm)<-c("HalfMaximum","max","min")
       width<-c(width[2]-width[1],width)
     }
   return(list(PQM=-log(as.numeric(width[1]/ph)),hm=hm,width=width))
}


####################################################
##Affine calibration.
#


##expects an lm object as second param
correctinternal <- function(object,...)
  UseMethod("correctinternal")

getintcalib <- function(object,...)
  UseMethod("getintcalib")

applyintcalib<-function(object,...)
  UseMethod("applyintcalib")

####################################################
##global calibration
#
getglobalcalib<-function(object,...)
  UseMethod("getglobalcalib")

globalcalib<-function(object,...)
    UseMethod("globalcalib")

globalcalib.default<-function(object,...)
    warning("need massvectorlist\n")


#
#get abundant masses
#
gamasses<-function(object,...)
  UseMethod("gamasses")

#
#getting and setting calibobjects
#
getcalib<-function(object,...)
  UseMethod("getcalib")

setcalib<-function(object,...)
  UseMethod("setcalib")

getdiff<-function(object,...)
  UseMethod("getdiff")

diffFilter<-function(object,...)
  UseMethod("diffFilter")


##reads a peaklist from a file.
readBruker<-function(object,...)
    UseMethod("readBruker")

peaks<-function(object,...)
  UseMethod("peaks")


peaks.default<-function(object, max=TRUE,na.rm=FALSE,...){
  ##t Find peaks
  ##- Finds peaks - neighborpeaks smaller than central.
  ##+ object : numeric array.
  ##+ max : TRUE find maxima, FALSE find minima
  ##+ na.rm : handling of na.rm values.
  ##v index : index of the peaks   
  x<-object
  if (na.rm)
   omit<-is.na(x)
  else
   omit<-FALSE
  if (max){
   rval<-1+which(diff(sign(diff(x[!omit])))<0)
  }else{
   rval<-1+which(diff(sign(diff(x[!omit])))>0)
  }
  if (na.rm)
  {
   rval<-rval+cumsum(omit)[rval]
  }
  rval
}


#method for filtering massvectors and massvectorlists
mvFilter<-function(object,...)
  UseMethod("mvFilter")

join<-function(object,...)
    UseMethod("join")

join.default<-function(object,sep,...)
{
  if(missing(sep))
    sep<-""
  res <- object[1]
  if(length(object)>1)
    {
      for(y in 2:length(object))
        {
          res<-paste(res,object[y],sep=sep)
        }
    }
  return(res)
}

"setParms<-"<-function(object,value)
  {
    UseMethod("setParms<-")
  }


experiment<-function(object, ... )
  UseMethod("experiment")


project<-function(object, ... )
  UseMethod("project")

                                        #gets an object field
mget<-function(object,...)
    UseMethod("mget")

writeF<-function(object,...)
  UseMethod("writeF")


readF<-function(object,...)
  UseMethod("readF")
