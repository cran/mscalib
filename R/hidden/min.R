
calib  <-  function()
  {
    ##t Massvector with Calibration Masses
    ##- Returns massvector with masses for internal calibration. (Mainly trypsin masses.)

    calibn<-c("calib1"
              ,"calib2"
              ,"calib3"
              ,"tr1"
              ,"tr2"
              ,"tr3"
              ,"tr4"
              ,"tr5"
              ,"tr6"
              ,"tr7"
              ,"tr8"
              ,"acth"
              ,"newc0"
              ,"newc1"
              ,"newc2"
              ,"newc3"
              ,"ab1"
              ,"ab2"
              ,"ab3"
              ,"ab4"
              ,"ab5"
              ,"ab6"
              ,"ab7"
              ,"ab8"
              ,"ab9"
              )

    calib  <- 
      c(1296.6853
        ,1672.9175
        ,3147.4715
        ,842.5099
        ,1045.5642
        ,2211.1046
        ,2225.1203
        ,2239.1359
        ,2283.1807
        ,2297.1964
        ,2313.1913
        ,2465.198
        ,1126.5478
        ,2807.3290
        ,2839.3048
        ,3346.712
        ,720.40 
        ,863.41 
        ,1131.55
        ,1232.54
        ,1352.70
        ,1638.78
        ,2069.93
        ,2085.04
        ,3846.12
        )

    names(calib)  <-  calibn 
    calib <- massvector("calib",calib)
    return(calib)
  }


#
#The same as above but makes a list
#

readBrukerFiles <- function(folderName,filename,func = peaklistxml.short)
{
    require(XML)
    ntemp <- character(0)
    otemp <- character(0)
    temp <- dir(folderName,pattern="*Ref$")
    tsb <- proc.time()
    times <- NULL
    for(x in 1:length(temp))
    {
        ntemp <- c(ntemp,paste(folderName,"/",temp[x],"/pdata/1/",filename,sep=""))
        otemp <- c(otemp,paste(folderName,"/",temp[x],sep=""))
    }
    plot(0,0,xlim=c(0,length(ntemp)),ylim=c(0,1),pch="|",col=2,axes=FALSE)
        retlist <- list(NULL)
        names(retlist)[1] <- temp[1]
        retlist[[1]] <- func(ntemp[1])
        points(1,0,pch="|",col=2)
        for(x in 2:length(ntemp))
        {
            prt <- proc.time()
            tmp <- func(ntemp[x])
            now <- proc.time()
            times <- c(times,now[1]-tsb[1])
            if(length(tmp)>0)
            {
                retlist[[length(retlist)+1]] <- tmp
                names(retlist)[length(retlist)] <- temp[x]
            }
            points(x,0,pch="|",col=2)
        }
        cat("length retlist",length(retlist),"\n")
        return(retlist)  
}


readMPIIBPeaklist.folder<-function(folder,pattern="*.PKM")
  {
    #reads a whole file folder.
    ntemp <- dir(folder,pattern=pattern)
    plot(0,0,xlim=c(0,length(ntemp)),ylim=c(0,1),pch="|",col=2,axes=FALSE)

    retlist <- list(NULL)
    names(retlist)[1] <- ntemp[1]
        cat("test",ntemp[1],"\n")
    retlist[[1]] <- readMPIIBPeaklist(ntemp[1])
    print("test")
    points(1,0,pch="|",col=2)

        for(x in 2:length(ntemp))
        {
            prt <- proc.time()
            tmp <- readMPIIBPeaklist(ntemp[x])
            now <- proc.time()
            if(length(tmp)>0)
            {
                retlist[[length(retlist)+1]] <- tmp
                names(retlist)[length(retlist)] <- ntemp[x]
            }
            points(x,0,pch="|",col=2)
        }
        cat("length retlist",length(retlist),"\n")
        return(retlist)  
  }

readMPIIBPeaklist<-function(file)
  {
    tt<-readLines(file)
    ta<-grep("^[0-9]+",tt)
    tt<-tt[ta]
    tt<-strsplit(tt," +")
    retstrut<-NULL
    for(x in 1:length(tt))
      {
                                        #   print(unlist(tt)[c(1,2)])
        retstrut<-rbind(retstrut,as.double(unlist(tt[[x]])[c(1,2)]))
      }
    retstrut<-data.frame(retstrut)
    names(retstrut)<-c("mass","int")
    return(retstrut)
  }

###
# reads from the peaklist xml file the entries. mass absi & ind.
#
##
peaklistxml.short<-function(fileName)
{
    require(XML)
    #print("the trap begins")
    mass <- NULL
    absi <- NULL
    ind <- NULL
    
    test<-function() { 
        #vars <- character(0) ;
        massf<-function(x,attrs){
            mass <<- c(mass, xmlValue(x[[1]])); 
            NULL
        }
        absif <- function(x,attr) {
            absi<<-c(absi,xmlValue(x[[1]]));
            NULL
        }
        indf <- function(x,attr){
            ind<<-c(ind,xmlValue(x[[1]]))
        }
        list(mass = massf,absi = absif, ind=indf)
    }
    
    ##################################################
    if(file.exists(fileName)){
    # retrieve the names of the file
             tt  <-  strsplit(fileName,"/")
             pkname <- tt[[1]] [(1:length(tt[[1]]))[tt[[1]]=="pdata"]-1 ]
        #cat("before the trap",fileName,"\n")
             xmlTreeParse(fileName,handler=test())
        if(length( mass ) == length( absi ) & length( absi ) == length(ind ) & length( absi )!=0)
        {
            tmp <- data.frame(as.character(rep( pkname, length(mass) ) ) , 1:length(mass), as.double(mass) , as.double(ind) , as.double(absi) )
        }else{
           cat("pkname = ", pkname , " mass = ", length(mass) , " absi = ", length(absi) , " ind = ", length(ind)  , "\n" )
           return(NULL)
        }
        names(tmp) <- c("pkname","index","mass","ind","absi")
        return(tmp)
       
     }
     else
     {
        return(NULL)
    }
}
#peaklist access and writing.
#writes a column in a dataframe from a list of dataframes
#to files called peptidepeaklist in subfolders named like the list entries. 

writeListBruker <- function(lod,folder,masst="massr")
{
    name <- names(lod)
    for(x in 1:length(lod))
    {
       mkdir(paste(folder,name[x],"Rpeaklist",sep="/"))
       pathfile  <-  paste( folder,"/",name[x],"/","Rpeaklist",sep="" )
       #print(pathfile)
       dataf <- lod[[x]]
       mind  <-  (1:length(dataf))[names(dataf)==masst]
       write.table(dataf[[mind]],file = pathfile, row.names = FALSE, col.names = FALSE ,quote=FALSE)
    }
}


#---------------------------------------------------


#---------------------------------------------------------------------
# Needed for internal calib to list.

getacc2<-function(peaklist,calib,error=500)
  {
    ret<-NULL
    plind<-NULL
    calind<-NULL
    for(x in 1:length(calib))
      {
                                        #calculates the deviation in ppms
        tmp <- abs((peaklist-calib[x]))/calib[x]*1e6
        ind<-which(tmp<error)
        plind<-c(plind,ind)
        calind<-c(calind,rep(x,length(ind)))
      }
    return(list(plind=plind,calind=calind))
  }

getacc<-function(peaklist,calib,error=500)
  {
    ret<-NULL
    plind<-NULL
    calind<-NULL
    for(x in 1:length(calib))
      {
                                        #calculates the deviation in ppms
        tmp <- abs((peaklist-calib[x]))/calib[x]*1e6
        ind<-which(tmp<error)
        if(length(ind)>0)
          {
            mins<-tmp[ind]
            nind<-ind[mins==min(mins)]
            plind<-c(plind,nind)
            calind<-c(calind,x)
          }
      }
    return(list(plind=plind,calind=calind))
  }

#--------------------------------------------------------------------------
# needed for getListforCalib.
#external calibrations


getPPGmasses <- function(start=10,end=100)
{
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
        return(n)
}


#----------------------------------------------
##
#calculates the spline function from the return value of getListforCalib
#it returns an spline object this spline object can be used for calibration of ppg spectra by the function
#returns a list with the spline prediction used for calibration the theoretical time and the averaged time
##

getCalibSpline.lod <- function(lodppg,calib,masst="massa",error=250,more=FALSE)
{
  
  mind  <-  (1:length(names(lodppg[[1]])))[ names(lodppg[[1]]) == masst]
  #ppg<-getPPGmasses()
  ppg<-calib
  theo<-NULL
  expd<-NULL
  for(x in 1:length(lodppg))
    {
      exp<-lodppg[[x]][[mind]]      
      match<-getacc(exp,ppg,error=error)
      expd<-c(expd,exp[match$plind])
      theo<-c(theo,ppg[match$calind])
    }
  mord<-order(theo)
  expd<-expd[mord]
  theo<-theo[mord]
  cat("exp",print(length(expd)),"theo",print(length(theo)),"\n")
  require(modreg)
  #calculate mass dependent error
  error<-(expd-theo)*1e6/theo
  #ispl<-lm(error~ns(theo,df=20))
  ispl<-smooth.spline(theo,error)

  if(more)
    {
      plot(theo,error,pch="*",xlab="thoretical mass [m/z]",ylab="error [ppm]")
      lines(predict(ispl,theo),col=2)
      return(list(theo=theo,error=error,ispl=ispl))
    }else{
      return(ispl)
    }
}

getCalibSpline<-function(peaklist,calib,error=250,more=FALSE)
  {
    match<-getacc(peaklist,calib,error=error)
    expd<-peaklist[match$plind]
    theo<-calib[match$calind]
    mord<-order(theo)
    expd<-expd[mord]
    theo<-theo[mord]
    cat("exp",print(length(expd)),"theo",print(length(theo)),"\n")
    require(modreg)
                                        #calculate mass dependent error
    error<-(expd-theo)*1e6/theo
                                        #ispl<-lm(error~ns(theo,df=20))
#    nknots<-round(0.8*length(expd))
    ispl<-smooth.spline(theo,error,all.knots=TRUE)
#    require(splines)
    ispl<-interpSpline(theo,error,period=4,bSpline=TRUE)

    if(more)
    {
      plot(theo,error,pch="*",xlab="thoretical mass [m/z]",ylab="error [ppm]")
      lines(predict(ispl,theo),col=2)
      return(list(theo=theo,error=error,ispl=ispl))
    }else{
      return(ispl)
    }
  }



#----------------------------------------------------------
# methods for apply spline to spectra.
######################################
# Functions for external calibration
######################################
#calibration
#package unifying all functions needed for peaklist calibration.
# parameter - a bunch of spectra to which the linear calibration have to be applied.
# the calibrateted values are appended as a new row.
# spList is the List representation of the peaklists. 

applySpline.lod <- function(lod , spline,massr="massa",massw="masse")
{
  #lod - list of peaklists
  #spline to use for predicting the error function.
  #massr- mass to read.
  #massw-mass to write.
  mind  <-  (1:length(names(lod[[1]])))[ names(lod[[1]]) == massr]
  require(splines)
  for(x in 1:length(lod))
  {
    lod[[x]][[massw]] <- applySpline(lod[[x]][[mind]],spline)
  }
  return(lod)
}

applySpline<-function(peaklist,spline)
{
  #peaklist- array of peaks masses.
  #spline - spline to predict the error.
  error <- predict(spline,peaklist)
  error <- error$y
  #error<-predict(spline,peaklist)
  masspred <- peaklist/(1+error/1e6)
  return(masspred)
}




#--------------------------------
#peaklist rekalibration. is working.
#Recalibration.Wool<-function(){}

recalibrate.lod<-function(lod, massr="mass", massw="massr")
  {
    #list of peaklists to reakalibrate.
    #massr - masst to read.
    #massw - masst to write
    dataf<-lod[[1]]
    mind  <-  (1:length(dataf))[names(dataf)==massr]
    
    for(x in 1:length(lod))
      {
        dataf<-lod[[x]]
        peaklist<-dataf[[mind]]
        peaklist<-recalibrate(peaklist)
        lod[[x]][[massw]]<-peaklist
      }
    return(lod)
    
  }


recalibrate<-function(peaklist,plot=FALSE)
  {
    if(length(peaklist)<2)
      {
        return(peaklist)
      }
    lambda<-seq(0.999495,1.001495,0.000001)
    omega<-2*pi/lambda
    test<-peaklist%*%t(omega)
    testsin<-sin(test)
    testcos<-cos(test)
    colsumsin <-apply(testsin,2,sum)
    colsumcos <-apply(testcos,2,sum)
    sumcol<-sqrt(colsumsin^2+colsumcos^2)
    sumcol <- cbind(lambda,sumcol)
    #determine the maximum and get its index
    mmax <- max(sumcol[,2])
    index<-sumcol[,2]==mmax
    #temporary only for testing purposes.
    if(plot){
      len<-length(peaklist)
      main <- paste("Peaklist length : ", len, sep="")
      plot(sumcol,type="l",ylab="Amplitude",xlab="Wavelength",xlim=c(0.9995,1.0015),las=2, main = main)
      points(sumcol[,1][index],sumcol[,2][index],col=2,pch="*")
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
                                        #cat("bmax : ", bmax,"\n");
                                        #cat("aplpha : ",alpha , "\n"); 
    peak<- (peaklist-bmax)*alpha
    return(peak)
  }



#------------------------------------------------
# constants
#------------------------------------------------

calibn  <-  c("calib1"
,"calib2"
,"calib3"
,"tr1"
,"tr2"
,"tr3"
,"tr4"
,"tr5"
,"tr6"
,"tr7"
,"tr8"
,"acth"
,"newc0"
,"newc1"
,"newc2"
,"newc3"
,"ab1"
,"ab2"
,"ab3"
,"ab4"
,"ab5"
,"ab6"
,"ab7"
,"ab8"
,"ab9"
)

calib  <- 
c(1296.6853
,1672.9175
,3147.4715
,842.5099
,1045.5642
,2211.1046
,2225.1203
,2239.1359
,2283.1807
,2297.1964
,2313.1913
,2465.198
,1126.5478
,2807.3290
,2839.3048
,3346.712
,720.40 
,863.41 
,1131.55
,1232.54
,1352.70
,1638.78
,2069.93
,2085.04
,3846.12
)

names(calib)  <-  calibn 
calib <- sort(calib)

#----------------------------------------
#Display functions.

imageplot<-function(lod,zoom = NULL,masst = "masslin",main="")
{
    allm  <- getcol(lod,masst=masst)
    print("test")
    mmin  <- min(allm)
    mmax <- max(allm)
    dataf<-lod[[1]]
    mind  <-  (1:length(dataf))[ names(dataf) == masst ]
    print("test")
    if(is.null(zoom))
    {
    
      plot(lod[[1]][[mind]], rep(1,length(lod[[1]][[mind]])) ,xlim = c(mmin,mmax),ylim = c(1,length(lod)),pch=15,cex=0.5,xlab="m/z",ylab="sample",main=main)
    }
    else{
        plot(lod[[1]][[mind]], rep(1,length(lod[[1]][[mind]])) ,xlim = zoom,ylim = c(1,length(lod)),pch=15,cex=0.5,xlab="m/z",ylab="sample",main=main)
    }
    for( x in 2:length(lod))
    {
        points(lod[[x]][[mind]], rep(x,length(lod[[x]][[mind]])),pch=15,cex=0.5)
    }
    abline(h=1:length(lod),lty=3)
}

imageplotadd<-function(lod,masst="massr",col=2)
  {
    dataf<-lod[[1]]
    mind  <-  (1:length(dataf))[ names(dataf) == masst ]
    for( x in 1:length(lod))
    {
        points(lod[[x]][[mind]], rep(x,length(lod[[x]][[mind]])),pch=15,cex=0.5,col=col)
    }
  }

##
# becaues the structure used in most functions is the list of spectra
# a funcition merging the list is needed. 
# its called getcol(colname)
##
getcol <- function(lod, masst="")
{
 retarr <- NULL
    if(is.data.frame(lod[[1]]))
    {
    mind  <-  (1:length(names(lod[[1]])))[ names(lod[[1]]) == masst]
   

  for(x in 1:length(lod))
  {
    retarr <- c(retarr,lod[[x]][[mind]])
  }
  return(retarr)
  }
 
  if(is.double(lod[[1]]))
  {
    for(x in 1:length(lod))
    {
        retarr <- c(retarr,lod[[x]])
    }
    return(retarr)
  }
 return(retarr)
}
#returns a vector with peaklist lengths.
reportListHist <- function(lod)
{
    lengths <- rep(0,length(lod))
    for(x in 1:length(lod))
    {
        if(class(lod[[x]])=="data.frame")
        {
            lengths[x] <- length(lod[[x]][[1]])
        }else
        {
            lengths[x] <- length(lod[[x]])
        }
    }
    hist(lengths)
    return(lengths)
}



#--------------------------------------------------------
new.functions.for.internal.calibration<-function(){}
###
#
#  
correctaffine<-function(peaklist,calib,error=500,func=getacc2)
  {
                                        #peaklist- peaklist.
                                        #calib-list with calibrants.
    match<-func(peaklist,calib,error=error)
    if(length(match$plind)>0)
      {
                                        #calculate mass dependent error function.
                                        #get peaklist matching
        smallerrpl<-peaklist[match$plind]
                                        #get calibrants matching
        masstheo<-calib[match$calind]
                                        #finally calculate error
        err  <-  ( smallerrpl -masstheo)*1e6 /masstheo
                                        #distinguishing 2 cases. if peaks quite close to each other only correct for offset.
        if(diff(range(masstheo))<300)
          {
            mymod<-lm(mean(err)~mean(masstheo))
          }else{
            mymod<-lm(err~masstheo)          }
        errp <- predict(mymod,data.frame(masstheo=peaklist ) )
                                        #print(masspred)
        peaklist <- peaklist/(1+errp/1e6)
                                        #return(predict(mymod,data.frame(peaklist)))
      }
    return(peaklist)
  }


#tpeaks<-spectra2[[1]]$mass
#tmp<-correctaffine(tpeaks,calib,error=500)
## Predictions

correctaffine.lod<-function(lod,calib,error=500,massr="mass",massw="massa",func=getacc2)
  {
    # list of peaklists.
    # calib - array with  masses to use for calibration
    # massr - mass to read.
    # massw - mass to write.
    
    dataf<-lod[[1]]
    mind  <-  (1:length(dataf))[names(dataf)==massr]
    for(x in 1:length(lod))
      {
        dataf<-lod[[x]]
        peaklist<-dataf[[mind]]
        peaklist<-correctaffine(peaklist,calib,error=error,func=func)
        lod[[x]][[massw]]<-peaklist
      }
    return(lod)
  }

