test<-function(x)
 {
 print("hello");
 	res<-.C("mytest",
                as.integer(x),
                )
	return(res)
 }


getaccC<- function(pl,cal,error=500,ppm=TRUE,uniq=FALSE)
{
  #print("test")
  lpl <- length(pl)
  lcal <- length(cal)
  
  #cat("lpl ",lpl, " lcal ",lcal," lmods ",lmods,"\n")
  #print(lmods)
  tmp <- max(lpl,lcal) # there can not be more matches than tmp.
  plind <- rep(0,tmp)
  calind <- rep(0,tmp)
  ind <- 0
  if(length(pl)>0 & length(cal)>0)
  {
    res <- .C("getaccD", 
              as.double(pl), 
              as.integer(lpl),  
              as.double(cal), 
              as.integer(lcal),
              as.double(error),
              as.integer(plind),
              as.integer(calind),
              as.integer(ind),
              as.integer(ppm),
              as.integer(uniq)
       )
   }
   else{
    return(list(plind=NULL,calind=NULL))
   }     

  if(res[[8]]>0)
  {
     if(uniq)
     {
        plind <- res[[6]][1:res[[8]]]+1
        calind <- res[[7]][1:res[[8]]]+1
        dpl <- which(diff(plind)==0)
        dcl <- which(diff(calind)==0)
        if(length(dpl)>0 | length(dcl)>0)
          {
    rdp <- NULL
    rcl <- NULL
            if(length(dpl) > 0)
              {
            
                for(dp in dpl)
                  {
                    if(abs(pl[plind[dp]]-cal[calind[dp]]) < abs(pl[plind[dp+1]] - cal[calind[dp+1]]))
                      {
                        rdp<-c(rdp,(dp+1))
                      }else
                    {
                      rdp<-c(rdp,dp)
                    }
                  }
              }
            if(length(dcl) > 0)
              {
                
                for(dp in dcl)
                  {
                    if(abs(pl[plind[dp]]-cal[calind[dp]]) < abs(pl[plind[dp+1]] - cal[calind[dp+1]]))
                      {
                        rcl<-c(rcl,(dp+1))
                      }else
                    {
                      rcl<-c(rcl,(dp+1))
                    }
                  }
              }
            pcl<-unique(c(rdp,rcl))
                                        #print(pcl)
            plind<-plind[-pcl]
            calind<-calind[-pcl]
          }
        
        test<-list(plind=plind,calind=calind)   
     }else{
      test <- list(plind=res[[6]][1:res[[8]]]+1,calind=res[[7]][1:res[[8]]]+1)
     }
  }
  else{
    test <- list(plind=NULL,calind=NULL)
  }
   return(test)
}
  

getaccT<-function(pl,cal,error=500,ppm=TRUE,uniq=TRUE)
{
  #print("test")
  lpl <- length(pl)
  lcal <- length(cal)
  
  #cat("lpl ",lpl, " lcal ",lcal," lmods ",lmods,"\n")
  #print(lmods)
  tmp <- max(lpl,lcal) # there can not be more matches than tmp.
  plind <- rep(0,tmp)
  calind <- rep(0,tmp)
  ind <- 0
    res <- .C("getaccD", 
              as.double(pl), 
              as.integer(lpl),  
              as.double(cal), 
              as.integer(lcal),
              as.double(error),
              as.integer(plind),
              as.integer(calind),
              as.integer(ind),
              as.integer(ppm),
              as.integer(uniq)
       )
      
return(res)
  if(res[[8]]>0)
  {
     test <- list(plind=res[[6]][1:res[[8]]]+1,calind=res[[7]][1:res[[8]]]+1)
  }
  else{
    test <- list(plind=NULL,calind=NULL)
  }
   return(test)
}
  
  
  

#---------------------------------------------------------------------
# Needed for internal calib to list.
#counts all matchings.


getaccmut <- function(peaklist , calib , modif , error = 500 , ppm=TRUE)
{
  #names(peaklist) <- 1:length( peaklist )
  #names(calib) <- 1:length( calib )
  #err2 <- 0.01 # if peaks are closer together then merge them.
    if(length(peaklist)>length(calib))
    {}else{
      tmp <- calib
      calib <- peaklist
      peaklist <- tmp
    }
    plt <- peaklist
    calt <- calib
  
    for(x in 1:length(modif))
    {
      plt<-c(plt,peaklist+modif[x],peaklist-modif[x])
    }
  plt<-plt[plt>0]
  plt <- sort(plt)
  #calt<-calt[calt>0]
  #calt <- sort(calt)

  #ttt <- (which(diff(plt)<err2)+1)
  #if(length(ttt)!=0)
  #{
  #	plt<-plt[-ttt]
  #}
  #remove peaks with small diffs.
  #ttt<-(which(diff(calt)<err2)+1)
  #if(length(ttt)!=0)
  #{
  #  	calt <- calt[-ttt]#remove peaks with small diffs.
  #}
  #cat("lpl = ",length(plt)," lcal = ", length(calt),"\n")

   # return(list(plt,calt))
#    cat("plt",length(plt),"calt",length(calt),"error",error,"ppm",ppm,"\n")
  test <- getaccC(plt,calt,error=error,ppm=ppm,uniq=TRUE)
    if(is.null(test$plind))
      {
        return(test)
      }
    else
      {
        test <- list(plind=unique(test$plind)
                     ,calind=unique(test$calind)
                     )
        return(test)
      }
}


#counts only the closes  matching.
getacc<-function(peaklist,calib,error=500,ppm=TRUE)
  {
    #ppm if true error is given in ppm otherwise in Da.
    ret<-NULL
    plind<-NULL
    calind<-NULL
    peaklist <- unique(peaklist)
    calib <- unique(calib)
    mungler <- 1
    if(length(peaklist)<length(calib))
      {
        tmp <- calib
        calib <- peaklist
        peaklist <- tmp
        mungler<-2
      }
    #print(length(peaklist))
    #print(length(calib))
    for(x in 1:length(calib))
      {
                                        #calculates the deviation in ppms
        if(ppm){
          tmp <- abs( (peaklist-calib[x]) )/calib[x]*1e6
          ind<-which( tmp<error )
        }else{
          tmp<-abs( peaklist-calib[x] )
          ind<-which( tmp<error )
        }
        if(length(ind)>0)
          {
            mins<-tmp[ ind ]
            nind<-ind[ mins==min(mins) ]
            plind<-c( plind,nind )
            calind<-c( calind,x )
          }
      }
    if(mungler==1)
      {
        return(list(plind=plind,calind=calind))
      }
    else
      {
        return(list(plind=calind,calind=plind))
      }
  }
  
