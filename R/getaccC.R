#Copyright 2004, W. Wolski, all rights reserved.
getaccC<- function(pl,cal,error=500,ppm=TRUE,uniq=FALSE)
{
  ##t Find Matching Masses
  ##- Returns the indices of masses matching to each other.
  ##+ pl : vector with masses
  ##+ cal : vector with masses
  ##+ error : error in Da or ppm.
  ##+ ppm : \code{TRUE} - then error has to be given as relative error in ppm. \code{FALSE} - error are absolute error in dalton.
  ##+ uniq : \code{TRUE} - return only closest match to mass in cal. \code{FALSE} - return all matches in the error range.
  ##e getaccC(1001:1010,1001:1010,error=300,ppm=TRUE,uniq=TRUE)
  lpl <- length(pl)
  lcal <- length(cal)
  #cat("lpl ",lpl, " lcal ",lcal," lmods ",lmods,"\n")
  tmp <- max(lpl,lcal) # there can not be more matches than tmp.
  plind <- rep(0,lpl*lcal)
  calind <- rep(0,lpl*lcal)
  ind <- 0
  if(length(pl)>0 & length(cal)>0)
  {
    if(uniq)
      {
      res <- .C("getaccU", 
                as.double(pl), 
                as.integer(lpl),  
                as.double(cal), 
                as.integer(lcal),
                as.double(error),
                as.integer(plind),
                as.integer(calind),
                as.integer(ind),
                as.integer(ppm),
                PACKAGE="mscalib"
                )
    }
    else
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
                  PACKAGE="mscalib"
                  )
      }
   }
   else{
    return(list(plind=NULL,calind=NULL))
   }     
  if(res[[8]]>0)
   {
      test <- list(plind=res[[6]][1:res[[8]]]+1,calind=res[[7]][1:res[[8]]]+1)
   }
  else
    {
      test <- list(plind=NULL,calind=NULL)
    }
  return(test)
}
  

