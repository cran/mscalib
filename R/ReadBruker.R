readBruker.massvectorlist <- function(object,folder,expname,project,filename="peaklist.xml",...)
{
  ##t Read Bruker XML File
  ##- Read Bruker XML peaklist.xml file format.
  ##+ object : object of class massvectorlist
  ##+ folder : path to folder with Bruker file structure.
  ##+ expname : name of the experiment. If missing the name are extracted from the folder path.
  ##+ project : name of the project.
  ##+ filename : default = "peaklist.xml"
  ##+ ... : further parameters
  
  if(missing(expname))
    {
      expname <- strsplit(folder,"/")
      expname <- expname[[1]]
      expname <- expname[length(expname)]
      object<-experiment(object,expname)
    }else{
      object<-experiment(object,expname)
    }
  if(!missing(project))
    {
      set(object,"project")<-project
    }

  require(XML)
  name <- character(0)
  temp <- dir(folder,pattern="*Ref$")
  for(x in 1:length(temp))
    {
      name <- c(name,paste(folder,"/",temp[x],"/pdata/1/",filename,sep=""))
    }

  for(x in 1:length(name))
    {
      object[[x]] <- readBruker(massvector(),name[x],temp[x])
      if(x%%10==0)
        cat(formatC(x,width=3)," ")
      if(x%%100==0)
        cat("\n")
    }
  cat("\nlength massvectorlist",length(object),"\n")
  return(object)
}


getcoorBruker<-function(infop)
  {
    ##t Coordinates on Sample support.
    ##- Gets the coordinates on sample support from the default infor parameter (e.g. \code{0_A1_1SRef}).
    ##+ infop : path to folder.
    X <- 1:16
    names(X) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
    nam<-infop
    nam <- sub("[0-9]_","",nam)
    nam <- sub("_1.+$","",nam)
    y <-as.numeric(sub("[A-Z]+","",nam))
    x <-sub("[0-9]+","",nam)
    tcoor <- c( X[ x ] , y )
    names(tcoor)<-c(x,y)
    tcoor
  }


readBruker.massvector<-function(object,filep,infop,...)
{
  ##t Read peaklist.xml
  ##- Read Bruker Daltonic peaklist.xml file format.
  ##+ object : object of class massvector.
  ##+ filep : path and filename of file to load.
  ##+ infop : info field. If missing the info field is extracted out of the file path.
  
  require(XML)
  mass <- NULL
  area <- NULL
  test<-function() { 
                                        #vars <- character(0) ;
    massf<-function(x,attrs)
      {
        mass <<- c(mass, xmlValue(x[[1]])); 
        NULL
      }
    areaf<-function(x,attr)
      {
        area<<-c(area,xmlValue(x[[1]]));
        NULL
      }
    list(mass = massf,area=areaf)
  }
##################################################
  tt  <-  strsplit(filep,"/")
  pkname <- tt[[1]] [(1:length(tt[[1]]))[tt[[1]]=="pdata"]-1 ]
  
  if(file.exists(filep))
    {
                                        # retrieve the names of the file
                                        #cat("before the trap",fileName,"\n")
      xmlTreeParse(filep,handler=test())
      if(length( mass )!=0)
        {
          mass<-as.double(mass)
          area<-as.double(area)
          mass<-cbind(mass,area)
          colnames(mass) <- c("mass","area")
          
          or<-order(mass[,1])
          mass[,1]<-mass[,1][or]
          mass[,2]<-mass[,2][or]
          rownames(mass)<-1:length(mass[,1])
          object <- peaks(object,mass)
          
          if(!missing(infop))
            {
              peaks <-info(object,infop)
            }
          else
            {
             peaks <- info(object,pkname)
            }
          mcor <- getcoorBruker(pkname)
          setParms(peaks)<-list(tcoor=mcor)
        }
      else
        {
          mcor <- getcoorBruker(pkname)
          setParms(object)<-list(tcoor=mcor)
          peaks <- info(object,pkname)
        }
      return(peaks)
    }
  else
    {
      warning(paste("file : ",filep,"not fount!\n" ))
      mcor <- getcoorBruker(pkname)
      setParms(object)<-list(tcoor=mcor)
      return(info(object,pkname))
    }
}
