#Copyright 2004, W. Wolski, all rights reserved.
mget.myobj <- function(object,attrn,...)
  {
    ##t Field Access
    ##- Acces fields in object of class myobj
    ##+ object : object of class myobj
    ##+ attrn : name of field (Attribute)
    ##e data(mv1)
    ##e mget(mv1)
    ##e mget(mv1,"info")

    allow <- attr(object,"allow")
    if(attrn %in% allow)
      {
          return(attr(object,attrn))
      }
    else
      {
        warning( attrn , " not in the attributes list :", join(allow,sep=" ") ,"\n")
      }
  }

"setParms<-.myobj" <- function(object,value)
  {
    ##t Field Access
    ##- Set attributes in object of class myobj
    ##+ object : object of class myobj
    ##+ value : a list where list names are attributes names.
    ##e data(mv1)
    ##e mv1
    ##e setParms(mv1) <- list(info="test")
    ##e mv1
    
    tmp<-value
    allow<-attr(object,"allow")
    for(u in names(tmp))
      {
        if(u %in% allow)
          {
            if(!(u %in% "mass"))
              {
                attr(object,u)<-tmp[[u]]
              }
            else
              {
                tp<-attributes(object)
                object<-sort(tmp[[u]])
                attributes(object)<-tp
                return(object)
              }
          }
        else
          {
            warning(u," are not a allowed attribute list:", print(object$allow) ,"\n")
          }
      }
    object
  }


info.myobj <- function(object,info,...)
  {
    ##t Info Acces
    ##- access to the info field of the massvector.
    ##+ object : massvector
    ##+ info : info character. If missing function returns the current info. If not missing function returns massvector with new info field content.
    ##e data(mv1)
    ##e info(mv1)
    ##e mv1<-info(mv1,"testname")
    if(!missing(info))
      {
        attr(object,"info")<-info
        return(object)
      }
    return(attr(object,"info"))
  }
