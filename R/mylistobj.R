hist.mylistobj<-function(x,...)
  {
    ##t hist
    ##- hist
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ x : object of class mylistobj
    warning("hist not implemented for object of class mylistobj!\n")
  }

image.mylistobj <- function(x,...)
  {
    ##t image
    ##- image
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ x : object of class mylistobj
    warning("image not implemented for object of calss mylistobj!\n")
  }


mget.mylistobj <- function(object,attrn,...)
  {
    ##t Field Access
    ##- Access to the fields in the mylistobj
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ object : mylistobj
    ##+ attrn : The of the field to access. If missing the fields of the objects are returned.
    ##+ ... : further parameters
    ##v xxx : depends which field in the massvector are accessed.
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e class(res)
    ##e mget(mv1,"lenghtmv")
    ##e mget(mv1,"Coef.Intercept")
    
    allow <- object$allow
    if(attrn %in% allow)
      {
          return(object[[attrn]])
      }
    else
      {
        stop( attrn , " not in the attributes list :", join(allow,sep=" ") ,"\n")
      }
  }


info.mylistobj <- function(object,info , ...)
  {
    ##t Info Acces
    ##- access to the info field of the mylistobj
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ object : mylistobj
    ##+ info : info character. If missing function returns the current info. If not missing function returns mylistobj with new info field content.
    ##e data(mv1)
    ##e res <- getrecalib(mv1)
    ##e info(res)
    ##e res<-info(res,"testname")

    if(missing(info))
      mget(object,"info")
    else
      setParms(object)<-list(info=info)
  }

plot.mylistobj<-function(x,...)
  {
    ##t plot
    ##- plot
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ x : object of class mylistobj
    stop(" plot not implemented for object of class mylistobj!\n")
  }


"setParms<-.mylistobj"<-function(object,value)
  {
    ##t Field Access
    ##- Set attributes in object of class mylistobj
    ##d mylistobj provides basic functionality for objects implemented using list to store the attributes.
    ##+ object : object of class myobj
    ##+ value : a list where list names are attributes names.
    ##sa "setParms<-.myobj",mget.mylistobj
    ##e data(mv1)
    ##e res<-getrecalib(mv1)
    ##e setParms(res) <- list(info="test")
    ##e res

    tmp<-value
    allow<-object$allow
    for(x in names(tmp))
      {
        if(x %in% allow)
          {
            object[[x]]<-tmp[[x]]
          }
        else{
          stop(x," are not a allowed attribute list:",join(object$allow,sep=" "),"\n")
        }
      }
    object
  }
