###############
# finding peaks
###############




                                        #extends subset family.
#subset.list<-function (x, select, ...) 
#{
#    nl <- as.list(1:length(x))
#    names(nl) <- names(x)
#    vars <- eval(substitute(select), nl, parent.frame())
#    x[vars, drop = FALSE]
#}


                                        #cheks form objekt field
#has<-function(x,...)
#  {
#    UseMethod("has")
#  }

#has.default <- function(x,attrname,...)
#  {
#    if(!is.null(x$allow))
#      {
#        name <- x$allow
#      }
#    else if(!is.null(attr(x,"allow")))
#      {
#        name <- attr(x,"allow")
#      }
#    return(attrname %in% name)
#  }


#compare.massvectorlist<-function(object,mvl2,plot=TRUE,error=150,ppm=TRUE,...)
#  {
#    mvl1 <- object
#    rm(object)
#    ##t Compares massvectorlists
#    ##- Compares the masses in the massvectors. Returns the indices of the matching peaks given an measurment error.
#    ##- Plots the relative or absolute error of matchin peaks.
#    ##+ object : massvectorlist
#    ##+ mvl2 : massvector
#    ##+ plot : True - plot the relatve or absolute error. \code{FALSE} - no plotting.
#    ##+ error : size of the measurment error (default 150 ppm)
#    ##+ ppm : \code{TRUE} - relative error in parts per million, \code{FALSE} - absolute error.
#    ##v plind : indices of masses in mv1 matching to masses in mv2
#    ##v calind : indices of masses in  mv2 matchin to masses in mv1
#    
#    print("test")
#    if(!inherits(mvl2,"massvectorlist"))
#      stop("Second arg are not a massvectorlist but should be!\n")
#    res<-NULL
#    for(x in names(mvl1))
#      {
#        res <- rbind(res,compare(mvl1[[x]],mvl2[[x]],plot=FALSE,error=error,ppm=ppm))
#      }
#    rownames(res)<-names(mvl1)
#    if(plot)
#      pairs(res)
#    invisible(res)
#  }
