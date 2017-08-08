#' @import data.table


# Inspired from http://gettinggeneticsdone.blogspot.com/2013/06/customize-rprofile.html 


### Create a new invisible environment for all the functions to go in so it doesn't clutter your workspace.

showNA <- function(x) {
  ## 15.6.15 als data.frame
  ## 7.7.15 apply statt sapply damit auch mit matrix funzend
  resi = apply(x,2, function(y) sum(is.na(y)))
  resi2 = data.frame(var = names(resi), NAs = as.vector(resi), vals = nrow(x)-as.vector(resi))
  resi2
  # if(is.data.table(x)) {
  # setDT(resi2)
  # return(resi2)} else return(resi2)
  
}

### Show the first an last  rows  of a data frame or matrix -> 3.8.17 provided as independent function
# ht <- function ( d, myrows=5 ) 
# { ## updated 11.3. to show all if dim 1 smaller than myrows*2
#   rows2show = min(dim(d)[1],myrows)
#   if(dim(d)[1] <= 2*rows2show) return(d) 
#   rbind ( head ( d ,  rows2show ), tail ( d ,  rows2show ))
# }


### matching und dabeiaufpassen, dass matchvariable unique ist
match_hk = function(x, y, testunique =T, makeunique = F,importcol = NULL, ...) {
  ##150122 makeunique = F statt T, na.omit bei duplicated y, fehlenden ok fall includiert
  ##160328 check auf gleiche laenge x und improtcol
  ##160616 match hk zeigt die duplikated zeilen statt mytabl falls ein Fehler kommt
  #   x = transkripte_eqtl$nuid
  #   y = ilmnAnnot013$nuid
  
  yname = deparse(substitute(y))
  
  # 150119 unique check auf schnelles duplicated umgestellt, auto makeuniuq
  if(testunique ==T){
    check = as.numeric(sum(duplicated(na.omit(y))))
    if(identical(check, 0)) return(match(x, y, incomparables=c(NA, NaN),...))   
    
    if(identical(check, 0)==F  & makeunique == F) {
      print(y[duplicated(y)])
      stop(paste(yname ,"ist nicht unique"))
    }
    
    if(identical(check, 0)==F  & makeunique == T) {
      
      ## try to make it nunique
      if(is.null(importcol)) stop("When asking for make unique, please provide vector with values to be imported")
      if(length(importcol) != length(y)) stop("When asking for make unique, please provide vector with values to be imported")
      
      matcher = unique(data.frame(index = y, importcol = importcol))
      matcher = matcher[ matcher$index %in% x,]
      matchercheck = as.numeric(sum(duplicated(na.omit(matcher$index))))
      if(identical(matchercheck, 0)==F  ) {
        print(matcher[allDuplicatedEntries(matcher$index),])
        stop(paste(yname ,"ist nicht unique after trying to make index and importcol unique..."))
      }
      return(match(x, y, incomparables=c(NA, NaN),...))           
      indinfo[match(x, y, incomparables=c(NA, NaN)),id_prepro_ge1]
      
    }
    
    
    
  }
  if(testunique ==F)  return(match(x, y, incomparables=c(NA, NaN),...))            
}

## Returns a logical vector TRUE for elements of X not in Y
`%nin%` <- function(x, y) !(x %in% y) 


### Plot a non-proportional 2-Way, 3-Way or 4-Way Venn Diagram 
### http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/vennDia.R
### NOTE from t.girke: This script has been replaced by overLapper.R, which provides much more  powerful and scalable utilities. The new overLapper.R script is available at:   http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn


## Define venndiagram function and three wrappers
venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE, ...) {
  ## Remove duplicates and NA fields in x, y, z and w
  if(unique==T) {
    x <- unique(x); x <- as.vector(na.omit(x))
    y <- unique(y); y <- as.vector(na.omit(y))
    if(!missing("z")) {
      z <- unique(z); z <- as.vector(na.omit(z))
    }
    if(!missing("w")) {
      w <- unique(w); w <- as.vector(na.omit(w))
    }
  }
  
  ## Check valid type selection
  if(!type %in% c("2", "2map", "3", "3map", "4", "4map", "4el", "4elmap")) {
    return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
  }
  
  ## Plot a 2-way venn diagram
  if(type=="2") {
    # Define ovelap queries
    q1 <- x[x %in% y]
    q2 <- x[!x %in% y]
    q3 <- y[!y %in% x]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep="")} else {mysub <- ""}
    if(plot==T) {
      ## Plot the 2-way venn diagram
      symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 2-way mapping venn diagram
  if(type=="2map") {
    olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
    symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
  }
  
  ## Plot a 3-way venn diagram
  if(type=="3") {
    ## Define ovelap queries
    q1 <- x[x %in% y & x %in% z]
    q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
    q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
    q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
    q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
    q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
    q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep="")} else { mysub <- "" }
    if(plot==T) {
      ## Plot the 3-way venn diagram
      symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 3-way mapping venn diagram
  if(type=="3map") {
    olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
    symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
  }
  
  ## Overlap queries for 4-way venn diagram
  if(type=="4" | type=="4el" | type=="4elmap") {
    ## Define ovelap queries
    xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
    q1 <- xy[xy %in% zw]
    q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
    q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
    q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
    q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
    q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
    q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
    q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w]
    q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
    q10 <- x[!x %in% c(y,z,w)]
    q11 <- y[!y %in% c(x,z,w)]
    q12 <- z[!z %in% c(x,y,w)]
    q13 <- w[!w %in% c(x,y,z)]
    q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
    q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)
    
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep="") } else { mysub <- "" }
    
    ## Plot 4-way venn diagram as circles
    if(plot==T & type=="4") {
      symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
      text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
      text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
    }
    
    ## Plot 4-way venn diagram as ellipses
    if(plot==T & (type=="4el" | type=="4elmap")) {
      olDF <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=countDF$count)
      ## Plot ellipse
      plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
        angles <- (0:segments) * 2 * pi/segments
        rotate <- rotate*pi/180
        ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
        ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
        ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
        plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
      }
      ## Plot ellipse as 4-way venn diagram
      ellipseVenn <- function(lines=lines, olDF, title=title, labels=labels, sub=mysub, main, lcol=lcol, tcex=1.3, ...) {
        split.screen(c(1,1))
        plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=title, sub=mysub, ...)
        screen(1, new=FALSE)
        plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, ...)
        screen(1, new=FALSE)
        plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, ...)
        screen(1, new=FALSE)
        plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, ...)
        text(olDF[1:15,1], olDF[1:15,2], olDF[1:15,3], col=tcol, ...)
        text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels, col=lcol, ...)
        close.screen(all=TRUE)
      }
      ## Plot 4-way ellipse venn diagram
      if(type=="4el") {
        ellipseVenn(olDF=olDF, lcol=lcol, lines=lines, labels=labels, title=title, ...)
      }
      
      ## Plot 4-way ellipse mapping venn diagram
      if(type=="4elmap") {
        olDFdebug <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=paste("q", 1:15, sep=""), ...)
        ellipseVenn(olDF=olDFdebug, lcol=lcol, lines=lines, labels=paste(labels, "=", c("x","y","z","w")), title="Mapping Venn Diagram", ...)
      }
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 4-way circle mapping venn diagram
  if(type=="4map") {
    olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
    symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
    text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
  }
  
}


venn2 = function(x1,y1, mytitle="2-Way Venn Diagram", mylabels = NA, plotte =T)
{
  # 28/2/13 plotte par    
  # 150119 vector check
  if(all(is.vector(x1) | is.factor(x1),is.vector(y1)|is.factor(y1))==F) stop("All input data must be vectors...")
  if(is.na(mylabels[1])) mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)))
  qlist <- venndiagram(x=x1, y=y1, unique=T, title=mytitle, labels= mylabels, plot=plotte, lines=c(2,3), lcol=c(2,3), tcol=c(1,1,1), lwd=3, cex=1.3, printsub=T, type="2")
  qlist
}


venn3 = function(x1,y1,z1, mytitle="3-Way Venn Diagram", mylabels = NA,  plotte =T)
{
  # 28/2/13 plotte par
  # 150119 vector check
  if(all(is.vector(x1)|is.factor(x1),is.vector(y1)|is.factor(y1),is.vector(z1)|is.factor(z1))==F) stop("All input data must be vectors...")
  
  if(is.na(mylabels[1])) mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)), deparse(substitute(z1)))
  qlist <- venndiagram(x=x1, y=y1, z=z1, unique=T, title=mytitle, labels= mylabels, plot=plotte, lines=c(2,3,4), lcol=c(2,3,4), tcol=c(1,1,1,1,1,1,1), lwd=3, cex=1.3, printsub=T, type="3")
  qlist
}

venn4 = function(x1,y1,z1,w1, mytitle="4-Way Venn Diagram", mylabels = NA,  plotte =T)
{
  # 13/07/03
  # 150119 vector check
  if(all(is.vector(x1)|is.factor(x1),is.vector(y1)|is.factor(y1),is.vector(z1)|is.factor(z1),is.vector(w1)|is.factor(w1))==F) stop("All input data must be vectors...")
  
  if(is.na(mylabels[1])) mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)), deparse(substitute(z1)), deparse(substitute(w1)))
  qlist <- venndiagram(x=x1, y=y1, z=z1, w=w1, unique=T, title=mytitle, labels=mylabels,
                       plot=plotte, lines=c(2,3,4,6), lcol=c(2,3,4,6), tcol=1, lwd=3, cex=1, printsub=T, type="4el")
  qlist
}



# ### Show the first 5 rows and first 5 columns of a data frame or matrix --> 3.8.17 provided as independent function
# hh = function ( d, mydims=5 ) {
#   # 29.1.15 data.table included
#   if("data.table" %in% class(d)) { 
#     d[1 : min(dim(d)[1],mydims), names(d)[ 1 : min(dim(d)[2],mydims)], with = F]
#   } else   d [ 1 : min(dim(d)[1],mydims) , 1 : min(dim(d)[2],mydims) , drop =F]
# }



### write delimeted tab
write.delim = function(x, y, writeColnames=T,writeRownames = F, createDir = F, ...) {
  ## create Dir option hinyugefuegt
  # 8.2. rownameparameter hinzugefuegt
  if(createDir ==T ){
    oldwd = getwd()
    pathname = unlist(stringr::str_split(y, pattern="/"))
    pathname = paste(pathname[1:(length(pathname)-1)], collapse="/")
    vortest = try(setwd(pathname), silent=T)
    test = identical(vortest , oldwd)
    setwd(oldwd)    
    if (test == F) {
      dir.create(pathname,recursive=T)
      message("\n...created directory ", pathname)
      
    } else message("\n... directory ", pathname, " already exists...")
    
  }
  
  write.table(x, y, quote=F, col.names=writeColnames, row.names= writeRownames, sep="\t")
}


### table categs
mytable = function (x, mydigits = 1, doprint = F, do_science_output = F) {
  res = c()
  res$num <- table(x, useNA="always")   
  res$output = data.frame(res$num)
  names(res$output) = c("category" , "freq")
  res$output$percent = res$num /sum(res$num)
  res$output$observed = paste0(format(res$output$freq, big.mark=",", scientific=do_science_output), " (", round(res$output$percent*100,mydigits), "%)")
  res$output$observed = stringr::str_trim(res$output$observed)
  zeilennamen = as.character(res$output$category)
  zeilennamen[is.na(zeilennamen)] = "NA"
  rownames(res$output) = zeilennamen
  res$output$category = zeilennamen
  if(doprint) print(res$output[,c('observed'), drop= F])
  res$output
  
}


### show allwahys NA when using xtabs
xtabs_hk = function(...) xtabs(... , exclude = NULL, na.action= "na.pass") # 12.6.15 zweites komma weggemacht, in .env verschoben

### Takes a dataframe and a column name, and moves that column to the front of the DF.
moveColFront <- function ( d = dataframe , colname = "colname" ) {
  ## 15.6.15 data.table auf setcolorder umgestellt
  ## multiple ohne warning
  stopifnot(all(colname %in% names(d)))
  index <- match ( colname , names ( d ))
  old_order <- 1:ncol(d)
  new_order <- c(index, old_order[old_order %nin% index])
  if(is.data.table(d)) setcolorder(d, new_order)
  if(is.data.table(d)==F)  d= d[,new_order]
  return(d)
}

### numbers as nice strings
huebsch = function(x, stellen =1) format(round(x,stellen), big.mark = ",")

### numbers as nice percentages
proz = function(x, stellen = 1) paste0(round(100*x, stellen), "%")



### Function for arranging ggplots. example:   multiplot(p1, p2, p3, p4, cols=2). It can take any number of plot objects as arguments, or if it can take a list of plot objects passed to plotlist.

multiplot <- function(..., plotlist=NULL,  cols=1, layout=NULL) {
  ## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #

  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## An T F index for alle duplicated entries, not only the duplicated one
allDuplicatedEntries <- function (vektor) {
  ## 150303 umgestellt auf datatable
  if(length(vektor)==0) return(0)

  vektab = data.table(myvektor = vektor, num = 1:length(vektor))
  duplicated_vals = vektab[duplicated(myvektor),myvektor]
  duplicated_entries = vektab[ myvektor %in% duplicated_vals]
  setkey(duplicated_entries, myvektor)
  duplicated_entries$num
  
}

### quickly-find-class-of-dataframe
showClassDF <- function(x) {
  ## 12.6.15 als data.frame
  resi = unlist(lapply(unclass(x),class))
  resi = data.frame(column = names(resi), class = as.vector(resi))
  resi$column[is.na(resi$column)] = "NA"
  rownames(resi) = as.character(resi$column)
  resi$column = NULL
  resi
} 


formateTimediff = function(timediff, mydigits = 3) paste0(format(unclass(timediff), digits = mydigits), " ", attr(timediff, "units"))
