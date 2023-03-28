plotCircular<-function(area1,area2=NULL,spokes=NULL,
                       scale=0.8,labels,stats=TRUE,dp=1,
                       clockwise=TRUE,spoke.col='black',lines=FALSE,
                       inlab = TRUE, str_cen = "", incex = 1,
                       inunit = "", inscale = 0.9, insrt = 0, centrecirc=0.03,
                       main="", xlab="", ylab="", pieces.col=c("white","gray"),
                       length=FALSE, legend=TRUE, pt = 0.8,
                       auto.legend=list(x="bottomright",fill=NULL, labels=NULL, title=""), ...){
  
  ## NB: need some serious argumnt checking
  ## NB2: can find list of available colours using colors()
  
  ## No legend if only one variable and check area vars same length
  if (is.null(area2)) {
    legend <- FALSE
  } else {
    if (length(area1)!=length(area2))
      cat("Warning: length of", deparse(substitute(area1)),
          "and", deparse(substitute(area2)),"not equal\n")
  }
  
  ##print(pieces.col)
  density.1 <- density.2 <- 0
  if (pieces.col[1] != "white") density.1 <- NA
  if (length(pieces.col)==2) {
    if (pieces.col[2] != "white") density.2 <- NA
  }
  
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  on.exit(par(op)) # restore graphic settings whenever function exits
  
  bins<-length(area1)
  clockstart=pi/2 # default clock start at 12 o'clock
  half<- 2*pi/(bins*2) # for moving text/spokes half-way round
  if (clockwise==TRUE) mult=-1 else mult=1
  
  ## First plot a circle (of radius 1) as a frame
  detail<-200 # number that controls graphical detail of cheeses
  circle<-matrix(nrow=detail+1,ncol=2,data=0)
  frac<-1/detail
  for (i in 1:(detail+1)){
    circle[i,1]<-1*cos(2*pi*i*frac)
    circle[i,2]<-1*sin(2*pi*i*frac)
  }
  plot(circle,type='l',col='black',bty='n',yaxt='n',main=main,
       xlab=xlab, ylab=ylab, xlim=c(-1,1),ylim=c(-1,1),xaxt='n', ...)
  
  text(0,0, str_cen)
  
  ## scale cheeses to their area
  aarea1<-area1
  if(is.null(area2)==FALSE){
    aarea2<-area2
  }
  if(length==F){
    aarea1<-sqrt(area1*12/pi)
    if(is.null(area2)==FALSE){
      aarea2<-sqrt(area2*12/pi)
    }
  }
  
  ## scale the area to the maximum multiplied by the user-defined scale
  ## draw the cheeses
  for (cheeseno in 1:bins){
    if(is.null(area2)==TRUE){
      scaled1<-scale*aarea1/max(aarea1)
      cheese<-matrix(nrow=102,ncol=4,data=0)
      start<-2*pi*((cheeseno-1)/bins)+clockstart
      frac<-1/100
      cheese[1,1]<-centrecirc*mult*cos((2*pi*frac/bins)+start)
      cheese[1,2]<-centrecirc*sin((2*pi*frac/bins)+start)
      cheese[102,1]<-centrecirc*mult*cos((2*pi*100*frac/bins)+start)
      cheese[102,2]<-centrecirc*sin((2*pi*100*frac/bins)+start)
      cheese[1,3]<-centrecirc*mult*cos((2*pi*frac/bins)+start)
      cheese[1,4]<-centrecirc*sin((2*pi*frac/bins)+start)
      cheese[102,3]<-centrecirc*mult*cos((2*pi*100*frac/bins)+start)
      cheese[102,4]<-centrecirc*sin((2*pi*100*frac/bins)+start)
      for (i in 1:100){
        cheese[i+1,1]<-mult*scaled1[cheeseno]*cos((2*pi*i*frac/bins)+start)
        cheese[i+1,2]<-scaled1[cheeseno]*sin((2*pi*i*frac/bins)+start)
        cheese[i+1,3]<-mult*scaled1[cheeseno]*inscale*cos((2*pi*i*frac/bins)+start)
        cheese[i+1,4]<-scaled1[cheeseno]*inscale*sin((2*pi*i*frac/bins)+start)
      }
      polygon(cheese,density=density.1,angle=0,lty=1,lwd=1,border="black", col=pieces.col[1]) 
      if(inlab == TRUE){text(cheese[51,3], cheese[51,4], labels = paste(round(area1[cheeseno],2),inunit), cex = incex, srt = insrt)}
    }
    
    ## plot with two segments #
    ## 1st pattern
    if(is.null(area2)==FALSE){
      allarea<-c(aarea1,aarea2)
      scaled1<-scale*aarea1/max(allarea)
      scaled2<-scale*aarea2/max(allarea)
      cheese1<-matrix(nrow=52,ncol=4,data=0)
      cheese2<-matrix(nrow=52,ncol=4,data=0)
      start<-2*pi*((cheeseno-1)/bins)+clockstart
      frac<-1/100
      
      
      ## centrecirc: do not start at c(0,0) to prevent a dense block
      cheese1[1,1]<-centrecirc*mult*cos((2*pi*0*frac/bins)+start)
      cheese1[1,2]<-centrecirc*sin((2*pi*0*frac/bins)+start)
      cheese1[52,1]<-centrecirc*mult*cos((2*pi*51*frac/bins)+start)
      cheese1[52,2]<-centrecirc*sin((2*pi*51*frac/bins)+start)
      
      cheese1[1,3]<-centrecirc*mult*cos((2*pi*0*frac/bins)+start)
      cheese1[1,4]<-centrecirc*sin((2*pi*0*frac/bins)+start)
      cheese1[52,3]<-centrecirc*mult*cos((2*pi*51*frac/bins)+start)
      cheese1[52,4]<-centrecirc*sin((2*pi*51*frac/bins)+start)
      
      for (i in 1:50){
        cheese1[i+1,1]<-mult*scaled1[cheeseno]*cos((2*pi*i*frac/bins)+start)
        cheese1[i+1,2]<-scaled1[cheeseno]*sin((2*pi*i*frac/bins)+start)
        cheese1[i+1,3]<-mult*inscale*scaled1[cheeseno]*cos((2*pi*i*frac/bins)+start)
        cheese1[i+1,4]<-scaled1[cheeseno]*inscale*sin((2*pi*i*frac/bins)+start)
      }
      polygon(cheese1,density=density.1,angle=0,lty=1,lwd=1,border="black", col=pieces.col[1]) 
      if(inlab == TRUE){text(cheese1[26,3], cheese1[26,4], labels = paste(round(area1[cheeseno],2), inunit), cex = incex, srt = insrt)}
      
      ## 2nd pattern
      start<-2*pi*((cheeseno-1)/bins)+clockstart
      cheese2[1,1]<-centrecirc*mult*cos((2*pi*50*frac/bins)+start)
      cheese2[1,2]<-centrecirc*sin((2*pi*50*frac/bins)+start)
      cheese2[52,1]<-centrecirc*mult*cos((2*pi*100*frac/bins)+start)
      cheese2[52,2]<-centrecirc*sin((2*pi*100*frac/bins)+start)
      
      cheese2[1,3]<-centrecirc*mult*cos((2*pi*50*frac/bins)+start)
      cheese2[1,4]<-centrecirc*sin((2*pi*50*frac/bins)+start)
      cheese2[52,3]<-centrecirc*mult*cos((2*pi*100*frac/bins)+start)
      cheese2[52,4]<-centrecirc*sin((2*pi*100*frac/bins)+start)
      
      for (i in 51:100){
        cheese2[i+1-50,1]<-mult*scaled2[cheeseno]*cos((2*pi*i*frac/bins)+start)
        cheese2[i+1-50,2]<-scaled2[cheeseno]*sin((2*pi*i*frac/bins)+start)
        cheese2[i+1-50,3]<-mult*inscale*scaled2[cheeseno]*cos((2*pi*i*frac/bins)+start)
        cheese2[i+1-50,4]<-scaled2[cheeseno]*inscale*sin((2*pi*i*frac/bins)+start)
      }
      polygon(cheese2,density=density.2,angle=0,lty=1,lwd=1,border="black", col=pieces.col[2]) 
      if(inlab == TRUE){text(cheese2[26,3], cheese2[26,4], labels = paste(round(area2[cheeseno],2), inunit), cex = incex, srt = insrt)}
    } 
  }
  ## add the text
  if (is.null(labels)==FALSE&stats==FALSE){
    for (cheeseno in 1:bins){
      x<-mult*pt*cos((2*pi*cheeseno/bins)+start+half)
      y<-pt*sin((2*pi*cheeseno/bins)+start+half)
      text(x,y,labels[cheeseno])
    }
  }
  
  ## add the labels with stats
  if (is.null(labels)==FALSE&stats==TRUE){
    clabel2<-formatC(area1, format="f", digits=dp) # convert to character
    for (cheeseno in 1:bins){
      x<-mult*pt*cos((2*pi*cheeseno/bins)+start+half)
      y<-ptsin((2*pi*cheeseno/bins)+start+half)
      label1<-labels[cheeseno]
      label<-paste(label1,"\n",clabel2[cheeseno])
      text(x,y,label)
    }
  }
  
  ## add spokes representing uncertainty
  if (is.null(spokes)==FALSE){
    scaleds<-scale*spokes/max(spokes)
    halfcheese<-(2*pi)/(bins*2);
    for (cheeseno in 1:bins){
      spokes<-matrix(data=0,nrow=2,ncol=2)
      spokes[1,1]<-centrecirc*mult*scaleds[cheeseno]*
        cos((2*pi*cheeseno/bins)+start+half)
      spokes[1,2]<-centrecirc*scaleds[cheeseno]*
        sin((2*pi*cheeseno/bins)+start+half)
      spokes[2,1]<-mult*scaleds[cheeseno]*cos((2*pi*cheeseno/bins)+start+half)
      spokes[2,2]<-scaleds[cheeseno]*sin((2*pi*cheeseno/bins)+start+half)
      lines(spokes,pch=0,type='l',col=spoke.col,lty=1,lwd=1.5) 
    }
  } # end of spokes
  
  ## add dotted lines to separate months
  if (lines==TRUE){
    halfcheese<-(2*pi)/(bins*2);
    for (cheeseno in 1:bins){
      breaks<-matrix(data=0,nrow=2,ncol=2)
      breaks[1,1]<-centrecirc*cos((2*pi*cheeseno/bins)+start)
      breaks[1,2]<-centrecirc*sin((2*pi*cheeseno/bins)+start)
      breaks[2,1]<-cos((2*pi*cheeseno/bins)+start)
      breaks[2,2]<-sin((2*pi*cheeseno/bins)+start)
      lines(breaks,pch=0,type='l',lty=3,lwd=1) 
    }
  } # end of lines
  
  ## add legend if set
  if (legend) {
    if (length(auto.legend$x)==0) {
      legend.x <- "bottomright"
    } else {
      legend.x <- auto.legend$x
    }
    
    ##if (length(auto.legend$title)==0) {
    ##  title <- NULL
    ##} else {
    ##  title <- auto.legend$title
    ##}
    
    if (length(auto.legend$labels)==0){
      labels <- c(deparse(substitute(area1)),
                  deparse(substitute(area2)))
    } else {
      labels <- auto.legend$labels
    }
    
    if (length(auto.legend$fill)==0){
      fill <- pieces.col
    } else {
      fill <- auto.legend$fill
    }
    
    ##print(legend.x)
    ##print(labels)
    ##print(fill)
    ##print(title)
    legend(x=legend.x, legend=labels, title=auto.legend$title, fill=fill, ...)
  }
  
}

### Univariate plots ###

lab = c("BM model (1 in 20)", "KLOS model (1 in 20)", "POT model (1 in 20)", "BM model (1 in 100)", "KLOS model (1 in 100)", "POT model (1 in 100)", "BM model (1 in 200)", "KLOS model (1 in 200)", "POT model (1 in 200)")
a_1 = c(5269.34, 5442.17, 5342.79, 5777.78, 6290.67, 5734.75, 5925.04, 6594.52, 5836.55)
a_2 = c(4533.15, 4614.68, 4720.78, 5712.09, 5829.51, 6342.96, 6215.26, 6366.87, 7183.29)

plotCircular(area1=a_1, area2=a_2, scale=0.7, pt = 0.82, str_cen = "MRD",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.07, inlab = TRUE, incex = 1.2 ,inscale = 0.85,
             main="Scenario comparison for the time series T_MRD", insrt = 5)

a_1 = c(49.53, 49.00, 54.94, 79.53, 70.02, 108.50, 96.98, 80.56, 148.20)
a_2 = c(29.72, 29.34, 28.95, 37.94, 37.15, 34.18, 41.45, 40.41, 41.45)

plotCircular(area1=a_1, area2=a_2, scale=0.76, pt = 0.82, str_cen = "MRS",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.07, inlab = TRUE, incex = 1 ,inscale = 0.83,
             main="Scenario comparison for the time series T_MRS", insrt = 5)


### Multivariate plots ###

lab = c("CWBM model (G) (1 in 20)", "GCB (C) model (1 in 20)", "GCB (F) model (1 in 20)", "CWBM (G) model (1 in 100)", "GCB (C) model (1 in 100)", "GCB (F) model (1 in 100)", "CWBM (G) model (1 in 200)", "GCB (C) model (1 in 200)", "GCB (F) model (1 in 200)")
a_1 = c(75.67, 75.67, 69.21, 188.72, 111.48, 111.30, 270.56, 135.77, 135.77)
a_2 = c(34.49, 53.13, 51.33, 40.55, 66.13, 64.36, 44.64, 72.96, 75.26)

plotCircular(area1=a_1, area2=a_2, scale=0.7, pt = 0.78, str_cen = "MRS (OR)",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.1, inlab = TRUE, incex = 0.9 ,inscale = 0.8,
             main="Scenario comparison for the time series T_MRS", insrt = 10)

a_1 = c(23.74, 20.99, 23.58, 24.52, 23.33, 23.23, 24.74, 24.09, 27.03)
a_2 = c(19.45, 18.25, 18.84, 20.99, 20.97, 21.89, 22.05, 22.03, 23.05)

plotCircular(area1=a_1, area2=a_2, scale=0.5, pt = 0.75, str_cen = "MRS (AND)",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.11, inlab = TRUE, incex = 1 ,inscale = 0.8,
             main="Scenario comparison for the time series T_MRS", insrt = 5)


a_1 = c(5682.65, 5682.65, 5651.26, 6054.54, 6054.54, 6054.54, 6186.18, 6186.19, 6186.18)
a_2 = c(5395.05, 5041.63, 5041.63, 6634.02, 6208.54, 6208.54, 6934.77, 6869.81, 6869.81)

plotCircular(area1=a_1, area2=a_2, scale=0.55, pt = 0.75, str_cen = "MRD (OR)",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.1, inlab = TRUE, incex = 1 ,inscale = 0.83,
             main="Scenario comparison for the time series T_MRD", insrt = 0)

a_1 = c(4048.76, 3699.13, 4087.05, 4537.42, 4176.42, 4597.58, 4721.31, 4365.90, 4779.14)
a_2 = c(2990.54, 2893.40, 2976.69, 3728.82, 3283.89, 3411.14, 3437.95, 3437.95, 3575.97)

plotCircular(area1=a_1, area2=a_2, scale=0.55, pt = 0.77, str_cen = "MRD (AND)",
             labels=lab, stats = FALSE, dp=0, lines=TRUE, pieces.col=c("green","red"), inunit = "",
             auto.legend=list(labels=c("P1","P2"), title="TH"),centrecirc = 0.11, inlab = TRUE, incex = 0.8 ,inscale = 0.83,
             main="Scenario comparison for the time series T_MRD", insrt = 5)

