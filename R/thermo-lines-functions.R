# draw lines in PV and Ts planes

path.lines <- function(x,plane='Pv',shade.between=FALSE,lab.cycle=FALSE,shade.cycle=FALSE){
 # in case a simple list, make it a list of lists
 if(any(sapply(x, class) == "list")==FALSE) x <- list(x)
 # number of lists in list
 nl <- length(x)

 if(plane == 'Pv'){
 xlab=~"Specific volume "*italic('v')*"("*m^3*"/kg)"
 ylab=~"Pressure "*italic('P')*"(bar)"
 }

 if(plane == 'Ts'){
 xlab=~"Specific entropy "*italic('s')*"((kJ/K)/kg)"
 ylab=~"Temperature "*italic('T')*"("*degree*"C)"
 }

 y <- list(); leg <- array(); linetype<- array()
 pr <- array(); nv <- array();nT <- array(); nP <- array()
 ns <- array();nxx <- array(); nyy <- array()

 vl <- matrix(nrow=2,ncol=nl)
 pl <- vl; lab <- vl
 
 for(i in 1:nl){
  y[[i]] <- path.calc(x[[i]])
  if(plane == 'Ts'){
    y[[i]]$xx <-  y[[i]]$s
    y[[i]]$yy <-  y[[i]]$T
  }
  if(plane == 'Pv'){
    y[[i]]$xx <-  y[[i]]$v
    y[[i]]$yy <-  y[[i]]$P
  }

  nv[i] <- length(y[[i]]$v)
  nP[i] <- length(y[[i]]$P)
  ns[i] <- length(y[[i]]$s)
  nT[i] <- length(y[[i]]$T)
  nxx[i] <-   length(y[[i]]$xx)
  nyy[i] <-   length(y[[i]]$yy)

  lab[,i] <- x[[i]]$lab
  
  dirpr <- sign(y[[i]]$v[1]- y[[i]]$v[nv[i]])
    if(dirpr==1)pr[i] <- " compr "
    if(dirpr==0)pr[i] <- " "
    if(dirpr==-1)pr[i] <- " expan "

  leg[i] <- as.expression(bquote(.(lab[1,i])*"-"*.(lab[2,i])*":"*.(x[[i]]$path)~.(pr[i])~
            .(round(y[[i]]$T[1],1))*degree*"C->"*.(round(y[[i]]$T[nT[i]],1))*
            degree*"C, W="*.(round(y[[i]]$WQtot[1],1))*"kJ/kg, Q="*
            .(round(y[[i]]$WQtot[2],1))*"kJ/kg"))

  vl[,i] <- c(min(y[[i]]$xx),max(y[[i]]$xx))
  pl[,i] <- c(min(y[[i]]$yy),max(y[[i]]$yy))
  linetype[i]=1
  if(nl>1){
   if(x[[i]]$path=='isotherm') linetype[i]=1
   if(x[[i]]$path=='adiabat' ) linetype[i]=2
  }
  else linetype[i]=1
 }
 #gastxt <- paste(names(y[[1]]$gas),"=",y[[1]]$gas,sep="",collapse=", ")
 xlim<- c(0.8*min(vl[1,]), 1.2*max(vl[2,]))
 ylim<- c(0.8*0, 1.2*max(pl[2,]))

 plot(y[[1]]$xx,y[[1]]$yy, xlab=xlab, ylab=ylab,
      xlim=xlim,ylim=ylim,lwd=2,cex.axis=0.8,cex.lab=0.9,type="n")

 for(i in 1:nl){

  lines(y[[i]]$xx,y[[i]]$yy,lwd=2,lty=linetype[i])
  x0 <- y[[i]]$xx[0.5*nxx[i]];y0=y[[i]]$yy[0.5*nyy[i]]
  x1 <- y[[i]]$xx[0.6*nxx[i]];y1=y[[i]]$yy[0.6*nyy[i]]
  arrows(x0,y0,x1,y1,length=0.15,lwd=2)

  if(is.null(x[[i]]$shade.under)==FALSE){
   if(x[[i]]$shade.under==TRUE)
   polygon(c(y[[i]]$xx,rev(y[[i]]$xx)), c(y[[i]]$yy,rep(ylim[1],nyy[i])),
           density=10)
  }

   pv <- matrix(nrow=2,ncol=2)
   # rows are coordinates of start and end
    pv[1,] <- c(y[[i]]$xx[1],y[[i]]$yy[1]) # start
    pv[2,] <- c(y[[i]]$xx[nxx[i]],y[[i]]$yy[nyy[i]]) # end
   if(pr[i]=='compr'){
    pv[c(2,1),] <- pv[c(1,2),]
    lab[,i] <- rev(lab[,i]) 
   }
   if(lab.cycle == FALSE){
    dirlab <- sign(pv[2,1]-pv[1,1]);if (dirlab==0) dirlab<-1
    xs <- 0.02*dirlab*(1+ 0.5*rnorm(1))
    ys <- 0.02*dirlab*(1+ 0.5*rnorm(1))
    text(pv[1,1]-xs*xlim[2],pv[1,2]+ys*ylim[2],lab[1,i],cex=0.7)
    text(pv[2,1]+xs*xlim[2],pv[2,2]+ys*ylim[2],lab[2,i],cex=0.7)
   }
   else{
    dirlab <- c(-1,-1,1,1)
    xs <- 0.02*dirlab[i]; ys <- 0.02*dirlab[i]
    text(pv[1,1]+xs*xlim[2],pv[1,2]+ys*ylim[2],lab[1,i],cex=0.7)
   }

 } # end loop lines

  if(shade.cycle==TRUE){
   xbox <- y[[1]]$xx; for(i in 2:nl) xbox <- c(xbox,y[[i]]$xx)
   ybox <- y[[1]]$yy; for(i in 2:nl) ybox <- c(ybox,y[[i]]$yy)
   polygon(c(xbox),c(ybox),density=10)
  }

  if(shade.between==TRUE){
   npts <- length(y[[1]]$xx)
   xleft <- c(y[[1]]$xx[1],y[[2]]$xx[1])
   yleft <- c(y[[1]]$yy[1],y[[2]]$yy[1])
   xright <- c(y[[1]]$xx[npts],y[[2]]$xx[npts])
   yright <- c(y[[1]]$yy[npts],y[[2]]$yy[npts])
   xbox <- c(xleft, y[[1]]$xx, xright, rev(y[[2]]$xx))
   ybox <- c(yleft, y[[1]]$yy, yright, rev(y[[2]]$yy))
   polygon(c(xbox),c(ybox),density=10)
  }


  legend("top", leg, cex=0.8)

  #pplt <- par("plt"); adjx <- pplt[1] / (pplt[2] - pplt[1])
  #mtext(gastxt, side=3, line=0, cex=0.7)
}


 
