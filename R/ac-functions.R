# -------- calculate waves
waves <- function(x,f=60,nc=2){
 # angular frequency (radians/sec)
 w=2*pi*f
 # period (sec)
 T= 1/f 
 # time sequence for a number n of cycles
 t = seq(0,nc*T,0.01*T); nt = length(t)
 # number of waves
 nw <- length(x)
 # matrix with all waves, arrays with amplitude and phase
 y <- matrix(ncol=nw,nrow=nt); ym <- array(); ang <- array()
 # for each wave rename mag and phase and calc to fill the matrix
 for(i in 1:nw){
  ym[i] <- x[[i]][1]; ang[i] <- x[[i]][2]
  y[,i]=ym[i]*cos(w*t+ang[i]*pi/180)
 }
 yrms <- ym/sqrt(2)
 return(list(w=w,t=t,nw=nw,ym=ym,ang=ang,y=y,yrms=yrms))
}

# ------- plot horizontal lines with labels for magnitude and rms
horiz.lab <- function(nw,ym,tmax,ymax,units,yrms,rms){
 for(i in 1:nw){
  abline(h=ym[i]+0.005*ymax,lty=i,col='gray')
  text(0.1*tmax,ym[i]+0.02*ymax,paste(round(ym[i],1),units[i],sep=""),cex=0.7)
  if(rms==TRUE){
   abline(h=yrms[i]+0.005*ymax,lty=i,col='gray')
   text(0.2*tmax,yrms[i]+0.02*ymax,paste("RMS=",round(yrms[i],0),units[i],sep=""),cex=0.7)
  }
 }
}

# -------------- write out waves for legend
wave.leg <- function(nw,ang,lab,ym,w,units){
 wave <- array(); s.sym <- "+"
 ym <- round(ym,1); w<- round(w,0);ang <- round(ang,1) 	 
 for(i in 1:nw){
  if(ang[i]<0) s.sym <-"" else s.sym <- "+"
  wave[i] <- as.expression(bquote(.(lab[i])*"="*.(ym[i])*"cos("*.(w)*"t"*.(s.sym)* .(ang[i])*degree*")"*.(units[i])))
 } 
 legend('topright',legend=wave,lty=1:nw,col=1,bg='white',cex=0.7)
}

# --------plot AC cosine wave ---------------------
ac.plot <- function(v.t,v.lab="v(t)",v.units="V",y.lab="v(t)[V]",rms=FALSE){
 # pad max with margin of 20%
 ymax <- 1.2*max(v.t$y); tmax <- max(v.t$t) 
 # plot graph
 matplot(v.t$t,v.t$y,type="l", ylim=c(-ymax,ymax), xlab="t[sec]",
         ylab=y.lab,lty=1:v.t$nw,col=1,lwd=1.5)
 abline(v=0,h=0,lty=1,col='gray')
 horiz.lab(v.t$nw,v.t$ym,tmax,ymax,v.units,v.t$yrms,rms)
 wave.leg(v.t$nw,v.t$ang,v.lab,v.t$ym,v.t$w,v.units)
}


# -------------- write out waves for legend
phas.leg <- function(np,mag,ang,lab,units,lty.p){
  phas <-array()
  mag <- round(mag,1);ang <- round(ang,1)
  for(i in 1:np){
   phas[i] <- as.expression(bquote(.(lab[i])*"="*.(mag[i])*.(units[i])*","*.(ang[i])*degree))
 }
 legend('topleft',legend=phas,lty=lty.p,col=1,bg='white',cex=0.7)
}

arc <- function(mag,ang){
 rd <- mag
 ar <- seq(ang[1],ang[2],0.01); ar <- ar*pi/180
 nar <- length(ar); fnar <- round(0.9*nar,0)
 lines(rd*cos(ar),rd*sin(ar),col=1,lwd=1)
 arrows(x0=rd*cos(ar[fnar]),y0=rd*sin(ar[fnar]),
        x1=rd*cos(ar[nar]),y1=rd*sin(ar[nar]),
        col=1,lwd=1,lty=1,length=0.05)
}

rot.fig <- function(vp,v.lab="wt"){
 # number of phasors
 np <- length(vp)
 # magnitude
 mag <- array(); ang <- array()
 for(i in 1:np){
  mag[i] <- vp[[i]][1]; ang[i] <- vp[[i]][2]
 }
 angr <- ang*pi/180
 # max of all magnitude and pad by 20%
 xmax <- max(mag); x <- seq(-xmax,xmax,0.1)
 # input conversion to rect
 v.r <- list(); for(i in 1:np) v.r[[i]] <- recta(vp[[i]])
 # plot no axis
 plot(x,x,type="n",axes=FALSE,ann=FALSE)
 abline(h=0,v=0,col='gray')
 # draw circle
 ar <- seq(0,2*pi,0.01)
 lines(xmax*cos(ar),xmax*sin(ar),col=1,lwd=1)
 # plot arrows 
 for(i in 1:np){
  # plot phasor and label
  arrows(x0=0,y0=0,x1=v.r[[i]][1],y1=v.r[[i]][2],length=0.1,lty=i,lwd=2) 
  #text(v.r[[i]][1]+0.05*xmax,v.r[[i]][2]+0.05*xmax,v.lab[i],cex=0.7)
  arc(mag[i]/2,c(0,ang[i]))
  text((mag[i]/2)*cos(angr[i])+0.2*xmax,(mag[i]/2)*sin(angr[i])-0.1*xmax,v.lab[i])
 }
 #phas.leg(np,mag,ang,v.lab,v.units)
}

sinplot <- function(xlab,ylab){
 a <- seq(0,360,0.1); ar <- a*pi/180
 y <- 1*sin(ar)
 plot(a,y,type="l",ylab=ylab,xlab=xlab)
 abline(h=0, col='gray')
}

gridcir <- function(rmax){
 x <- rmax
 # max of all magnitude and pad by 20%
 xmax <- 1.2*x
 xd <- pretty(c(0,xmax)); nd <- length(xd)
 xs <- seq(-xd[nd],xd[nd],0.1)
 xl <- 1.01*c(min(xs),max(xs))
 # plot axis
 plot(xs,xs,type="n",xlab="Re",ylab="Im",xlim=xl,ylim=xl,cex.axis=0.6)
 abline(h=0,v=0,col='gray')
 # circles and labels
 rd <- xd; nr <- nd
 ar <- seq(0,2*pi,0.01)
 for(i in 1:nr){
  lines(rd[i]*cos(ar),rd[i]*sin(ar),col='gray')
  text(rd[i],-0.05*x,rd[i],cex=0.6, col='gray')
 }
 # radial lines and labels
 ang <- seq(0,360,30);ng <- length(ang); angr <- ang*pi/180
 for(i in 1:(ng-1)){
  coord <- c(rd[nr]*cos(angr[i]),rd[nr]*sin(angr[i]))
  lines(c(0,coord[1]),c(0,coord[2]),col='gray')
  coordtxt <- coord+sign(coord)*0.05*x
  text(coordtxt[1],coordtxt[2],as.expression(bquote(.(ang[i])*degree)),cex=0.6,col='gray')
 }
}

# ---------plot phasors -----------------
phasor.plot <- function(v.p,v.lab="V",v.units="V",lty.p=NULL){
 # number of phasors
 np <- length(v.p)
 if(is.null(lty.p)) lty.p <- c(1:np)
 # magnitude and phase
 mag <- array(); ang <- array()
 for(i in 1:np) {mag[i] <- v.p[[i]][1];ang[i] <- v.p[[i]][2]}
 # input conversion to rect
 v.r <- list(); for(i in 1:np) v.r[[i]] <- recta(v.p[[i]])
 xmax <- max(mag)
 gridcir(xmax)
 # plot arrows 
 for(i in 1:np){
  # plot phasor and label
  arrows(x0=0,y0=0,x1=v.r[[i]][1],y1=v.r[[i]][2],length=0.1,lty=lty.p[i],lwd=1.8) 
  xt <- v.r[[i]][1]; if(xt==0) xt <- 0.001
  yt <- v.r[[i]][2]; if(yt==0) yt <- 0.001
  text(xt+0.06*xmax*sign(xt),yt+0.06*xmax*sign(yt),v.lab[i],cex=0.7)
 }
 phas.leg(np,mag,ang,v.lab,v.units,lty.p)
}

 
# ------------ complex conversion
polar <- function(rec){
 x <- rec[1]; y <- rec[2]
 mag <- sqrt(x^2+y^2)
 if(is.nan(y/x)) ang <- 0 else ang <- atan(y/x)*180/pi
 return(round(c(mag,ang),3))
}

recta <- function(pol){
 mag <- pol[1]; ang <- pol[2]
 a <- mag * cos(ang*pi/180)
 b <- mag * sin(ang*pi/180)
 return(round(c(a,b),3))
}

mult.polar <- function(x1,x2){
 z.m <- x1[1]*x2[1]; z.ang <- x1[2]+x2[2]
 return(round(c(z.m,z.ang),3))
}

div.polar <- function(x1,x2){
 z.m <- x1[1]/x2[1]; z.ang <- x1[2]-x2[2]
 return(round(c(z.m,z.ang),3))
}

# ---------- admittance
admit <- function(Z.r){
 Z.p <- polar(Z.r)
 Y.p <- div.polar(c(1,0),Z.p)
 Y.r <- recta(Y.p)
 return(round(Y.r,3))
}

vector.phasor <- function(V,I){
 nV <- length(V); VI <- list()
 V.p <- matrix(c(Re(V),Im(V)),ncol=nV)
 for(i in 1:nV) VI[[i]] <- polar(V.p[i,])
 nI <- length(I)
 I.p <- matrix(c(Re(I),Im(I)),ncol=nI)
 for(i in 1:nI) VI[[i+nV]] <- polar(I.p[i,])
 return(list(V.p=V.p,I.p=I.p,VI=VI))
}
