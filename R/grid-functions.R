# grid or electrical power systems functions

infinite.bus <- function(x){

 E <- x$E; nE<- length(E)
 delta <- x$delta; nd <- length(delta)

 P <- matrix(nrow=nd,ncol=nE); Q <- P
 I.p.m <- P; I.p.a <- P; cross <- array()

 for(kE in 1:nE){

 I.r.re <- E[kE]*sin(delta*pi/180)/x$X
 I.r.im <-(x$V-E[kE]*cos(delta*pi/180))/x$X
 #I.p.m <- I.r.re; I.p.a <- I.r.re
 for(i in 1:nd){ 
  I.p <- polar(c(I.r.re[i],I.r.im[i])) 
  I.p.m[i,kE] <- I.p[1]
  I.p.a[i,kE] <- I.p[2]
 }  
  
 P[,kE] <- ((E[kE]*x$V)/x$X)*sin(delta*pi/180)
 Q[,kE] <- (x$V/x$X)*(E[kE]*cos(delta*pi/180)-x$V)

 if(length(which(Q[,kE]>0))==0) cr <- 1 else cr <- which(Q[,kE]>0)[length(which(Q[,kE]>0))]
 cross[kE] <- delta[cr]
 }

 delta.lim <- c(delta[1],delta[nd])
 panels(6,6,2,2,pty="m")
 xlabel <- as.expression(bquote(delta*"["*degree*"]"))
 
 leglabels <- array()
 for(kE in 1:nE){ 
   leglabels[kE] <- as.expression(bquote("E="*.(E[kE])*"kV"*"Cr="*.(cross[kE])*degree))
 }
 
 matplot(delta,P,type="l",xlab=xlabel,ylab="Real power P[MW]",
         col=1,lty=1:nE,xlim=delta.lim)
 legend('topleft',legend=leglabels,lty=1:nE,cex=0.8)
 abline(v=cross,lty=2,col='gray')

 matplot(delta,Q,type="l",xlab=xlabel,ylab="Reactive power Q[MVAR]",
         col=1,lty=1:nE,xlim=delta.lim)
 legend('bottomleft',legend=leglabels,lty=1:nE,cex=0.8)
 abline(h=0,lty=2,col='gray'); abline(v=cross,lty=2,col='gray')

 matplot(delta,I.p.m,type="l",xlab=xlabel,ylab="Mag Current[kA]",
         col=1,lty=1:nE,xlim=delta.lim)
 legend('topleft',legend=leglabels,lty=1:nE,cex=0.8)
 abline(v=cross,lty=2,col='gray')

 matplot(delta,I.p.a,type="l",xlab=xlabel,ylab=as.expression(bquote("Angle current [ "*degree*" ]")),
         col=1,lty=1:nE,xlim=delta.lim)
 legend('bottomright',legend=leglabels,lty=1:nE,cex=0.8)
 abline(h=0,lty=2,col='gray'); abline(v=cross,lty=2,col='gray')
 

 #return(list(cross=cross))
}

# ---------------- demand ERCOT -------------

seg.ts <- function(Xts,dh0,dhf,var){ 
 X <- Xts
 x.t <- ts(round(X[,2:10],2),start=c(2010,1,1,1),end=c(2011,1,1,0),frequency=365*24)
 X[,1]<- as.character(X[,1])
 st <- which(X[,1]==dh0)
 ed <- which(X[,1]==dhf)
 x.t.seg <- x.t[st:ed,var]
 return(list(x.t=x.t,x.t.seg=x.t.seg))
}

series.plot <- function(x.t){
 ymax <- 1.2*max(x.t)
 ymin <- 0.8*min(x.t)
 ns <- dim(x.t)[2]
 lty.set <- c(1:ns)
 leg.lab <- colnames(x.t)
 plot.ts(x.t[,1],type="l",ylab="Demand (MW)",ylim=c(ymin,ymax),
         lty=lty.set[1],lwd=2)
 if (ns >= 2){
  for(i in 2:ns){
   lines(x.t[,i],type="l",lty=lty.set[i],lwd=2)
  }
 } 
  legend('top',legend=leg.lab,cex=0.7,col=1, lty=lty.set,horiz=T,adj = c(-0.0,0.1))
}

daily.min.max <- function(x.t.y,ylabel,inst.cap,inst.lab){ 
 ncap <- length(inst.cap)
 acc.inst <- array()
 for(i in 1:ncap) acc.inst[i] <- sum(inst.cap[1:i]) 
 x.d.max <- array();x.d.min <- array()
 for(i in 1:365){
  h0 <- (i-1)*24+1; h24 <- (h0)+23
  x.d.max[i] <- max(x.t.y[h0:h24]) 
  x.d.min[i] <- min(x.t.y[h0:h24]) 
 }
 maxlab <- paste(ylabel,"Daily Max")
 minlab <- paste(ylabel,"Daily Min")
 x.d.max.min <- ts(cbind(x.d.max,x.d.min),
      frequency=1,names=c(maxlab,minlab))
 series.plot(x.d.max.min)
 first <- 20
 offset <- 0.02* max(x.d.max)
 for(i in 1:ncap){
  abline(h=acc.inst[i],lty=i+1);text(first,acc.inst[i]+offset,inst.lab[i],cex=0.8,font=2)
 }
}

load.duration <- function(x.t.y, inst.cap,inst.lab){
 load <- sort(as.numeric(x.t.y),decreasing=T)
 if(length(load) >(24*365)) load <- load[1:(24*365)] 
 hours <- seq(1,24*365,1); hrs.max <- max(hours)

 ncap <- length(inst.cap)
 acc.inst <- array()
 for(i in 1:ncap) acc.inst[i] <- sum(inst.cap[1:i])  

 max.load <- max(load)
 ymax <- 1.1*max(load,inst.cap[ncap]) 
 plot(hours,load,type="l",ylim=c(0,ymax),ylab="Demand(MW)",xlab="Hours/year",lwd=2)
 abline(v=0, col='gray');abline(h=0, col='gray')

 acc.load <- acc.inst;acc.load[ncap] <- max.load; acc.hrs <- array()
 deno <- array();nume <- array();CF <- array()
 for(i in 1:ncap){
  acc.hrs[i] <- max(hours[which(load>=acc.load[i])])
  if(i<ncap){
   abline(h=acc.load[i],lty=2);abline(v=acc.hrs[i],lty=2)
   offx <-0.03*hrs.max;offy <- 0.015*max.load 
   text(acc.hrs[i]+offx,-offy,acc.hrs[i],cex=0.7)
  }
  deno[i] <- inst.cap[i]*hrs.max
  if(i==1) {hrs.fin <- hrs.max; bott <- 0
  } else {hrs.fin <- acc.hrs[i-1];bott<- acc.load[i-1]}
  nume[i] <- acc.load[i]*acc.hrs[i] + sum(load[acc.hrs[i]:hrs.fin])
  nume[i] <- nume[i] - bott*hrs.fin
  xx <- c(1:acc.hrs[i],(acc.hrs[i]+1):hrs.fin)
  yy <- c(rep(acc.load[i],acc.hrs[i]),load[(acc.hrs[i]+1):hrs.fin])
  nyy <- length(yy)
  CF[i] <- round(nume[i]/deno[i],2)
  polygon(c(xx,rev(xx)), c(yy,rep(bott,nyy)),density=10*i)
  offset <- 0.02* max.load
  abline(h=acc.inst[i],lty=2);text(8500,acc.inst[i]+offset,inst.lab[i],cex=0.8)
  text(8500,acc.inst[i]-2*offset,paste("CF=",CF[i],sep=''),cex=0.8)
 }
 energy <- "TWh"; possible <- round(deno/10^6,2); utilization <- round(nume/10^6,2)
 return(list(Hours=acc.hrs,energy=energy, possible=possible,utilization=utilization,CF=CF))
}



