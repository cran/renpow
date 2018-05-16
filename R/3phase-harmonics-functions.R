# 3-phase  functions

generator <- function(x){

 #S3p apparent in VA; Vl.rms line vol; pf 
 Il.rms <- x$S3p/(sqrt(3)*x$Vl.rms) # in A
 Il.p <- c(Il.rms,sign(x$lead.lag)*acos(x$pf)*180/pi)
 Vp.rms <- x$Vl.rms/sqrt(3) # in V
 Vp.p <- c(Vp.rms,0)
 Zs.p <- polar(x$Zs.r)
 E.r <- recta(mult.polar(Il.p,Zs.p))+recta(Vp.p)
 E.p <- polar(E.r)
 E.p.plot <- c(E.p[1],E.p[2])
 Vp.p.plot <- c(Vp.p[1],Vp.p[2])
 Il.p.plot <- c(Il.p[1]*10,Il.p[2])
 v.p <- list(E.p.plot,Vp.p.plot,Il.p.plot)
 v.units <- c("V","V","/10A")
 v.lab <- c("E","Vp","Iline")
 phasor.plot(v.p,v.lab,v.units,lty.p=c(1,2,3))
 return(list(E.p=E.p,Vp.p=Vp.p,Il.p=Il.p))
}

# power quality harmonics
harmonic <- function(x,harm.odd,lab.units){
 f=60 # frequency Hz
 w=2*pi*f #angular frequency (radians/sec)
 T= 1/f # period (sec)
 # time sequence for a number n of cycles
 n = 2; t = seq(0,n*T,0.001*T); nt = length(t)

 # number of fundamental waves
 nw <- length(x)
 Itot <- matrix(ncol=nw,nrow=nt)
 Isum <- Itot[,1]; Isum[] <- 0

 for(kw in 1:nw){
  # number of harmonics
  nharm<-length(harm.odd[[kw]])
  # sequence of harmonics
  Nharm <- seq(3,nharm*2+1,2)
  Im <- x[[kw]][1]; Iang <- x[[kw]][2] 
  # AC current
  # define arrays and matrices
  Iharm <- matrix(ncol=nharm+1,nrow=nt)
  Iharm[,1] <- Im*cos(w*t+Iang*pi/180)
  Itot[,kw] <- Iharm[,1]
  # loop thru harmonics, calc current and total
  for(k in 1:nharm){
   Iharm[,k+1] <- Im*harm.odd[[kw]][k]*cos(Nharm[k]*w*t+Nharm[k]*Iang*pi/180)
   Itot[,kw] <- Itot[,kw] + Iharm[,k+1]
 }
 THD <- round(sqrt(sum(harm.odd[[kw]]^2)),3)
 DF <- round(sqrt(1/(1+THD^2)),3)
 # plot graph of all harmonics
 y <- Iharm
 matplot(t,y,type="l", lty=1:(nharm+1),col=1, xlab="t [sec]",
         ylab=paste("Fundamental and harmonics",lab.units),
         ylim=1.2*c(min(y),max(y)))
 abline(v=0,h=0,lty=2)
 legend("topright",legend=c(1,Nharm),lty=1:(nharm+1),col=1,cex=0.7)
 mtext(paste("THD=",THD," DF=",DF),3,-1,cex=0.7)

 # plot graph of sum of all harmonics
 y <- cbind(Itot[,kw],Iharm[,1])
 matplot(t,y,type="l", col=1,lty=1:2, xlab="t [sec]",
         ylab=paste("Fundamental and total",lab.units),
         ylim=1.2*c(min(y),max(y)))
 abline(v=0,h=0,lty=2)
 mtext(paste("THD=",THD," DF=",DF),3,-1,cex=0.7)
 legend("topright",legend=c("Total","Fundam"),lty=1:2,col=1,cex=0.7)

 Isum <- Isum + Itot[,kw]
 }
 y <- Itot
 matplot(t,Itot,type="l",xlab="t [sec]", ylab=paste("Total per phase",lab.units),
         lty=1:nw,col=1, ylim=1.2*c(min(y),max(y)))
 mtext(paste("THD=",THD," DF=",DF),3,-1,cex=0.7)
 legend("topright",legend=c("a","b","c"),lty=1:nw,col=1,cex=0.7)
 plot(t,Isum,type="l",xlab="t [sec]", ylab=paste("Total neutral",lab.units),
      ylim=1.2*c(min(y),max(y)))
 mtext(paste("THD=",THD," DF=",DF),3,-1,cex=0.7)

return(list(t=t,Itot=Itot,Isum=Isum,THD=THD))
}




