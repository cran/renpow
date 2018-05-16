# power electronics functions


rectifier <- function(v.t,full=FALSE){
 x <- v.t
 ny <- length(x$y)
 # for practical circuit
 if(full == TRUE){
  x$y <- abs(x$y)
  yavg <- 2*x$ym/pi
 } else {
  for(i in 1:ny){
   if (x$y[i] <= 0) x$y[i]<-0
  }
  yavg <- x$ym/pi
 }
 return(list(w=x$w,t=x$t,nw=x$nw,ym=x$ym,ang=x$ang,y=x$y,yrms=x$yrms,yavg=yavg))
}

ac.plot.rect <- function(V.t,v.lab="v(t)",v.units="V",y.lab="v(t)[V]",rms=FALSE){
 x <- V.t
 # pad max with margin of 20%
 ymax <- 1.2*max(x$y); tmax <- max(x$t) 
 # plot graph
 matplot(x$t,x$y,type="l", ylim=c(-ymax,ymax), xlab="t[sec]",
         ylab=y.lab,lty=1:x$nw,col=1,lwd=1.5)
 abline(v=0,h=0,lty=1,col='gray')
 horiz.lab(x$nw,x$ym,tmax,ymax,v.units,x$yrms,rms)
 wave.leg(x$nw,x$ang,v.lab,x$ym,x$w,v.units)
 abline(h=x$yavg,lty=2);text(0.001,x$yavg,"Avg",cex=0.7)

}


inverter <- function(x){

 w <- 2*pi*x$f; Tw <- 1/x$f
 t<- seq(0,x$nc*Tw,Tw/100); nt<- length(t)
 pad <- 1.5; hifreq.carrier <- 15
 vm <- cos(w*t)
 vc <- cos(hifreq.carrier*w*t)
 v <- cbind(vm,vc)
 wd=7; ht=7; cex1=0.8 
 panels(wd,ht,2,1,pty='m')

 matplot(t,v,type="l",col=1,ylim=c(-pad,pad),xlab="Time(s)",ylab="V (V)")
 abline(h=0,col='gray')
 legend('top',legend=c('Vmodulating','Vcarrier'),lty=1:2,col=1,cex=0.7) 
 vpwm <- array()
 for(i in 1:nt){
  if(vm[i]>vc[i]) vpwm[i] <- 1 else vpwm[i] <- -1
 }
 vpwm.S1.S2 <- vpwm; vpwm.S3.S4 <- -vpwm

 #plot(t,vpwm.S1.S2,type="l")
 #plot(t,vpwm.S3.S4,type="l")

 vout <- array()
 for(i in 1:nt){
  if(vpwm.S1.S2[i] > 0) vout[i] <- x$vin else vout[i] <- -x$vin
 }

 vout.filt <- x$vin*cos(w*t)
 vof <- cbind(vout.filt,vout)
 matplot(t,vof,type="l", xlab="Time(s)",col=1, ylab="Vout (V)",
        ylim=c(pad*min(vof),pad*max(vof)),lwd=c(2,1))
 abline(h=0,col='gray')
 legend('top',legend=c('Vout filtered','Vout'),lty=1:2,col=1,lwd=c(2,1),cex=0.7) 

}
 


