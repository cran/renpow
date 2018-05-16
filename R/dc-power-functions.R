
resistor <- function(V,R){
 I <- V/R
 P <- abs(V*I)
 txt <- bquote(italic(R)*"="*.(R)*Omega)
 return(list(V=V,I=I,P=P,txt=txt))
}

diode <- function(V){
 # Diode model
 I0 = 10^(-9) # 1 nA
 q.kT = 38.9
 # I-V model
 I <- I0* (exp(q.kT*V) -1)
 P <- abs(V*I)
 txt <- ""
 return(list(V=V,I=I,P=P,txt=txt))
}

ivplane <- function(x,x0=FALSE,y0=FALSE){

 vlab <- expression("Voltage "*italic("V")*" [V]")
 ilab <- expression("Current "*italic("I")*" [A]")
 plab <- expression("Power "*italic("P")*" [W]")

 X <- list(x=x$V,y1=x$I,y2=x$P,xlab=vlab,y1lab=ilab, y2lab=plab,
           y1simb="I",y2simb="P",txt=x$txt)  
 plot2yaxis(X,x0,y0)
 mtext(X$txt,side=3,line=-1,cex=0.7)
}

vsource <- function(Voc,Rs){
 # x is Voc and Rs
 V <- seq(0,Voc,0.1)
 I <- (Voc -V)/Rs
 P <- V*I
 txt <- bquote(italic(R[s])*"="*.(Rs)*Omega*","*italic(V[oc])*"="*.(Voc)*"V")
 return(list(V=V,I=I,P=P,txt=txt))
}

isource <- function(Isc,Rp){
 # x is Isc and Rp
 I <- seq(0,Isc,0.1)
 V <- (Isc - I)*Rp
 P <- V*I
 txt <- bquote(italic(R[p])*"="*.(Rp)*Omega*","*italic(I[sc])*"="*.(Isc)*"A")
 return(list(V=V,I=I,P=P,txt=txt))
}

eff.pow <- function(x.eff.pow){
 x <- x.eff.pow
 # arg Rth,Voc
 RL <- x$Rth*seq(0.1,10,0.1)
 VL <- x$Voc*RL/(RL+x$Rth);I <- x$Voc/(RL+x$Rth)
 Pin<- x$Voc*I;Pout <- VL*I
 nup <- RL/(RL+x$Rth)

 X <- list(x=RL/x$Rth,y1=nup,y2=cbind(Pin,Pout),
          xlab=expression(italic(R[L])*"/"*italic(R[th])),y1lab="Efficiency", y2lab="Power [W]",
          y1simb=expression(eta[p]),y2simb=c("Pin","Pout"))
 plot2yaxis(X)
 abline(v=1,col='gray'); abline(h=0.5,col='gray')
 mtext(bquote(italic(R[th])*"="*.(x$Rth)*Omega*","*italic(V[oc])*"="*.(x$Voc)*"V"),side=3,line=-1,cex=0.7)
}

transient<- function(ys,tau,ylabel,yslabel){
 tmax <- 5*tau
 t = seq(0,tmax,0.01) 
 y.d <- ys*exp(-t/tau)
 y.c <- ys -y.d
 ymax <- max(y.d)

 plot(t,y.c,type="l",xlab="Time (s)",ylab=ylabel,xaxs="i",yaxs="i",
     ylim=c(0,1.2*ymax), cex.lab=0.8,cex.axis=0.8)
 lines(c(tau,tau),c(0,y.c[which(t==tau)]),lty=2) 
 lines(c(tau,0),c(y.c[which(t==tau)],y.c[which(t==tau)]),lty=2) 
 text(tau+0.05*tmax,0.05*ymax, expression(tau))
 text(0.05*tmax,0.63*ymax+0.05*ymax, "63%",cex=0.7)
 abline(h=ys,lty=2) 
 text(0.1*tmax,ymax+0.05*ymax, yslabel,cex=0.7)


 plot(t,y.d,type="l",xlab="Time (s)",ylab=ylabel,xaxs="i",yaxs="i",,
     ylim=c(0,1.2*ymax),cex.lab=0.8,cex.axis=0.8) 
 lines(c(tau,tau),c(0,y.d[which(t==tau)]),lty=2) 
 lines(c(tau,0),c(y.d[which(t==tau)],y.d[which(t==tau)]),lty=2) 
 text(tau+0.05*tmax,0.05*ymax, expression(tau))
 text(0.05*tmax,0.37*ymax+0.05*ymax, "37%",cex=0.7)


}

PVcell <- function(x.PVcell){
 # arg Isc Area Rs Rp Light
 # PV model
 x <- x.PVcell 
 Vd <- seq(-0.1,0.6,0.01); nd <- length(Vd) # volts
 I0 <- x$I0.A*10^-12*x$Area
 q.kT = 38.9
 Vt <- matrix(nrow=nd,ncol=length(x$Light))
 It <- matrix(nrow=nd,ncol=length(x$Light))
 pos <- array()
 # I-V model
 Isc <- x$Isc.A*10^-3*x$Area*x$Light
 Id <- I0*(exp(q.kT*Vd) -1)
 for(i in 1:length(x$Light)){
  It[,i] <- Isc[i] - Id - Vd/x$Rp
  Vt[,i] <- Vd - It[,i]*x$Rs
  pos[i] <- length(which(Vt[,i]<0))
 }
  st <- max(pos)
  #V <- round(Vt[-c(1:st),],3); I <- round(It[-c(1:st),],3)
  V <- round(Vt,3); I <- round(It,3)
  P <- round(V*I,3)
  txt <- bquote(italic(I[sc])*"="*.(Isc)*" A,"
         *italic(R[p])*"="*.(x$Rp)*Omega*", "*italic(R[s])*"="*.(x$Rs)*Omega)
 return(list(V=V,I=I,P=P,Light=x$Light,txt=txt))
}

PVcell.plot <- function(y.PVcell){

x <- y.PVcell

wd=7;ht=7;panels(wd,ht,2,1,pty="m")

nl <- length(x$Light)
ymax <- 1.2*max(x$I)
matplot(x$V,x$I,type="l",xlim=c(0,0.7),ylim=c(0,ymax),xaxs="i",yaxs="i",
            xlab="V (V)",ylab="I [A]",lty=1:nl, col=1,lwd=1)
legend("topright", paste(as.character(x$Light*100),"% Full Sun"),
     lty=1:nl,col=1,lwd=1,cex=0.8)

# power
ymax <- 1.2*max(x$P)
matplot(x$V,x$P,type="l", ylim=c(0,ymax), xlim=c(0,0.7),
        xaxs="i",yaxs="i", xlab="V (V)",ylab="P (W)",lty=1:nl, lwd=1,col=1)
legend("topright", paste(as.character(x$Light*100),"% Full Sun"),
     lty=1:nl,col=1,lwd=1,cex=0.8)
}

fuel.cell<- function(x.fcell){
 x <- x.fcell
 #arguments
 area.cm2 <- x$area.cm2; Rload.ohm <- x$Rload.ohm
 # voltage
 V <- seq(0,1.2,0.01)
 # parameters
 V0 <- 0.6; a <- 0.11; Jmax=2; Js <- 1.6; Ja <-0.25
 # current
 I <- area.cm2*Jmax/(1+exp((V-V0)/a))
 # ohmic approx
 I.o <- area.cm2*(V-0.85)/(-0.25)
 # power	
 P <- V*I
 VmaxP <- V[which(P==max(P))]
 # max power point

 par(mar = c(5,5,2,5))

 plot(V,I,type="l", xlab="Cell Voltage (V)",ylab="Current (A)",
     ylim=c(0,1.05*max(I)),xlim=c(0,1.1*max(V)),xaxs="i",yaxs="i",lwd=1.8)
 lines(V,I.o,lty=3,lwd=1.8)

 abline(v=VmaxP,lty=3,col='gray')
 text(0.7,1.02*Js*area.cm2,"Mass transfer loss")
 abline(h=Js*area.cm2,lty=3,col='gray')
 text(0.9,area.cm2*1.0,"Ohmic loss")
 abline(h=area.cm2*Ja,lty=3,col='gray')
 text(1.0,0.98*area.cm2*Ja,"Activation loss")

 leg1 <- c("Current","Power","Ohmic")

 if(is.na(Rload.ohm)==FALSE){ 
  Iload <- area.cm2*V/Rload.ohm
  lines(V,Iload, lty=4)
  leg.txt <- c(leg1,"Rload")
  leg.lty <- c(1:4) 
 } else {
  leg.txt <- leg1
  leg.lty <- c(1:3)
 }

 par(new=T)
 ymax <- 2.5*max(P)
 plot(V,P, axes=F, type="l", ylim=c(0,ymax), xlim=c(0,1.1*max(V)),xlab=NA,
        ylab=NA,lty=2,col=1,xaxs="i",yaxs="i",lwd=1.8)
 axis(side=4); mtext(side=4, line=3, "Power (W)")

 legend('topright',legend=leg.txt,lty=leg.lty)

} 




