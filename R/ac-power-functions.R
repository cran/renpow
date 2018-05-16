# ------------- calculate v,i and power
inst.pow.calc <- function(x,freq=60,nc=2){
 # time sequence for a number n of cycles
 t <- seq(0,nc/freq,0.01/freq); nt = length(t)
 # matrix v and i, array for power
 v.i <- matrix(ncol=2,nrow=nt)
 p <- array()
 # angular frequency (rad/sec)
 w=2*pi*freq
 vm <- c(x[[1]][1],x[[2]][1]);  ang <- c(x[[1]][2],x[[2]][2])
 vm[3] <- (1/2)*vm[1]*vm[2]/1000; ang[3] <- ang[1]-ang[2]
 pf <- cos(ang[3]*pi/180)
 pavg <- vm[3]*pf
 # v, i and p
  v.i[,1] <- vm[1]*cos(w*t+ang[1]*pi/180)
  v.i[,2] <- vm[2]*cos(w*t+ang[2]*pi/180)
  p <- v.i[,1]*v.i[,2]/1000
 return(list(w=w,vm=vm,ang=ang,pf=pf,pavg=pavg,t=t,v.i=v.i,p=p))
}

# ---------- legend for waves and power
inst.pow.leg <- function(ang,lab,ym,w,units,pf){
 wave <- array();s.sym <- "+"
 ym <- round(ym,1); w<- round(w,0);ang <- round(ang,1) 	 
 pf <- round(pf,2)	 
 for(i in 1:2){
  if(ang[i]<0) s.sym <-""
  wave[i] <- as.expression(bquote(.(lab[i])*"="*.(ym[i])*"cos("*.(w)*"t"*.(s.sym)*.(ang[i])*degree*")"*.(units[i])))
 }
 s.sym <- "+" 	 
 if((ang[1]+ang[2])<0) s.sym <-""
 wave[3] <- as.expression(bquote(.(lab[3])*"="*.(ym[3])*"("*.(pf)*"+"*
               "cos("*.(2*w)*"t"*.(s.sym)*.(ang[1]+ang[2])*degree*")"*.(units[3])))
 legend('topright',legend=wave,lty=1:3,col=1,bg='white',cex=0.7)
}

# --------plot AC cosine wave ---------------------
inst.pow.plot <- function(x,rms=FALSE,freq=60,nc=2){
 v.lab <- c("v(t)","i(t)","p(t)"); v.units <- c("V","A","kW")
 y.lab <- paste(v.lab,"[", v.units,"]",sep='')
 y1.lab <- paste(y.lab[1],"or",y.lab[2]); y2.lab <- y.lab[3]
 y <- inst.pow.calc(x)
 # pad max with margin of 20%
 ymax <- 1.2*max(y$vm[1:2]); tmax <- max(y$t)
 vrms <- y$vm[1:2]/sqrt(2)

 # plot graph v and i
 par(mar = c(5,5,2,5))
 matplot(y$t,y$v.i,type="l", ylim=c(-ymax,ymax), xlab="t[sec]",
        ylab=y1.lab,lty=1:2,col=1,lwd=1.7)
 abline(v=0,h=0,lty=1,col='gray')
 horiz.lab(2,y$vm,tmax,ymax,v.units,vrms,rms)

 # plot power right hand side axes
 par(new=T)
 pmax <- 1.4*(y$vm[3]+y$pavg)
 plot(y$t,y$p, axes=F, type="l", ylim=c(-pmax,pmax), xlab=NA,
        ylab=NA,lty=3,col=1,lwd=1.7)
 axis(side=4); mtext(side=4, line=3, y2.lab)
 abline(h=y$pavg+0.005*pmax,lty=3,col='gray')
 text(0.1*max(y$t),y$pavg+0.02*pmax,paste('Pavg=',round(y$pavg,2),v.units[3],sep=""),cex=0.7)

 # legend
 inst.pow.leg(y$ang,v.lab,y$vm,y$w,v.units,y$pf)
}

# -----------------
complex.pow.calc <- function(xc,dig=2,res=TRUE){
 theta <- xc[[3]]; thetar <- theta*pi/180
 V <- xc[[1]]; I <- xc[[2]] 
 S <- V*I/1000 
 P <- S*cos(thetar)
 Q <- S*sin(thetar)
 pf <- cos(thetar)
 y <- round(c(S,theta,P,Q,pf),dig)
 units<-c("kVA","kW","kVAR","A")
 prnt <- c(paste0("S=",y[1],units[1],",theta=",y[2],"deg"),
           paste0("P=",y[3],units[2],", Q=",y[4],units[3]),
           paste0("pf=",y[5])
          )
		  
 if(res==TRUE) cat(prnt,sep="\n")
 return(list(units=units, S=y[1], theta=y[2],P=y[3], Q=y[4], pf=y[5],prnt=prnt))
}


complex.pow.plot <- function(cp){
 # input is from complex.pow.calc
 # magnitude and phase
 mag <- cp$S; ang <- cp$theta; angr <- ang*pi/180
 # conversion to rect
 v.r <- c(cp$P,cp$Q)
 # max of all magnitude
 xmax <- mag 
 # plot axis
 gridcir(xmax)
 # plot arrows 
  # plot apparent
  arrows(x0=0,y0=0,x1=v.r[1],y1=v.r[2],length=0.1,lty=1,lwd=1.5) 
  text(v.r[1]+0.05*xmax,v.r[2]+0.05*xmax,"S",cex=0.7)
  # plot real
  arrows(x0=0,y0=0,x1=v.r[1],y1=0,length=0.1,lty=1,lwd=1) 
  text(v.r[1]-0.05*xmax,0+0.05*xmax,"P",cex=0.7)
  # plot imaginary
  arrows(x0=0,y0=0,x1=0,y1=v.r[2],length=0.1,lty=1,lwd=1) 
  text(0+0.05*xmax,v.r[2]-0.05*xmax,"Q",cex=0.7)
  # plot projection imag
  arrows(x0=v.r[1],y0=0,x1=v.r[1],y1=v.r[2],length=0.1,lty=3,lwd=1.5) 
  # plot projection real
  arrows(x0=0,y0=v.r[2],x1=v.r[1],y1=v.r[2],length=0.1,lty=3,lwd=1.5) 
  # arc for angle
  amag <- mag/3; cang <- 0.9*amag*c(cos(angr/2),sin(angr/2))
  arc(amag,c(0,ang)); alab<- expression(theta)
  text(cang[1],cang[2],alab,cex=0.8)

  leglab <- c(as.expression(bquote("S="*.(mag)*.(cp$units[1])*", "*theta*"="*.(ang)*degree)),
               paste0("P=",cp$P,cp$units[2],", Q=",cp$Q,cp$units[3]),
               paste0("pf=",cp$pf)
              )		  

  legend('topleft',legend=leglab,bg='white',cex=0.7)

}

complex.pow.tri <- function(cp){
 # magnitude and phase
 mag <- cp$S; ang <- cp$theta; angr <- ang*pi/180
 # max of magnitude
 xmax <- mag; xs <- seq(0,xmax,0.1)
 # input conversion to rect
 v.r <- c(cp$P,cp$Q)
 # plot axis
 plot(xs,xs,type="n",xlab="Re(S) = P[kW]",ylab="Im(S) = Q[kVAR]",
      pty="s",cex.axis=0.6)
 abline(h=0,v=0,col='gray')
 # plot arrows 
  # plot apparent
  arrows(x0=0,y0=0,x1=v.r[1],y1=v.r[2],length=0.1,lty=1,lwd=1.5) 
  text(v.r[1]/2,v.r[2]/2+0.05*xmax,"S[kVA]",cex=0.8)
  # plot real
  arrows(x0=0,y0=0,x1=v.r[1],y1=0,length=0.1,lty=2,lwd=1.5) 
  text(v.r[1]-0.05*xmax,0+0.05*xmax,"P",cex=0.8)
  # plot imaginary for triangle
  arrows(x0=v.r[1],y0=0,x1=v.r[1],y1=v.r[2],length=0.1,lty=3,lwd=1.5) 
  text(v.r[1]+0.05*xmax,v.r[2]-0.05*xmax,"Q",cex=0.8)
  amag <- mag/3; cang <- 0.9*amag*c(cos(angr/2),sin(angr/2))
  arc(amag,c(0,ang)); alab<- expression(theta)
  text(cang[1],cang[2],alab,cex=0.8)
  
  leglab <- c(as.expression(bquote("S="*.(mag)*.(cp$units[1])*", "*theta*"="*.(ang)*degree)),
               paste0("P=",cp$P,cp$units[2],", Q=",cp$Q,cp$units[3]),
               paste0("pf=",cp$pf)
              )		  

  legend('topleft',legend=leglab,bg='white',cex=0.7)
}

pf.corr <- function(P,V,pf,pfc,w=377,dig=2){
 theta <- acos(pf)*180/pi
 # P argument is in kW and V in volt
 I <- P*1000/(V*pf)
 cp <- complex.pow.calc(list(V,I,theta),dig,res=FALSE)
 # use P in kW from cp
 Sc <- P/pfc
 thetac <- acos(pfc)*180/pi
 Qc <- Sc*sin(acos(pfc))
 kWf <- tan(acos(pf)) - tan(acos(pfc))
 Qcap <- kWf*P
 Qcap <- cp$Q - Qc
 C.uF <- Qcap*1000*10^6/(w*V^2)
 Ic <- P*1000/(V*pfc)
 cpc <- complex.pow.calc(list(V,Ic,thetac),dig,res=FALSE)
 y <- round(c(V,I,Ic,kWf,Qcap,C.uF),dig)
 prnt <- array()
 prnt[1] <- paste0("P=",cp$P, cp$units[2],",  V=",y[1],"V",  ", pf=",cp$pf, ", pfc=",cpc$pf)
 prnt[2] <- paste0("S=",cp$S, cp$units[1],", theta=",cp$theta,"deg,",
            "  Q=",cp$Q, cp$units[3],",  I=",y[2],cp$units[4])
 prnt[3] <- paste0("Sc=",cpc$S, cpc$units[1],", thetac=",cpc$theta,"deg,",
            "  Qc=",cpc$Q, cpc$units[3],",  Ic=",y[3],cp$units[4])
 prnt[4] <- paste0("kWf=", y[4], ",  Qcap=",y[5],cpc$units[3], ",  Cap=",y[6],"uF")
 cat(prnt,sep="\n")

 return(list(cp=cp,cpc=cpc,y=y,prnt=prnt)) 
}

pf.corr.tri <- function(xpfc){
 # magnitude and phase
 mag <- array(); ang <- array()
 mag[1] <- xpfc$cp$S; ang[1] <- xpfc$cp$theta
 mag[2] <- xpfc$cpc$S; ang[2] <- xpfc$cpc$theta
 angr <- ang*pi/180

 # max of magnitude
 xmax <- max(mag); xs <- seq(0,xmax,0.1)
 # input conversion to rect
 v.r <- list() 
 v.r[[1]] <- c(xpfc$cp$P,xpfc$cp$Q)
 v.r[[2]] <- c(xpfc$cpc$P,xpfc$cpc$Q)

 Qlab <- c("Q","Qc"); Slab <- c("S","Sc")
 alab<- c(expression(theta),expression(theta*c))
 
 # plot axis
 plot(xs,xs,type="n",xlab="Re(S) = P[kW]",ylab="Im(S) = Q[kVAR]",
      pty="s",cex.axis=0.6)
 abline(h=0,v=0,col='gray')
 # plot arrows 
  # plot apparent
  for(i in 1:2){
   arrows(x0=0,y0=0,x1=v.r[[i]][1],y1=v.r[[i]][2],length=0.1,lty=1,lwd=1.5) 
   text(v.r[[i]][1]/2,v.r[[i]][2]/2+0.05*xmax,Slab[i],cex=0.8)
   # plot real
   arrows(x0=0,y0=0,x1=v.r[[i]][1],y1=0,length=0.1,lty=2,lwd=1.5) 
   text(v.r[[i]][1]-0.05*xmax,0+0.05*xmax,"P",cex=0.8)
   # plot imaginary for triangle
   arrows(x0=v.r[[i]][1],y0=0,x1=v.r[[i]][1],y1=v.r[[i]][2],length=0.1,lty=3,lwd=1.5) 
   text(v.r[[i]][1]+0.05*xmax,v.r[[i]][2]-0.05*xmax,Qlab[i],cex=0.8)
   amag <- mag[i]/3; cang <- 0.9*amag*c(cos(angr[i]/2),sin(angr[i]/2))
   arc(amag,c(0,ang[i]))
   text(cang[1],cang[2],alab[i],cex=0.8)
 }

   leglab <- c(paste0("P=",xpfc$cp$P, xpfc$cp$units[2],",  V=",xpfc$y[1],"V",  ", pf=",xpfc$cp$pf, ", pfc=",xpfc$cpc$pf),
               as.expression(bquote("S="*.(xpfc$cp$S)*.(xpfc$cp$units[1])*", "*theta*"="*.(xpfc$cp$theta)*degree*", Q="*
			                .(xpfc$cp$Q)*.(xpfc$cp$units[3])*",  I="*.(xpfc$y[2])*.(xpfc$cp$units[4]))),
               as.expression(bquote("Sc="*.(xpfc$cpc$S)*.(xpfc$cpc$units[1])*", "*theta*"c="*.(xpfc$cpc$theta)*degree*
                             "  Qc="*.(xpfc$cpc$Q)*.(xpfc$cpc$units[3])*",  Ic="*.(xpfc$y[3])*.(xpfc$cp$units[4]))),
               paste0("kWf=", xpfc$y[4], ",  Qcap=",xpfc$y[5],xpfc$cpc$units[3], ",  Cap=",xpfc$y[6],"uF")
              )		  

 legend('topleft',legend=leglab,bg='white',cex=0.7)
}

