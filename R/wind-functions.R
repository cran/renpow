# wind power functions

pow.rho.v3.table <- function(x){
 g=9.8
 P <- (1/2)*x$rho*x$A*x$v^3
 
 Pow <- P*10^-6; Punits <- "(MW)"
 if(Pow < 10) {Pow <- Pow*10^3; Punits <- "(kW)"}
 if(Pow < 10) {Pow <- Pow*10^3; Punits <- "(W)"}

 X <- t(c(x$rho,x$v,x$A,Pow)); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(3,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("Density(kg/m3)","Vel(m/s)","Area(m2)", paste("Power",Punits,sep=""))
 return(list(X=X,P=P))
}

pow.rho.v3 <- function(xw){
 x <- xw
 g=9.8
 nrho <- length(x$rho)
 nv <- length(x$v)
 Pow <- matrix(nrow=nv,ncol=nrho)
 for(i in 1:nv){
  for(j in 1:nrho){
   Pow[i,j] <- (1/2)*x$rho[j]*x$A*x$v[i]^3
  }
 }
 return(list(rho=x$rho,v=x$v,Pow=Pow))
}

pow.v3.plot <- function(x){
 ny <- length(x$y)
 nv <- length(x$v)
 matplot(x$v,x$Pow,type="l",ylab="Specific Power (W/m2)",
          xlab="Wind speed(m/s)",col=1,cex.axis=0.7,cex=0.7,cex.lab=0.7)
 legend('top',legend=paste(x$yleg,"=",x$y),lty=1:ny,cex=0.7)
 levels <- c(1,2,5)*10^0; for(i in 1:5) levels <- c(levels,c(1,2,5)*10^i)
 contour(x$v,x$y,x$Pow,levels=levels,ylab=x$ylabel,
         xlab="Wind speed(m/s)",col=1,cex.axis=0.7,cex.lab=0.7)
 mtext(side=3,line=-1,"Specific Power (W/m2)",cex=0.7)
}

rho.pT.air <- function(pT){
 x <- pT
 M = 28.97*10^-3; R = 8.314; Tref= 273
 if(x$punit=='bar') p.kPa <- x$p*100 else p.kPa <- x$p*101.325 
 p.Pa <- p.kPa*10^3
 rho <- p.Pa*M/(R*(x$T.C+Tref))

 X <- t(c(p.kPa,x$T.C,rho)); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(4,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("Pressure(kPa)","Temp(C)","Density(kg/m3)")
 return(X)
}

 rho.zT.air <- function(zT){
  x <- zT
  g=9.8
  if(x$punit=='bar') {kappa <- 348.45; beta <- kappa*g/(100*10^3)
  } else { kappa <- 353.06; beta <- kappa*g/(101.325*10^3)}
  if (is.null(x$lapse)==F){ 
   lapse.m <- x$lapse*10^-3
   T.C.z <- x$T.C - lapse.m*x$z
   T.K <- T.C.z+273; T0.K <- x$T.C+273
   p <- ((T0.K-lapse.m*x$z)/T0.K)^(beta/lapse.m)
  } else {
   T.C.z <- x$T.C
   T.K <- T.C.z+273
   p <- exp(-beta*x$z/T.K)
  }
  # density
  rho <- kappa*p/T.K

  X <- t(c(x$z,p,T.C.z,rho)); nX <- length(X)
  X <- data.frame(X)
  rdig <- rep(4,nX)
  for (i in 1:nX) X[i] <- round(X[i],rdig[i])
  names(X) <-c("Elevation (m)",paste("Pressure ","(",x$punit,")",sep=""),"Temp(C)","Density(kg/m3)")
  return(list(X=X,rho=rho))
}

pow.wind <- function(pw){
 x <- pw
 rho <- rho.zT.air(x)$rho
 xd <- list(rho=rho,v=x$v,A=1)
 X <- pow.rho.v3(xd)
 if(length(x$z)==1) y <- x$T.C else y <- x$z 
 vz <- list(v=X$v,y=y,Pow=X$Pow,yleg=x$yleg,ylabel=x$ylabel)
 pow.v3.plot(vz)
}

betz <- function(){
 x <- seq(0,1,0.01)
 nu <- 0.5*(1-x^2+x-x^3)
 plot(x,nu,type="l",xlab="v2/v1",ylab="Efficiency",ylim=c(0,0.7))
 abline(v=1/3,lty=2)
 abline(h=0.593,lty=2)
 text(1/3,0.01,"1/3")
 text(0.05,0.593+0.02,"Betz limit") 
 return()
}
 
v.H <- function(vh){
 x <- vh
 H <- seq(10,200,10); nh<-length(H)
 H0 = 10 
 nalpha <- length(x$alpha)
 v.v0.exp <- matrix(nrow=nh,ncol=nalpha)
 for(i in 1:nalpha){
 v.v0.exp[,i] <- (H/H0)^x$alpha[i]
 }
 matplot(v.v0.exp,H,type="l",col=1, cex.axis=0.7, cex.lab=0.7,xlim=c(1,3),
         xlab="Ratio v/v0",ylab="Height (m)",lwd=1.7)
 legend("bottomright", legend=paste("alpha",x$alpha),lty=1:nalpha,col=1,cex=0.7,lwd=1.7)

 nrough <- length(x$rough)
 v.v0.log <- matrix(nrow=nh,ncol=nrough)
 for(i in 1:nrough){
 v.v0.log[,i] <- log(H/x$rough[i])/log(H0/x$rough[i])
 }

 matplot(v.v0.log,H,type="l",col=1,cex.axis=0.7, cex.lab=0.7,xlim=c(1,3),
         xlab="Ratio v/v0",ylab="Height (m)",lwd=1.7)
 legend("bottomright", legend=paste("rough",x$rough),lty=1:nrough,col=1,cex=0.7,lwd=1.7)

 # power
 P.P0.exp <- v.v0.exp^3; P.P0.log <- v.v0.log^3
 matplot(P.P0.exp,H,type="l",col=1, cex.axis=0.7, cex.lab=0.7,xlim=c(1,12),
         xlab="Ratio P/P0",ylab="Height (m)",lwd=1.7)
 legend("bottomright", legend=paste("alpha",x$alpha),lty=1:nalpha,col=1,cex=0.7,lwd=1.7)
 matplot(P.P0.log,H,type="l",col=1, cex.axis=0.7, cex.lab=0.7,xlim=c(1,12),
         xlab="Ratio P/P0",ylab="Height (m)",lwd=1.7)
 legend("bottomright", legend=paste("rough",x$rough),lty=1:nrough,col=1,cex=0.7,lwd=1.7)

 return(list(v.v0.exp=v.v0.exp,v.v0.log=v.v0.log)) 
}


cal.vH <- function(calvh){
 x <- calvh
 # later look at intercept
 #cal.lm <- lm(x$v1.v2[,2]~x$v1.v2[,1])
 #a <- cal.lm$coef[1]; b<-cal.lm$coef[2]

 # for now just do without intercept
 cal.lm0 <- lm(x$v1.v2[,2]~0+ x$v1.v2[,1])
 a0 <- as.numeric(cal.lm0$coef[1]) 
 alpha <- round(log(a0)/log(x$H2/x$H1),3)
 rough <- round((x$H1^a0/x$H2)^(1/(a0-1)),3)
 plot(x$v1.v2[,1],x$v1.v2[,2],xlab="v1 (m/s)",ylab="v2 (m/s)",col='grey')
 abline(cal.lm0,lwd=2)
 mtext(side=3,line=-1,paste("alpha=",alpha,"  rough=",rough,sep=""),cex=0.8) 
 return(list(alpha=alpha, rough=rough))
}

 cdf.plot <- function(rv,xlab,ylab){
 x <- rv 
 # number of observations
 nx <- length(x)
 # data sorted	
 vx <- sort(x)
 # ranks as fractions (quantiles)
 Fx <- seq(nx)/nx
 # empirical cumulative plot	
 plot(vx, Fx, xlab=xlab, ylab=ylab,cex=0.1,
      main="Empirical Cumulative Distribution (ECDF)",cex.main=0.7)
}

weibull.plot <- function(xmax,scale,shape){
 wd=7; ht=3.5; 
 panels(wd,ht,1,2,pty="m")

 nk <- length(shape)
 quant <- seq(0,xmax,0.1); nq <- length(quant)
 dens <- matrix(nrow=nq,ncol=nk);cum <- dens
 
 for(i in 1:nk){
  dens[,i] <- dweibull(quant, shape=shape[i], scale=scale)
  cum[,i] <-  pweibull(quant, shape=shape[i], scale=scale)
 }
 scale <- round(scale,2); shape <- round(shape,2)

 matplot(quant,dens, type="l",xlab="Wind Speed(m/s)",ylab="Prob Density", 
     col=1, main=paste("Weibull pdf Scale=",scale),cex.main=0.7,cex.axis=0.7,cex.lab=0.7)
 legend('topright',legend=paste("shape=",shape),col=1,lty=1:nk,cex=0.7)
 matplot(quant,cum,type="l", xlab="Wind Speed(m/s)",ylab="Cumulative Prob",
     col=1,main = paste("Weibull cdf Scale=",scale), cex.main = 0.7,cex.axis=0.7,cex.lab=0.7)
 legend('bottomright',legend=paste("shape=",shape),col=1,lty=1:nk,cex=0.7)
}

fit.wind <-function (xd, vlabel = "Wind Speed (m/s)") { 
 x <- xd
 # calculate parameters Rayleigh pdf 
 Wv.avg <- mean(x,na.rm=T)
 # Use Eqn 13.31
 scale <- round(Wv.avg*2/sqrt(pi),2)
 shape <- round(2.0,2)
 Wv.avg <- round(Wv.avg,2)
 
 wd=7; ht=7; 
 panels(wd,ht,2,2,pty="m")

 hist(x, main = paste("Histogram"), freq=F, xlab = vlabel, cex.main = 0.7)
 #plot(density(x,na.rm=T), main = paste("Density", yrlabel), xlab = vlabel,cex.main = 0.7)
 cdf.plot(x, xlab = vlabel, ylab= 'Cumulative Prob')

 # check fit of pdf using graphs
 Wv.d <- hist(x, plot=F)$density
 quant <- hist(x, plot=F)$mids
 plot(quant, Wv.d, type="p", xlab="Wind Speed(m/s)",ylab="Prob Density", 
     main="Fit to Rayleigh pdf",cex.main=0.7)
 Wv.W <- dweibull(quant, shape=shape, scale=scale)
 lines(quant,Wv.W, type="l")
 mtext(side=1, line=-1,paste("Avg=",Wv.avg,"Shape=",shape,"Scale=",scale),cex=0.7)

 plot(ecdf(x), xlab = vlabel, lty=2,main = paste("Fit to Rayleigh cdf",""), cex.main = 0.7)
 lines(quant,pweibull(quant, shape=shape, scale=scale),lty=1)
 mtext(side=1, line=-1,paste("Avg=",Wv.avg,"Shape=",shape,"Scale=",scale),cex=0.7)
}

# re write these for classes
pow.class <- function(wc){
 x <- wc
 g=9.8; rho=1.225; v=seq(0,10);A=1
 p.int = c(0,100,150,200,250,300,400,1000)
 v.int <- c(0,4.4,5.1,5.6,6.0,6.4,7.0,9.4)
 nv <- length(v)
 Pow <- array()
 for(i in 1:nv){
   Pow[i] <- 1.91*(1/2)*rho*A*v[i]^3
 }
 plot(v,Pow,type="l",ylab="Rayleigh winds power density (W/m2)",
       xlab="Average Wind speed at 10 m (m/s)",col=1,cex.axis=0.7,
       cex=0.7,cex.lab=0.9,lwd=2)

 # classes
 for(i in 2:8){
  lines(c(v.int[i],v.int[i]),c(0,p.int[i]),lty=2)
  lines(c(0,v.int[i]),c(p.int[i],p.int[i]),lty=2)
 }
 text(v.int,0,v.int,cex=0.7)
 text(p.int,0,v.int,cex=0.7)
 mid.int <- array()
 for(i in 1:7){
  mid.int[i] <- (p.int[i]+p.int[i+1])/2
 }
 text(0.5,mid.int,paste("Class",1:7),cex=0.8)
 #Point x
 Pow.x <- round(1.91*(1/2)*rho*A*x^3,2)
 lines(c(x,x),c(0,Pow.x),col='gray',lty=3)
 lines(c(0,x),c(Pow.x,Pow.x),col='gray',lty=3)
 text(x,0,x,cex=0.7,col='gray')
 text(x,Pow.x+20,Pow.x,cex=0.7,col='gray')
 for(i in 1:7){
  if(Pow.x >p.int[i]&&Pow.x <= p.int[i+1]) class.x <- i
 }

 return(list(Pow.x=Pow.x,class.x=class.x))
}

power.curve <- function(pc){
 x <- pc
 g=9.8; rho=1.225; nv <- length(x$v)
 tune=8;shift=1.05; nu=0.4
 vrated <- 0.8*x$vrated
 P.rated <- nu*(1/2)*rho*x$A*x$vrated^3
 z.cutin <- nu*(1/2)*rho*x$A*x$cutin^3
 P.cutin <- P.rated/(1+exp(-(x$cutin-vrated/shift)/(vrated/tune)))
 Pow <- array(); z <- array()
 for(i in 1:nv){
   z[i] <- nu*(1/2)*rho*x$A*x$v[i]^3 - z.cutin
   Pow[i] <- (P.rated+P.cutin)/(1+exp(-(x$v[i]-vrated/shift)/(vrated/tune)))- P.cutin
   if(x$v[i]<= x$cutin)  {Pow[i] <- 0; z[i] <-0}
   if(x$v[i]>= x$vrated) z[i] <- P.rated
   if(x$v[i]>= x$cutout) {Pow[i] <- 0; z[i] <- 0}
 }
 Pow <- Pow*10^-3 # kW
 z <- z*10^-3
 P.rated <- P.rated*10^-3
 plot(x$v,Pow,type="l",ylab="Power (kW)",
       xlab="Wind speed (m/s)",col=1,cex.axis=0.7,
       cex=0.7,cex.lab=0.9,lwd=1)
 #lines(x$v,z,type="l")

 text(x$cutin-2,max(Pow)*0.05,paste("Cutin=",x$cutin),cex=0.7)
 text(x$vrated,0.2,paste("Rated=",x$vrated),cex=0.7)
 text(x$cutout,0.2,paste("Cutout=",x$cutout),cex=0.7)
 return(list(Pow=Pow, z=z,P.rated=P.rated))
} 

prob.power.curve <- function(pc,avg){ 
 x <- pc
 scale <- 2*avg/sqrt(pi)
 plot(x$v,pweibull(x$v,scale=scale,shape=2),type="l",cex=0.7,
     xlab="Wind Speed (m/s)", ylab="Probability",cex.lab=0.8,cex.axis=0.7)

 p.cutin <- pweibull(x$cutin,scale=scale,shape=2)
 h.cutin <- p.cutin*8760
 p.cutout <- pweibull(x$cutout,scale=scale,shape=2)
 h.cutout <- (1- p.cutout)*8760
 p.rated <- pweibull(x$vrated,scale=scale,shape=2)
 h.rated <- (1-p.rated)*8760
 h.rated.nstop <- h.rated-h.cutout
 h.run.below.rated <- (p.rated - p.cutin)*8760

x1 <- x$cutin; y1 <- p.cutin
lines(c(x1,x1),c(0,y1),lty=2)
lines(c(0,x1),c(y1,y1),lty=2)

x1 <- x$cutout; y1 <- p.cutout
lines(c(x1,x1),c(0,y1),lty=2)
lines(c(0,x1),c(y1,y1),lty=2)

x1 <- x$vrated; y1 <- p.rated
lines(c(x1,x1),c(0,y1),lty=2)
lines(c(0,x1),c(y1,y1),lty=2)

text(x$cutin-2,0,paste("Cutin=",x$cutin),cex=0.7)
text(x$vrated,0,paste("Rated=",x$vrated),cex=0.7)
text(x$cutout,0,paste("Cutout=",x$cutout),cex=0.7)
return(list(h.cutin=h.cutin,h.cutout=h.cutout,h.rated.nstop=h.rated.nstop,
            h.run.below.rated=h.run.below.rated))
}

wind.energy <- function(pc,Pow,avg){
 x <- pc
 scale <- 2*avg/sqrt(pi)
 nv <- length(x$v)
 Prob <- array()
 Prob[1]<-0
 for(i in 2:nv){
  Prob[i] <- pweibull(x$v[i],scale=scale,shape=2) - pweibull(x$v[i-1],scale=scale,shape=2)
 }
 hrs <- 8760*Prob
 energy<- round(sum(Pow$z*hrs),2)
 CF <- round(energy/(Pow$P.rated*8760),2) 
 return(list(energy=energy,unit="kWh/yr",CF=CF))
}



