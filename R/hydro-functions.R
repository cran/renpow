P.hA <- function(x){
 rho=1000; g=9.8
 v <- sqrt(2*g*x$h)
 Q <- v*x$A
 P <- (1/2)*rho*x$A*v^3*10^-6; Punits <- "(MW)"
 if(P < 10) {P <- P*10^3; Punits <- "(kW)"}
 X <- t(c(x$h,v,x$A,Q,P)); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(2,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("Head(m)","Vel(m/s)","Area(m2)", "Flow(m3/s)",paste("Power",Punits,sep=""))
 return(X)
}

P.Qh <- function(x){
 rho=1000; g=9.8
 Pw <- rho*g*x$Q*x$h*10^-6; Punits <- "(MW)"
 if(Pw < 10) {Pw <- Pw*10^3; Punits <- "(kW)"}
 X <- t(c(x$h,x$Q,Pw)); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(2,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("Head(m)","Flow(m3/s)",paste("Power",Punits,sep=""))
 return(X)
}
 
Pmax.Qh <- function(x){
  rho=1000; g=9.8
  net.h <- x$h - x$h/3
  Pw <- rho*g*x$Q*net.h*10^-6; Punits <- "(MW)"
  if(Pw < 10) {Pw <- Pw*10^3; Punits <- "(kW)"}
  X <- t(c(x$h,net.h,x$Q,Pw)); nX <- length(X)
  X <- data.frame(X)
  rdig <- rep(2,nX)
  for (i in 1:nX) X[i] <- round(X[i],rdig[i])
  names(X) <-c("Gross head (m)","Net head (m)","Flow (m3/s)",paste("Power",Punits,sep=""))
  return(X)
}  

Pe.Pw <- function(x){
 rho=1000; g=9.8
 net.h <- x$h - x$h/3
 Pw <- rho*g*x$Q*net.h*10^-6; Punits <- "(MW)"
 if(Pw < 10) {Pw <- Pw*10^3; Punits <- "(kW)"}
 Pe <- x$nu*Pw
 press <- rho*g*net.h/1000
 X <- t(c(x$h,net.h,x$Q,press,x$nu,Pw,Pe)); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(2,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("GrossHead(m)","NetHead(m)","Flow(m3/s)","Press(kPa)","Eff", paste("PowWater",Punits,sep=""),paste("PowGen",Punits,sep=""))
 return(X)
}

Pmax.Qh.plot <- function(x){
 rho=1000; g=9.8
 Q <- c(0.1,0.5,1,5,10,50,100,500,1000,2000);nQ <- length(Q)
 P <- c(0.01,0.1,1,10,100,1000); nP <- length(P)
 h <- matrix(nrow=nQ,ncol=nP)
 for(i in 1:nP){
  h[,i] <- (P[i]*10^6)/(rho*g*Q)
 }
 par(pty='s')
 matplot(Q,h,type="l", col=1, lty=1, log='xy',ylim=c(1,2000),xlim=c(0.1,2000),
        xaxs='i',yaxs='i',xaxt='n',yaxt='n',ylab="Net Head (m)",xlab="Flow (m3/s)",
        lwd=1.5,cex.axis=0.8)
 tck <- c(0.1,0.2,0.5); xtck <- c(tck,10*tck,100*tck,1000*tck,10000*tck)
 ytck <- c(10*tck,100*tck,1000*tck,10000*tck)

 axis(1,at=xtck,labels=format(xtck,scientific=FALSE),cex.axis=0.6)
 axis(2,at=ytck,labels=format(ytck,scientific=FALSE),cex.axis=0.6)
 points(x$Q,x$h-x$h/3,type="p", pch=x$plab,cex=1,font=2)

 slab <- c(0.5,1.6,5,16,50,600)
 for(i in 1:nP){
  xt <- slab[i]; yt <- slab[i]*8
  text(xt,yt,paste(P[i],"MW",sep=""),srt=-55,cex=0.8)
 }
 grid(col='gray')
 
}

turbine.regions <- function(type){
 # kaplan
  if(type=='kaplan'){
  polygon(x=c(1,1,10,200,1000,80,1),y=c(1,20,80,80,15,1,1),lwd=2,border='gray')
  text(50,70,"Kaplan")
 }
 # francis
 if(type=='francis'){
  polygon(x=c(1,5,100,1000,1000,5,1),y=c(80,800,800,80,10,10,80),lwd=2,border='gray')
  text(200,200,"Francis",srt=-45)
 }
 # pelton
 if(type=='pelton'){
  polygon(x=c(1,1,20,50,40,1),y=c(80,1000,1000,800,600,80),lwd=2,border='gray')
  text(2,900,"Pelton")
 }
 # crossflow
 if(type=='crossflow'){ 
  polygon(x=c(0.5,0.8,20,15),y=c(3,5,5,3),lwd=2,border='gray')
  text(5,4,"Crossflow")
 }
 # SLH
 if(type=='slh'){ 
  polygon(x=c(0.5,0.8,20,15),y=c(3,5,5,3),lwd=2,border='gray')
  text(5,4,"SLH")
 }
}

turbine.regions.all <- function(){
 # kaplan
  polygon(x=c(1,1,10,200,1000,80,1),y=c(1,20,80,80,15,1,1),lwd=2,border='gray',lty=1)
  text(50,70,"Kaplan")
 # francis
  polygon(x=c(1,5,100,1000,1000,5,1),y=c(80,800,800,80,10,10,80),lwd=2,border='gray',lty=2)
  text(200,200,"Francis",srt=-45)
 # pelton
  polygon(x=c(1,1,20,50,40,1),y=c(80,1000,1000,800,600,80),lwd=2,border='gray',lty=3)
  text(2,900,"Pelton")
 # SLH
  polygon(x=c(0.5,0.8,20,15),y=c(3,5,5,3),lwd=2,border='gray',lty=4)
  text(5,4,"SLH")
}

pipe.loss <- function(pipe){
 x <- pipe
 k1=10.67; k2=1.852; k3=4.87 
 if (x$mat=='pvc') C <-c(150,150)
 if (x$mat=='concrete') C <-c(100,140)
 if (x$mat=='steel') C <-c(90,110)
 if (x$mat=='galvanized') C <-c(120,120)
 if (x$mat=='poly') C <-c(140,140)
 h.loss <- k1*(x$L/x$d^k3)*(x$Q/mean(C))^k2 
 X <- t(c(h.loss,mean(C))); nX <- length(X)
 X <- data.frame(X)
 rdig <- rep(2,nX)
 for (i in 1:nX) X[i] <- round(X[i],rdig[i])
 names(X) <-c("Head loss(m)","Roughness")
 return(X)
}

#   - ----- hydrology -----------
# exceedance prob
exceed <- function(flow){
 y <- flow; z <- sort(y)
 proby <- sort(rank(y),decreasing=T)/length(y)
 quant <- c(0.5,0.95)
 yp <- list(); yp.ea <- array()
 for(i in 1:length(quant)){
  yp[[i]] <- z[which(proby<=quant[i])]
  yp.ea[i] <- yp[[i]][1]
 }
 ym <- mean(y)
 pym <- proby[which(z>=ym)][1]

 prob <- round(quant,2)
 Q <- round(yp.ea,2)
 Prob.Qmean <- round(c(pym,ym),2)
 X <- rbind(c(quant,pym),c(yp.ea,ym))
 X <- data.frame(X)
 X <- round(X,2)
 names(X) <-c(c("Q50","Q95","Qmean"))
 row.names(X) <- c("Prob","Q")
 return(list(y=z,proby=proby, prob=prob, Q=Q, Prob.Qmean=Prob.Qmean, prob.Q=X))
}

# hypothetical stream flow
model.flow <- function(mf){
  x <- mf
  days <- seq(1,365)
  seasonal <- exp(-((days-x$day.peak)/x$length.season)^2) 
  delta <- x$peak.flow - x$base.flow
  noise <- rnorm(365,delta*x$variab[1],x$variab[2])
  flow <- array()
  flow[1:3] <- seasonal[1:3]   
  for (i in 4:365){
   flow[i] <-0
   for (j in 1:3) flow[i] <- flow[i]+seasonal[i]*x$coef[j]*flow[i-j] 
   flow[i] <- flow[i]+seasonal[i]*x$coef[4]*noise[i]
  }
  flow <- flow*x$peak.flow + x$base.flow
  for (i in 1:365)if(flow[i]<=0) flow[i] <-0
 return(flow=flow)
}

flow.plot <- function(flow,label){
 Qmean <- mean(flow)
 days <- seq(1,365)
 plot(days,flow,type="l",lty=1, ylim=c(0,max(flow)),
      ylab=label, xlab="Days", cex.lab=0.8)
 abline(h=Qmean,lty=2,lwd=0.7)
 text(20,Qmean,"Qmean",cex=0.7)
}

annual.avg <- function(mf,nyrs){
 x <- mf
 Xt <- matrix(nrow=365,ncol=nyrs)
 for(i in 1:nyrs)Xt[,i] <- model.flow(x)
 Xtm <- array()
 for (j in 1:365) Xtm[j] <- mean(Xt[j,])
 return(Xtm)
}

# hypothetical stream flow
flow.exc.plot <- function(flow,exc,label){
 levels <- list(lev=exc$Q,lev.lab=names(exc$prob.Q)[-length(names(exc$prob.Q))])
 Qmean <- exc$Prob.Qmean[2]; prob.Qm <-exc$Prob.Qmean[1]
 days <- seq(1,365)
 plot(days,flow,type="l",lty=1, ylim=c(0,max(flow)),
      ylab=label, xlab="Days", cex.lab=0.8)
 for(k in 1:length(levels$lev)){
  abline(h=levels$lev[k],lty=2,lwd=0.7)
  text(2,levels$lev[k],levels$lev.lab[k],cex=0.7)
 }
  abline(h=Qmean,lty=2,lwd=0.7)
  text(20,Qmean,"Qmean",cex=0.7)

 plot(exc$proby,exc$y,type="l",ylim=c(0,max(exc$y)),xlim=c(0,1),lty=1,
      ylab=label,xlab="Exceedance probability", cex.lab=0.8)
 for(k in 1:length(levels$lev)){
  abline(h=levels$lev[k],lty=2,lwd=0.7)
  text(0.01,levels$lev[k],levels$lev.lab[k],cex=0.7)
 }
  abline(h=Qmean,lty=2,lwd=0.7)
  text(0.1,Qmean,"Qmean",cex=0.7)
  abline(v=prob.Qm,lty=2,lwd=0.7)
  text(prob.Qm,0,paste("ProbQmean=",prob.Qm),cex=0.7)
}

area.vol <- function(xav){
 x <- xav
 # simplistic bathymetry V cross sectional profile and axial ramp 
 h <- seq(0,1,0.01)*(x$H-x$B)
 tail <- (x$L*1000/x$H)*h
 area <- (x$W*1000/x$H)*h*tail*10^-4
 vol <- (area/2)*h
 vmax <- max(vol)
 amax <- max(area)

 par(mar = c(5,5,5,2))
 plot(vol,h+x$B,type="l",xlim=c(0,vmax),lty=1,col=1, 
        xlab="Volume (ha.m)", ylab="Lake elevation asl (m)")
 abline(h=x$H,col='gray')
 text(vmax/2,x$H,"Cons. pool elevation (m)",cex=0.8)
 par(new=T)
 plot((area),h,type="l",xlim=c(amax,0),axes=F,xlab=NA,
        ylab=NA,lty=2,col=1)
 axis(side=3); mtext(side=3, line=3, "Area (ha)")
 legend('bottom',legend=c("Volume","Area"),lty=1:2,col=1,bg='white',cex=0.7)
}

