read.tide <- function(file){
 x <- read.table(file, sep=",",skip=1,header=T)
 spec <- scan(file, what=character(),sep=",",nlines=1)
 site <- paste(spec[1],", ",spec[2],sep="")
 t0 <- paste(spec[4],spec[5],sep="")
 return(list(x=x,t0=t0,site=site))
}

harmonics.tide <- function(x,days,ylabel="Tide wrt MSL (m)",plot=TRUE){
 nc <- dim(x$x)[1]
 t <- seq(0,days*24,1); nt <- length(t)

 mu <- array()
 for(i in 1:nc){
  # convert to rad
  mu[i] <- x$x[,4][i]*pi/180
 }

 y <- matrix(nrow=nt,ncol=nc)
 for(i in 1:nc){
  # convert to rad
  ang.rad <- x$x[,5][i]*pi/180
  y[,i] <- x$x[,3][i]*cos(ang.rad*t+mu[i])
 }

 z <- array()
 for(it in 1:nt){
  z[it] <- sum(y[it,])
 }
 if(plot==TRUE){
  plot(t,z,type="l",xlab="Time (Hours)", ylab=ylabel)
  mtext(side=3, line=-1,paste(x$site,"t0=",x$t0),cex=0.8)
 }
 return(list(t=t,z=z,x=x$x,t0=x$t0,site=x$site))
}

find.peaks <- function(tz,band){
  x <- tz$z; t <- tz$t
  # this function finds all peaks within threshold bands
  xmin <- min(x,na.rm=T);xmax <- max(x,na.rm=T)
  thresh <- xmax*band
  # shift x up and down a position (for peak identification)
  xu <- c(tail(x, -1), NA)
  xd <- c(NA, head(x, -1))
  # identify peaks that are in the correct range 
  # where x is higher than the point before and after
  peak <- which(x - xu >= 0 & x - xd >= 0 & x >= thresh[1] & x <= thresh[2])
  xp <- x[peak]
  tp <- t[peak]
  plot(tp,xp,ylab="Peaks (m)",xlab="Time (Hours)")
  range <- 2*xp
  return(list(xp=xp,tp=tp,range=round(range,2)))  
}

power.barrage.cycle <- function(xba){
 x <- xba
 #x is list(a,A,z,nu)
 T.h.min <- c(12,24); T.dec= 12.4
 rho=1025; g=9.8
 T.sec <- T.h.min[1]*3600+ T.h.min[2]*60
 pow <- (x$a/T.sec)*rho*g*x$A*1000^2*x$z^2
 pow.tide.MW <- round(pow*10^-6,2) # MW
 pow.gen.MW <- x$nu*pow.tide.MW
 gen.MWh <- round(pow.gen.MW*T.dec,2)
 return(list(pow.tide.MW=pow.tide.MW,pow.gen.MW=pow.gen.MW,gen.MWh=gen.MWh))
}

tide.current.abs <- function(tz, ylabel="Current abs (m/s)", plot=TRUE){
 x <- tz
 z <- abs(x$z)
 if(plot==TRUE){
  plot(x$t,z,type="l",xlab="Time (Hours)", ylab=ylabel)
  mtext(side=3, line=-1,paste(x$site,"t0=",x$t0),cex=0.8)
 }
 return(list(t=x$t,z=z))
}

tidal.power <- function(tz,Aflow){
 x <- tz; A <- Aflow
 rho=1025; g=9.8
 v <- x$z
 Q <- v*A
 P <- (1/2)*rho*A*v^3
 Pavg<- mean(P)
 Pavg <- Pavg*10^-6; Punits <- "(MW)"
 if(Pavg < 10) {Pavg <- Pavg*10^3; Punits <- "(kW)"}
 return(list(Pavg=round(Pavg,2),Punits=Punits))
}

