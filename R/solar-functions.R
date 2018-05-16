# solar power functions

I0.orbit <- function(leap=FALSE,plot=TRUE){
 if(leap==TRUE){ny=366; npx=4} else {ny=365; npx=3} 
 np=182+npx # July 4 perihelion
 na=1+2 # January 3 aphelion

 # Solar constant
 SC = 1.377 # kW/m2
 # ET solar radiation I0
 n <- seq(1,ny)
 I0 <- SC*(1+0.033*cos((2*pi/ny)*(n-np))) # kW/m2
 ylabel <- "I0 ET radiation (kW/m2)"
 if(plot==TRUE){
  y <- I0;x <- n
  ymin <- min(y); ymax <- max(y)
  yr <- (ymax - ymin) 
  shift <- 0.02*yr
  yseq <- round(seq(ymin,ymax,0.1*yr),2)
  ypad <- c(ymin-shift, ymax+shift) 

  plot(x,y,type="l",xlab="Number of Day",ylab=ylabel,lwd=2,ylim=ypad)
  abline(h=SC,lty=2);text(20,SC+shift,"Average",cex=0.7)
  abline(v=np,lty=2);text(np+2,ymin-shift,"Perihelion",cex=0.7)
  abline(v=na,lty=2);text(na+2,ymin-shift,"Aphelion",cex=0.7)
}
return(list(I0=I0))
}

days.mo <- function(day,leap=FALSE){
 #  day for the month
 if (leap==FALSE)daymo <- c(31,28,31,30,31,30,31,31,30,31,30,31)
 else  daymo <- c(31,29,31,30,31,30,31,31,30,31,30,31)
 day.mo <-array()
 day.mo[1] <- day 
 for(i in 2:12) day.mo[i] <- day.mo[i-1]+ daymo[i-1]
 mid.mo = c(15,46,74,105,135,166,196,227,258,288,319,349)
 return(list(day.mo=day.mo))
}

declination <- function(leap=FALSE,plot=TRUE){
 # declination
 if(leap==TRUE){ny=366; ne=81} else {ny=365; ne=80} 
 n=seq(1,ny)
 dec <- 23.45*sin((2*pi/ny)*(n-ne))
if(plot==TRUE){
  plot(n,dec,type="l",xlab="Day number",ylab=expression("Declination [ "*degree*"]"),lwd=2)
  abline(h=0,lty=2);text(20,0,"Average",cex=0.7)
  abline(v=ne,lty=2);text(90,-20,"March Equinox",cex=0.7)
  abline(v=ne+183,lty=2);text(260,-20,"September Equinox",cex=0.7)
}
return(list(ne=ne,dec=dec))
}

sun.elev <- function(xdec,lat,plot=TRUE){
 elev <- xdec$dec + 90 - lat
 me <- round(mean(elev),2)
 mine <- min(elev); maxe <- max(elev)
 shift <- 0.02*(maxe - mine)
 ypad <- c(mine, 1.05*maxe) 
 n <- seq(1,length(elev),1)
 if(plot==TRUE){
  plot(n,elev,type="l",xlab="Day number",ylab=expression("Sun Angle Elev [ "*degree*"]"),lwd=2,ylim=ypad)
  abline(h=me,lty=2);text(15,me+shift,as.expression(bquote("Average= "*.(me)*degree)),cex=0.7)
  abline(v=xdec$ne,lty=2);text(xdec$ne,mine,"Equinox",cex=0.7)
  abline(v=xdec$ne+183,lty=2);text(xdec$ne+183,mine,"Equinox",cex=0.7)
  mtext(side=3,line=-1,as.expression(bquote("Lat= "*.(lat)*degree)),cex=0.7)  
}
 return(list(ne=xdec$ne,elev=elev,lat=lat))
}

read.tau <- function(file){
 loc<- read.csv(file,header=FALSE,sep=",",nrows=1)
 lat.long.elev <- read.table(file,sep=",",skip=1,header=TRUE,nrows=1)
 tau <- read.table(file,skip=3,header=TRUE,sep=",")
 return(list(loc=loc, lat.long.elev=lat.long.elev, tau=tau))
}

beam.diffuse <- function(dat,plot=TRUE){

 kab <- c(1.454,-0.406,-0.268,0.021)
 kad <- c(0.507,0.205,-0.08,-0.19)
 
 lat <- as.numeric(dat$lat.long.elev[1])
 tau <- dat$tau

 y <- sun.elev(declination(plot=FALSE),lat=lat)
 z <- I0.orbit(plot=FALSE)

 ab <- array(); ad <- array()
 Ib <- array(); Id <- array()
 day21.mo <- days.mo(21)$day.mo
 elev <- y$elev[day21.mo]
 I0 <- z$I0[day21.mo]
 air.mass0 <- 1/sin(elev*2*pi/360)
 air.mass1 <- sqrt((708*sin(elev*2*pi/360))^2+1417)-708*sin(elev*2*pi/360)
 #air.mass2 <- 1/(sin(elev*2*pi/360)+0.50572*(6.07995+elev)^-1.6364)
 #matplot(day21.mo,cbind(air.mass0,air.mass1,air.mass2),type='b') 
 air.mass <- air.mass1

 for(i in 1:12){ 
  ab[i] <- kab[1]+kab[2]*tau[i,2]+kab[3]*tau[i,3]+kab[4]*tau[i,2]*tau[i,3]
  ad[i] <- kad[1]+kad[2]*tau[i,3]+kad[3]*tau[i,2]+kad[4]*tau[i,2]*tau[i,3]
  Ib[i] <- I0[i]*exp(-tau[i,2]*air.mass[i]^ab[i]) 
  Id[i] <- I0[i]*exp(-tau[i,3]*air.mass[i]^ad[i])
}
 Id.Ib <- round(100*(Id/Ib),2)

 if(plot==TRUE){
  panels(7,7,3,2,pty="m")
  plot(day21.mo,tau[,2],xlab="Day Number",ylab="Beam Pseudo Optical Depth",type='b')
  plot(day21.mo,I0,xlab="Day Number",ylab="I0 (kW/m2)",type='b')
  plot(day21.mo,tau[,3],xlab="Day Number",ylab="Diffuse Pseudo Optical Depth",type='b')
  plot(day21.mo,Ib,xlab="Day Number",ylab="Direct Beam Ibn (kW/m2)",type='b')
  plot(day21.mo,air.mass,xlab="Day Number",ylab="Air Mass",type='b')
  plot(day21.mo,Id,xlab="Day Number",ylab="Diffuse Id (kW/m2)",type='b')
 }
 return(list(day21.mo=day21.mo,tau=tau,air.mass=air.mass,lat=lat,Ib=Ib,Id=Id,Id.Ib=Id.Ib))
}

I0.blackbody <- function(T.sun,wl.nm,plot=TRUE){

 # solar spectrum from Plancks law
 lambda <- wl.nm
 lambda.m <- lambda*10^-9 # m

 c.light <- 2.99*10^8 # m/s 
 nu <- c.light/lambda.m

 h.planck <- 6.62606957*10^-34 # Js
 k.boltz <- 1.3806488*10^-23 # J/K
 radius.earth <- 1.496*10^8; radius.sun <- 6.95*10^5
 f.rad <- (radius.sun/radius.earth)^2

 hnu.kT <- (h.planck*nu)/(k.boltz*T.sun)
 mult <- 2*pi*h.planck*c.light^2/lambda.m^5 
 I.sun.m <- f.rad* mult /(exp(hnu.kT)-1)
 I.sun.nm <- I.sun.m *10^-9 
 
 I0 <- sum(I.sun.nm)
 ymax <- max(I.sun.nm)

 if(plot==TRUE){ 
  plot(lambda,I.sun.nm,type="l",
     xlab="Wavelength (nm)",ylab="Spectral irradiance (W/m2/nm)",
     xaxs="r",yaxs="r",ylim=c(0,ymax),lwd=2,cex.axis=0.8)
  abline(v=c(400,700,1400),lty=2)
  text(c(300,550,900,1500),rep(0.05,4),c("UV","Visible","NIR","IR"),cex=0.8) 
  legend('topright', legend=c(paste("Blackbody ",T.sun,"K")),cex=0.8)
 }
 return(list(lambda=lambda,I.sun.nm=I.sun.nm,I0=I0))
}

spectral <- function(X,label,wl.lim.nm,T.sun,plot.surf=FALSE){

lambda <- X[,1] 
I0.bb.list <- I0.blackbody(T.sun=T.sun,wl.nm=lambda,plot=FALSE)
I0.bb <- I0.bb.list$I.sun.nm

I0 <- X[,2]
I <- X[,3]
ymax <- max(X[,2:4])
labelI0 <- paste(label, "I0"); labelI <- paste(label, "I");
label.bb <- paste("Blackbody ",T.sun,"K")
matplot(lambda,cbind(I0.bb,I0),type="l", col=1, lty=c(1,1),
     xlab="Wavelength (nm)",ylab="Spectral irradiance (W/m2/nm)",
     xaxs="r",yaxs="r",ylim=c(0,ymax),xlim=wl.lim.nm,lwd=c(2,0.5),cex.axis=0.8)
legend('topright', legend=c(label.bb,labelI0),
     col=1, lty=c(1,1),lwd=c(2,0.5),cex=0.8)
abline(v=c(400,700,1400),lty=2)
text(c(300,550,900,1500),rep(0.05,4),c("UV","Visible","NIR","IR"),cex=0.8) 

if(plot.surf==TRUE){
 matplot(lambda,cbind(I0.bb,I0,I),type="l", col=c('gray','gray','black'), lty=1,
     xlab="Wavelength (nm)",ylab="Spectral irradiance (W/m2/nm)",
     xaxs="r",yaxs="r",ylim=c(0,ymax),xlim=wl.lim.nm,lwd=c(2,0.5,0.5),cex.axis=0.8)
 legend('topright', legend=c(label.bb,labelI0,labelI),col=c('gray','gray','black'),
        lty=1,lwd=c(2,0.5,0.5),cex=0.8)
 abline(v=c(400,700,1400),lty=2)
 text(c(300,550,900,1500),rep(0.05,4),c("UV","Visible","NIR","IR"),cex=0.8) 
}
}

useful.waste <- function(I0.bb){
 I0 <- I0.bb
 I.sun.nm <- I0$I.sun.nm; lambda <- I0$lambda
 lambda.m <- lambda*10^-9
 c.light <- 2.99*10^8 # m/s 
 h.planck <- 6.62606957*10^-34 # Js
 k.boltz <- 1.3806488*10^-23 # J/K
 radius.earth <- 1.496*10^8; radius.sun <- 6.95*10^5
 f.rad <- (radius.sun/radius.earth)^2

 y.1100 <- I.sun.nm[which(lambda==1100)]
 x.y.1100 <- lambda[min(which(round(I.sun.nm,2)== round(y.1100,2)))]
 lines(c(1100,1100),c(0,y.1100),lty=2)
 #lines(c(x.y.1100,1100),c(y.1100,y.1100),lty=2)

 # eV/J
 eV.J <- 1.6* 10^(-19)

 # Energy of photon
 Ep <- (h.planck*c.light/lambda.m)/eV.J
 plot(lambda,Ep, type="l",xlab="Wavelength(nm)",ylab="Photon energy (eV)",lwd=2)
 abline(h=1.12,lty=2); text(400,1.3,"Eg=1.12 eV",cex=0.7)
 abline(v=1110,lty=2); text(1300,2,"1110 nm",cex=0.7)
 arrows(1110,3,600,3,length=0.1); text(800,3.2,"Enough energy",cex=0.7)
 arrows(1110,0.4,2000,0.4,length=0.1); text(1500,0.6,"Not Enough energy",cex=0.7)
 
 y <- Ep - 1.12
 for(i in 1: length(y)) if(y[i] <0) y[i] <- 0 
 waste <- sum(y)/sum(Ep)
}

sun.path <- function(lat,nday,plot=TRUE){

# declination
dec <- 23.45*sin((2*pi/365)*(nday-81))

# hr angle and trig
hr.noon <- seq(-6,+6,1) 
hr.angle <- 15*hr.noon
sin.H <- sin((pi/180)*hr.angle)
cos.H <- cos((pi/180)*hr.angle)

# elevation and angle
sin.elev <- cos((pi/180)*lat)*cos((pi/180)*dec)*cos.H + 
            sin((pi/180)*lat)*sin((pi/180)*dec)
elev <- asin(sin.elev)*180/pi

# azimuth
sin.azi <- cos((pi/180)*dec) * sin.H/cos((pi/180)*elev)

# value for testing azimuth
test.tan <- tan(dec*pi/180)/tan(lat*pi/180)

kn <- which(hr.noon==0)
nn <- length(cos.H)
azi <- array()
for(k in 1:nn){
 if(cos.H[k] >= test.tan){
   azi[k] <- asin(sin.azi[k])*180/pi
   }
   else{
   if(k <= kn) {corr = -180} else {corr = 180} 
   azi[k] <- corr  - asin(sin.azi[k])*180/pi
   }
}

if(plot==TRUE){
 wd=7; ht=7
 panels(wd,ht,2,1,pty="m")
 plot(hr.noon,elev,type="b",xlab="Hour before noon",ylab=expression("Sun angle ["*degree*"]"),lwd=2)
 text(0,10,as.expression(bquote("Lat= "*.(lat)*degree*" Day= "*.(nday))),cex.main=0.9)
 plot(hr.noon,azi,type="b",xlab="Hour before noon",ylab=expression("Azimuth ["*degree*"]"),lwd=2)
 text(-4,10,as.expression(bquote("Lat= "*.(lat)*degree*" Day= "*.(nday))),cex.main=0.9)
 plot(azi,elev,type="b",xlab=expression("Azimuth ["*degree*"]"),ylab=expression("Sun Elevation ["*degree*"]"),lwd=2)
 text(0,35,"Noon");text(-40,12,"-4h < noon"); text(40,12,"+4h > noon");
 title(as.expression(bquote("Lat="*.(lat)*degree*" Day= "*.(nday))),cex.main=0.9)
}
return(list(nday=nday, hr.noon=hr.noon, azi=azi, elev=elev))
}

sun.diagram <- function(lat){

# allocate memory for lists
# later we will put results in lists
azi.path <- list()
elev.path <- list()

# number of days in each month
daymo = c(31,28,31,30,31,30,31,31,30,31,30,31)

# select days to display in diagram
# jan 21, mar 21, june 21
display = c(0,sum(daymo[1:2]),sum(daymo[1:5])) + 21
dayname <- c("Jan 21", "Mar 21", "Jun 21")
# loop for selected days
for(id in 1:length(display)){

 # day in the set
 nday <- display[id]

 # declination
 dec <- 23.45*sin((2*pi/365)*(nday-81))

 # hr angle and trig
 hr.noon <- seq(-12,+12,1) 
 hr.angle <- 15*hr.noon
 sin.H <- sin((pi/180)*hr.angle)
 cos.H <- cos((pi/180)*hr.angle)

 # elevation and angle
 sin.elev <- cos((pi/180)*lat)*cos((pi/180)*dec)*cos.H + 
             sin((pi/180)*lat)*sin((pi/180)*dec)
 elev <- asin(sin.elev)*180/pi

 # azimuth
 sin.azi <- cos((pi/180)*dec) * sin.H/cos((pi/180)*elev)

 # value for testing azimuth
 test.tan <- tan(dec*pi/180)/tan(lat*pi/180)

 kn <- which(hr.noon==0)
 nn <- length(cos.H)
 azi <- array()
 #loop for hours 
 for(k in 1:nn){
  if(cos.H[k] >= test.tan){
   azi[k] <- asin(sin.azi[k])*180/pi
   }
   else{
   if(k <= kn) {corr = -180} else {corr = 180} 
   azi[k] <- corr  - asin(sin.azi[k])*180/pi
   }
 } # end of hour loop 

 # store results in list
 azi.path[[id]] <- azi
 elev.path[[id]] <- elev

} # end of day loop

wd=7; ht=7
panels(wd,ht,1,1,pty="m")

plot(azi.path[[1]],elev.path[[1]],type="b",
   xlab=expression("Azimuth ["*degree*"]"),ylab=expression("Elevation ["*degree*"]"),lwd=2,
   ylim=c(0,90),xlim=c(-120,120))
for(id in 2:length(display)){
 lines(azi.path[[id]],elev.path[[id]],type="b",lwd=2)
}
 text(0,30,dayname[1]) 
 text(0,50,dayname[2])
 text(0,75,dayname[3])
}


collector <- function(Ibd,sunpath,tilt,azi.c,fr,label=""){

 elev <- sunpath$elev
 azi <- sunpath$azi
 hr.noon <- sunpath$hr.noon
 nday <- sunpath$nday
 cos.inci <- cos(elev*pi/180)*cos((azi-azi.c)*pi/180)*sin(tilt*pi/180)+
             sin(elev*pi/180)*cos(tilt*pi/180)

 daymo = c(31,28,31,30,31,30,31,31,30,31,30,31)
 sum.mo <-array()
 sum.mo[1] <- daymo[1]; if (nday <= sum.mo[1]) n.month <- 1
 for(i in 2:12) {
  sum.mo[i] <- sum.mo[i-1]+ daymo[i-1]
  if (nday <= sum.mo[i] && nday > sum.mo[i-1]) n.month <- i
 }
 
 Ib <- Ibd$Ib[n.month]
 Id <- Ibd$Id[n.month]
 
 Ibc <- round(Ib*cos.inci,3)
 Idc <- round(Id*(1+cos(tilt*pi/180))/2,3)

 # reflected
 IBH <- Ib*sin(elev*pi/180)
 IDH <- Id
 Irc <- round(fr*(IBH+IDH)*(1-cos(tilt*pi/180))/2,3)
 # total

 for(i in 1:length(Ibc)) if(Ibc[i]<=0) Ibc[i] <-0 
 for(i in 1:length(Irc)) if(Irc[i]<=0) Irc[i] <-0 
 for(i in 1:length(Idc)) if(Idc[i]<=0) Idc[i] <-0 

 Ic <- Ibc+Idc+Irc

 Ibc.h <- sum(Ibc)
 Idc.h <- sum(Idc)
 Irc.h <- sum(Irc)
 Ic.h <-  sum(Ic)
 I.h <-  c(Ibc.h,Idc.h,Irc.h,Ic.h)

 ymax <- 1.2*max(Ic)

 par(mar = c(4,4,1,4))
 matplot(hr.noon,cbind(Ibc,Ic),type="l",xlab="Hour noon",
     ylab="Direct and Total (kW/m2)",ylim=c(0,ymax),lwd=1,col=1,lty=1:2,cex.axis=0.8)
 mtext(side=3,line=-1,label)
 # plot y2 right hand side axes
 par(new=T)
 matplot(hr.noon,cbind(Idc,Irc), axes=F, type="l", ylim=c(0,max(Idc+Irc)), xlab=NA,
        ylab=NA,lty=3:4,col=1,cex.lab=0.8,cex.axis=0.7,lwd=1.8)
 axis(side=4,cex.axis=0.8); mtext(side=4, line=2,"Diffuse and Reflected (kW/m2)")
 legend('topright',legend=c("Direct","Total","Diffuse","Reflected"),
       lty=1:4,col=1,lwd=1,cex=0.7)
 
 return(list(Ib=Ib,Id=Id,Ibc=Ibc,Idc=Idc,Irc=Irc,Ic=Ic,I.h=I.h))
}

tilt.adj <- function(lat,days,labels){
 tilt <- lat - round(declination(plot=F)$dec[days],2)
 barplot(tilt,names.arg=labels,ylab=expression("Tilt ["*degree*"]"),ylim=c(0,1.2*max(tilt)),
         legend=c(as.expression(bquote("Lat= "*.(lat)*degree))),args.legend=list(x='top'))
 return(list(days=days,tilt=tilt))
}

month.prod <- function(dat){

 months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
 Ibd <- beam.diffuse(dat,plot=FALSE)
 days <- days.mo(21)$day.mo
 tilt <- array();Ep <- array(); Ea <- array()
 for(i in 1:12){
 wd=7; ht=7; panels(wd,ht,2,1,pty="m")
  sunpath <- sun.path(lat=Ibd$lat,nday=days[i],plot=FALSE)
  label.day <- paste("Day",days[i], ", 21 of",months[i]) 
  Ic <- collector(Ibd, sunpath, tilt=Ibd$lat,
       azi.c=0,fr=0.2,label=paste("Polar Mount",label.day))
  Ep[i] <- Ic$I.h[4]
  tilt[i] <- 90 - max(sunpath$elev)
  Ic <- collector(Ibd, sunpath, tilt=tilt[i],
        azi.c=0,fr=0.2,label=paste("Tilt Adjusted",label.day))
  Ea[i] <- Ic$I.h[4]
 }
 E <- rbind(Ep,Ea); Emax=1.5*max(E)
 
 wd=7; ht=7; panels(wd,ht,2,1,pty="m")
 barplot(tilt,names.arg=months,ylab=expression("Tilt ["*degree*"]"),ylim=c(0,1.2*max(tilt)))
 barplot(E,ylim=c(0,Emax),names.arg=months,beside=TRUE,ylab="E (kWh/m2/d)",
         legend=c("Polar Mount","Tilt_Adjusted"))

 Eadd <- round(100*(Ea-Ep)/Ep,2)
 return(list(tilt=tilt, Ea=Ea,Ep=Ep, Eadd=Eadd))
}


two.axis.tracking <- function(dat){

 Ibd <- beam.diffuse(dat,plot=FALSE)
 days <- days.mo(21)$day.mo
 months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

 tilt <- array();Ep <- array(); Ea <- array()
 for(i in 1:12){
  wd=7; ht=7; panels(wd,ht,2,1,pty="m")
  sunpath <- sun.path(lat=Ibd$lat,nday=days[i],plot=FALSE)
  label.day <- paste("Day",days[i], ", 21 of",months[i]) 
  Ic <- collector(Ibd, sunpath, tilt=Ibd$lat,
             azi.c=0,fr=0.2,label=paste("Polar Mount",label.day))
  Ep[i] <- Ic$I.h[4]
  Ic <- collector(Ibd, sunpath, tilt=90-sunpath$elev,
             azi.c=sunpath$azi,fr=0.2,label=paste("Two-axis Tracking",label.day))
  Ea[i] <- Ic$I.h[4]
 }
 E <- rbind(Ep,Ea); Emax=1.5*max(E)
 
 months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
 wd=7; ht=3.5; panels(wd,ht,1,1,pty="m")
 barplot(E,ylim=c(0,Emax),names.arg=months,beside=TRUE,ylab="E (kWh/m2/d)",
         legend=c("Polar Mount","Two-axis Tracking"))

 Eadd <- round(100*(Ea-Ep)/Ep,2)
 mean(Eadd)
 return(list(tilt=tilt, Ea=Ea,Ep=Ep, Eadd=Eadd))
}

one.axis.tracking <- function(dat,mode='PNS'){

 Ibd <- beam.diffuse(dat,plot=FALSE)
 days <- days.mo(21)$day.mo
 months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

 tilt <- array();Ep <- array(); Ea <- array()
 for(i in 1:12){
  wd=7; ht=7; panels(wd,ht,2,1,pty="m")
  sunpath <- sun.path(lat=Ibd$lat,nday=days[i],plot=FALSE)
  label.day <- paste("Day",days[i], ", 21 of",months[i]) 
  Ic <- collector(Ibd, sunpath, tilt=Ibd$lat,
             azi.c=0,fr=0.2,label=paste("Polar Mount",label.day))
  Ep[i] <- Ic$I.h[4]
  Ic <- collector(Ibd, sunpath, tilt=90-max(sunpath$elev),
             azi.c=sunpath$azi,fr=0.2,label=paste("One-axis PNS Tracking",label.day))
  Ea[i] <- Ic$I.h[4]
 }
 E <- rbind(Ep,Ea); Emax=1.5*max(E)
 
 months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
 wd=7; ht=3.5; panels(wd,ht,1,1,pty="m")
 barplot(E,ylim=c(0,Emax),names.arg=months,beside=TRUE,ylab="E (kWh/m2/d)",
         legend=c("Polar Mount","One-axis PNS Tracking"))

 Eadd <- round(100*(Ea-Ep)/Ep,2)
 mean(Eadd)
 return(list(tilt=tilt, Ea=Ea,Ep=Ep, Eadd=Eadd))
}



