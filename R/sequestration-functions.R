# Richards

forest.seq <- function(x,y){

 A <- (x$Bmax/x$B0)^x$nu-1
 B.t <- x$Bmax/(1+A*exp(-x$r*x$nu*x$t))^(1/x$nu)
 dB.t <- x$r*B.t*(1-(B.t/x$Bmax)^x$nu)
 max.dB <- max(dB.t)
 pos.max <- which(dB.t==max.dB)
 B.max <- B.t[pos.max]
 t.max <- x$t[pos.max]


 # from closed formula
 B.max.f <- x$Bmax/(1+x$nu)^(1/x$nu)
 max.dB.f <- x$Bmax*x$r*x$nu/(1+x$nu)^((x$nu+1)/x$nu)

 CO2.t <- B.t*0.5*44/12
 dCO2.t <- dB.t*0.5*44/12
 max.dCO2 <- 0.5*(44/12)*max.dB 
 CO2.max <- 0.5*(44/12)* B.max

 tCO2.emiss <- y$P*10^-3*y$kgCO2.kWh*8760*y$C*10^-3 
 area <- tCO2.emiss/max.dCO2

 ymax <- max(CO2.t)
 par(mar = c(5,5,2,5))
 matplot(x$t,CO2.t,type="l", ylim=c(0,ymax), 
         xlab="Years",ylab="CO2 sequestered (t/ha)",lty=1,col=1)

 # plot rate right hand side axes
 par(new=T)
 ymax <- 1.4*max(dCO2.t)
 plot(x$t,dCO2.t, axes=F, type="l", ylim=c(0,ymax), xlab=NA,
        ylab=NA,lty=2,col=1)
 axis(side=4); mtext(side=4, line=3, "Seq Rate (tCO2/(yr ha))")

 legend('topright',legend=c("Seq CO2","Seq Rate"),lty=1:2,col=1,bg='white',cex=0.7)

 return(list(yr.max=t.max,B.i=B.max.f,max.dCO2=max.dCO2,
        tCO2.emiss=tCO2.emiss,area=area))

}

