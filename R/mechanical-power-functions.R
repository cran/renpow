# calculate power given a time sequence
pow.work <- function(t,pow="const",p){
 # time sequence
 nt = length(t)
 const <- function(t,c) {y<- rep(c,length(t));return (y)}
 linear <- function(t,c){y<- c*t; return (y)}
 parab <- function(t,c){y<-(1/2)*c*t^2; return (y)}
 if(pow=="const") model <- list(pow=const,wrk=linear) 
 if(pow=="linear") model <- list(pow=linear,wrk=parab) 
 # matrix p and w, array for power
 p.w <- matrix(ncol=2,nrow=nt)
 p.w[,1] <- model$pow(t,p)
 p.w[,2] <- model$wrk(t,p)
 return(list(t=t,p.w=p.w))
}

# plot mechanical power and work
pow.work.plot <- function(x){
 # pad max with margin of 20%
 ymax <- 1.2*max(x$p.w[,1]); tmax <- max(x$t)

 tlab <- expression(italic("t")*"[sec]")
 plab <- expression(italic("p")*"("*italic(t)*")[W]")
 wlab <- expression(italic("w")*"("*italic(t)*")[J]")

 # plot graph p
 par(mar = c(5,5,2,5))
 matplot(x$t,x$p.w[,1],type="l", ylim=c(0,ymax), xlab=tlab,
        ylab=plab,lty=1,col=1)

 # plot work right hand side axes
 par(new=T)
 ymax <- 1.4*max(x$p.w[,2])
 plot(x$t,x$p.w[,2], axes=F, type="l", ylim=c(0,ymax), xlab=NA,
        ylab=NA,lty=2,col=1)
 axis(side=4); mtext(side=4, line=3, wlab)

 legend('topright',legend=c("Power","Work"),lty=1:2,col=1,bg='white',cex=0.7)
}

