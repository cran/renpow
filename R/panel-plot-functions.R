# panel functions

# generic
panels <- function (wd,ht,rows,cols,pty,int="r"){
 np <- rows*cols # number of panels
 mat <- matrix(1:np,rows,cols,byrow=T) # matrix for layout
 layout(mat, widths=rep(wd/cols,cols), heights=rep(ht/rows,rows), TRUE)
 par(mar=c(4,4,1,.5),xaxs=int,yaxs=int,pty=pty)
}

# plot functions

# generic
plot2yaxis <- function(x,x0=FALSE,y0=FALSE){
 # argument list (x,y1,y2,xlab,y1lab,y2lab,y1symb,y2symb)
 if(x0==FALSE) xlim <- c(min(x$x),max(x$x)) else xlim <- c(0,max(x$x))
 # pad max with margin of 20%
 if(y0==FALSE){
  y1lim <-c(0.8*min(x$y1),1.2*max(x$y1))
  y2lim <-c(0.8*min(x$y2),1.2*max(x$y2))
 }else {
  y1lim <-c(0,1.2*max(x$y1))
  y2lim <-c(0,1.2*max(x$y2))
 }

 if(is.null(dim(x$y1)[2])) ny1 <- 1 else ny1 <- dim(x$y1)[2]
 if(is.null(dim(x$y2)[2])) ny2 <- 1 else ny2 <- dim(x$y2)[2]
 lty.y1 <- 1:ny1; lty.y2 <- (ny1+1):(ny1+ny2)
 
 # plot graph y1
 par(mar = c(4,4,1,4))
 matplot(x$x,x$y1,type="l", ylim=y1lim, xlim=xlim, xlab=x$xlab,
        ylab="",lty=lty.y1,col=1,cex.lab=0.8,cex.axis=0.7,lwd=1.8)
 title(ylab=x$y1lab,line=2,cex.lab=0.8)
 abline(h=0,v=0,lty=1,col='gray')


 # plot y2 right hand side axes
 par(new=T)
 matplot(x$x,x$y2, axes=F, type="l", xlim=xlim,ylim=y2lim, xlab=NA,
        ylab=NA,lty=lty.y2,col=1,cex.lab=0.8,cex.axis=0.7,lwd=1.8)
 axis(side=4,cex.axis=0.7); mtext(side=4, line=2, x$y2lab,cex=0.8)
 abline(h=0,v=0,lty=2,col='gray')

 legend('topright',legend=c(x$y1simb,x$y2simb),lty=1:(ny1+ny2),col=1,bg='white',cex=0.7,lwd=1.8)
}

