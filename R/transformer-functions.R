# transformer  functions

xformer.ckt <- function(x,dig=2){
 n<- x$N[2]/x$N[1]
 # impedances
 Z2.r<-x$Zo.r+x$Zl.r
 Z1.r<-Z2.r/n^2
 Zi.r<-x$Zs.r+Z1.r
 # currents	
 I1.p <- div.polar(x$Vs.p,polar(Zi.r))
 I2.p <- mult.polar(I1.p,c(1/n,0)) 
 # voltages
 V1.p <- mult.polar(I1.p,polar(Z1.r))
 V2.p <- mult.polar(V1.p,c(n,0))
 # load voltage and current
 Vl.p<- mult.polar(I2.p,x$Zl)
 Il.p<- I2.p
 # conjugate of currents
 I1.conj <- c(I1.p[1],-I1.p[2]) 
 Il.conj <- c(Il.p[1],-Il.p[2])
 # complex power 
 Ss.p <- mult.polar(x$Vs.p,I1.conj); Ss.r <- recta(Ss.p)
 Sx.p <- mult.polar(V1.p,I1.conj); Sx.r <- recta(Sx.p)
 Sl.p <- mult.polar(Vl.p,Il.conj); Sl.r <- recta(Sl.p)
 S1.p <- mult.polar(V1.p,I1.conj); S1.r <- recta(S1.p)
 
 # output
   
 ys <- list(Vs.p=x$Vs.p, Is.p=I1.p, Ss.p=Ss.p, Ss.r=Ss.r)
 yx <- list(n=n,Z2.r=Z2.r,Z1.r=Z1.r,I1.p=I1.p,V1.p=V1.p,I2.p=I2.p,
             V2.p=V2.p,S1.p=S1.p,S1.r=S1.r)
 yl <- list(Vl.p=Vl.p, Il.p=Il.p, Sl.p=Sl.p,Sl.r=Sl.r)

 lapply(ys, function(ys) ys <- round(ys,dig))
 lapply(yx, function(yx) yx <- round(yx,dig))
 lapply(yl, function(yl) yl <- round(yl,dig))

 prntys <- paste(names(ys),"=",ys,sep="",collapse=", ")
 prntyx <- paste(names(yx),"=",yx,sep="",collapse=", ")
 prntyl <- paste(names(yl),"=",yl,sep="",collapse=", ")
 prnt <- c(paste("Source: ",prntys),
           paste("Transformer: ",prntyx),
           paste("Load: ",prntyl)
          )
         
 cat(prnt, sep="\n")

 yout <- list(ys=ys,yx=yx,yl=yl,prnt=prnt)
 return(yout)
}

