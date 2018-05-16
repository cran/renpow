generate.duration <- function(file){
tp <- seq(5,20); nt <- length(tp)
Hs <- seq(1,10); nh <- length(Hs)

x <- dnorm(Hs,6,1)
y <- 8760*x/sum(x)
tpm <- seq(4,13); toth <- array()
X <- matrix(nrow=nh,ncol=nt)
for(i in nh:1){
 z <- dnorm(tp,tpm[i],3)
 X[i,] <- y[i]*z/sum(z)
 toth[i] <- sum(X[i,])
}
Xr <- round(X,0)
for(i in nh:1){
 toth[i] <- sum(Xr[i,])
}
sum(toth)
Xr[9,10] <- Xr[9,10]-2
Xr[9,11] <- Xr[9,11]-1
Xr[9,12] <- Xr[9,11]-1
Xr[9,13] <- Xr[9,13]-1

Hrs <- data.frame(Xr)
names(Hrs) <- c(tp)
rownames(Hrs) <- rev(Hs)
write.table(Hrs, file, sep=",",row.names=T,col.names=T)
}


powflux.wave <- function(Hs,Tp){
 nt <- length(Tp);nh <- length(Hs)
 Pflux <- matrix(nrow=nh,ncol=nt)
 for(i in 1:nh){
  for(j in 1:nt){
   Pflux[i,j] <- 0.42*rev(Hs)[i]^2*Tp[j]
  }
 }
 return(list(Hs=Hs,Tp=Tp,X=Pflux))
}

duration.wave <- function(datafile,file=TRUE){
if(file==TRUE){ 
 Y <- read.table(datafile,sep=",")
 Tp <- Y[1,1:length(Y[1,])][-1]
 Hs <- Y[1:length(Y[,1]),1][-1]
 X <- Y[1:length(Y[,1]),1:length(Y[1,])][-1,-1]
} else{
 Tp <- datafile[1,-1]
 Hs <- datafile[-1,1]
 X <- datafile[-1,-1]
}
return(list(Hs=c(rev(Hs)),Tp=as.numeric(Tp),X=as.matrix(X)))
}

wave.contour <- function(X,label="",sum=FALSE,sumlabel=""){
  Pr <- t(apply(X$X, 2, rev))
  contour(x=X$Tp,y=X$Hs,Pr,xlab="Tp(s)",ylab="Hs(m)",pty="s")
  mtext(side=1,line=-3,label,cex=0.7)
  if (sum==TRUE){
   mtext(side=1,line=-2,paste("Total=",round(sum(Pr),2),sumlabel),cex=0.7)
  }
}

energy.wave <- function(Pflux,D){
  X <- Pflux
  Hs <- X$Hs; Tp <- X$Tp
  Pr <- X$X
  Dr <- D$X
  Er <- round(Pr*Dr*10^-3,2)
  return(list(Hs=Hs,Tp=Tp,X=Er))
} 

energy.gen <- function(Ew,L,nu){
  X <- Ew
  Hs <- X$Hs; Tp <- X$Tp
  Eg <- L*nu*X$X
  return(list(Hs=Hs,Tp=Tp,X=Eg))
} 



 


