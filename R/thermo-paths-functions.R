# calc paths in PV and Ts planes

cpcv.cal <- function(datafile, plots=FALSE){

 X <- read.table(datafile,header=TRUE,sep=",")[-1,]
 names(X) <- c("TK","cp","cv","cpcv")
 # cp and cv in kJ/(kgK)
 M<-28.97 # multiply by M to get J/(gK)
 TK.ref <- X$TK[1]; TK <- X$TK - TK.ref
 cp.ref <- X$cp[1]*M; cp <- X$cp*M - cp.ref
 cv.ref <- X$cv[1]*M; cv <- X$cv*M - cv.ref
 cpcv.ref <- X$cpcv[1]; cpcv <- X$cpcv - cpcv.ref 

 cv.lm <- lm(cv~0+TK+I(TK^2)+I(TK^3))
 cv.lm$coef
 #         TK       I(TK^2)       I(TK^3) 
 #2.985186e-03  6.329640e-06 -3.881567e-09 
 if(plots==TRUE){
  plot(TK,cv+cv.ref,type="l")
  y <- cv.ref + cv.lm$coef[1]*TK + cv.lm$coef[2]*TK^2 + cv.lm$coef[3]*TK^3   
  lines(TK,y,type="l",lty=2)
 }

 cpcv.lm <- lm(cpcv~0+TK+I(TK^2)+I(TK^3))
 cpcv.lm$coef
 #           TK       I(TK^2)       I(TK^3) 
 #-6.587279e-05 -7.156517e-08  5.381396e-11 
 if(plots==TRUE){
  plot(TK,cpcv+X$cpcv[1],type="l")
  y <- cpcv.ref + cpcv.lm$coef[1]*TK + cpcv.lm$coef[2]*TK^2 + cpcv.lm$coef[3]*TK^3   
  lines(TK,y,type="l",lty=2)
 }
  return(list(TK.ref=TK.ref,cv.ref=cv.ref,cpcv.ref=cpcv.ref, cv.coef=cv.lm$coef,cpcv.coef=cpcv.lm$coef))
}

RefCoefAir <- list(TK.ref=300, cv.ref=20.80046, cpcv.ref=1.4,
cv.coef = c(2.992416e-03, 6.344972e-06, -3.890969e-09), 
cpcv.coef= c(-6.587279e-05, -7.156517e-08,  5.381396e-11)) 

cp.cv <- function(TC,ref=RefCoefAir){
 # TC is temp in deg C
 TK <- TC+273 - ref$TK.ref
 R <- 8.314 # universal gas constant J/(mol K)
 M <- 28.97 # g/mol
 cv.coef <- unname(ref$cv.coef)
 cpcv.coef <- unname(ref$cpcv.coef)
 cv <- ref$cv.ref + cv.coef[1]*TK^1 + cv.coef[2]*TK^2 + cv.coef[3]*TK^3   
 cpcv <- ref$cpcv.ref + cpcv.coef[1]*TK^1 + cpcv.coef[2]*TK^2 + cpcv.coef[3]*TK^3
 cp <- cpcv*cv
 cv.kg=cv/M;cp.kg=cp/M
 cv <- round(cv,3); cp <- round(cp,3); cpcv <- round(cpcv,3)
 cv.kg <- round(cv.kg,3); cp.kg <- round(cp.kg,3)
 return(list(cv=cv,cp=cp,cp.cv=cpcv,cv.kg=cv.kg,cp.kg=cp.kg))
}

fcv <- function(TC){
 y <- cp.cv(TC)$cv
 return(y)
}

fcp <- function(TC){
 y <- cp.cv(TC)$cp
 return(y)
}

fcpcv <- function(TC){
 y <- cp.cv(TC)$cp.cv
 return(y)
}

fpv <- function(V,p){
 y <- p*V
 return(y)
}

simpson <- function(fun, a, b, pts=100){
 h <- (b-a)/pts
 x <- seq(a,b,by=h)
 if(pts==2){
  s <- fun(x[1]) + 4*fun(x[2]) +fun(x[3])
 } else{
 s <- fun(x[1])+ fun(x[pts+1]) +
      2*sum(fun(x[seq(2,pts,by=2)]))+
      4*sum(fun(x[seq(3,pts-1, by=2)]))
 }
 s <- s*h/3
 return(s)
}

simpson.pv <- function(nRT, a, b, pts=100){
 h <- (b-a)/pts
 x <- seq(a,b,by=h)
 fun <- function(x) nRT/x
 if(pts==2){
  s <- fun(x[1]) + 4*fun(x[2]) +fun(x[3])
 } else{
 s <- fun(x[1])+ fun(x[pts+1]) +
      2*sum(fun(x[seq(2,pts,by=2)]))+
      4*sum(fun(x[seq(3,pts-1, by=2)]))
 }
 s <- s*h/3
 return(s)
}


path.calc <- function(x){
 if(x$path=='isotherm') y <- isotherm(x)
 if(x$path=='adiabat')  y <- adiabat(x)
 if(x$path=='isochor')  y <- isochor(x)
 if(x$path=='isobar')   y <- isobar(x)
 return(y) 
}

path.summary <- function(y){
 # return from any path function is 
 # list(v=v,V=V,P=P,T=T,s=s,S=S,W=W,Q=Q,cv=cv,cp=cp,gamma=gamma,WQtot=WQtot,call=call,nM=nM)
 nvars <- 11
 pts <- length(y$v)
 start.end <- data.frame(matrix(nrow=2,ncol=nvars))
 rdig <- c(3,1,1,1,3,1,1,1,3,3,3)
 for (i in 1:nvars) start.end[,i] <- round(c(y[[i]][1],y[[i]][pts]),rdig[i])
 names(start.end) <-c("v(m3/kg)","V(l)","P(b)","T(C)","s(J/Kkg)",
                      "S(J/K)","W(kJ/kg)","Q(kJ/kg)","cv(kJ/Kkg)","cp(kJ/Kkg)","cp/cv")
 WQtot <- data.frame(t(y$WQtot)); names(WQtot) <- c("Wtot(kJ/kg)","Qtot(kJ/kg)")
 nM <- data.frame(t(y$nM)); names(nM) <- c("n(mol)","M(g/mol)", "m(kg)")
 return(list(start.end=start.end,WQtot=WQtot,pts=pts,call=y$call,nM=nM))
}
 
# ------------- detailed calc -------------

isochor <- function(x){
 mode <- 'none'
 res <- 0.001
 # arguments: V(l), P(bar), n(mol), T(C), M(g/mol)
 if(is.null(x$V)) {mode='V';Vab=c(NA,NA)} else Vab <- x$V
 if(is.null(x$P)) {mode='P';Pab=c(NA,NA)} else Pab <- x$P
 if(is.null(x$T)) {mode='T';Tab=c(NA,NA)} else Tab <- x$T
 if(is.null(x$S)) Sab=c(1,NA) else Sab <- x$S
 if(is.null(x$n))   n<-1 else n<- x$n
 if(is.null(x$M))   M<-28.97 else M <- x$M

 R <- 8.314 # universal gas constant J/(mol K)
 liter.m3 <- 0.001; bar.pa <- 10^5
 Vab <- Vab*liter.m3; Pab <- Pab*bar.pa; Tab <- Tab +273
 if(mode=='V') Vab[1] <- n*R*Tab[1]/Pab[1]
 if(mode=='P') Pab <- n*R*Tab/Vab[1]
 if(mode=='T') Tab <- Pab*Vab[1]/(n*R)

 dP <- (Pab[2]-Pab[1])*res
 P <- seq(Pab[1],Pab[2],dP); np <- length(P)
 V <- rep(Vab[1],np);dV=diff(V)
 T <- P*V/(n*R); TC <- T-273
 W <- array();Q <- array(); S <- array() 
 W[1] <- 0; Q[1] <-0; S[1] <- Sab[1]
 for(i in 2:np){
  dW <- -P[i]*dV[i-1] # increasing V makes W -
  dQ <- n*simpson(fcv,TC[i-1],TC[i],pts=2) #dT + makes Q +
  dS <- dQ/T[i]
  W[i] <- W[i-1] + dW
  Q[i] <- Q[i-1] + dQ
  S[i] <- S[i-1] + dS 
 }

 # specific volume; mass kg from mol and molar mass
 v <- round(V/(n*M*0.001),5)	
 V <- round((V/liter.m3),5)
 P <- round((P/bar.pa),5)
 T <- round(TC,5)
 S <- round(S,5)
 # specific entropy
 s <- round(S/(n*M),5)	#(kJ/K)/kg
 W <- round(W/(n*M),5)  #kJ/kg
 Q <- round(Q/(n*M),5)  #kJ/kg
 cv <- round(cp.cv(TC)$cv.kg,3)
 cp <- round(cp.cv(TC)$cp.kg,3)
 gamma <- round(cp.cv(TC)$cp.cv,3)
 WQtot <- round(c(W[length(W)],Q[length(Q)]),1)  #kJ
 call <- noquote(unlist(x))
 nM <- c(n,M,round(0.001*n*M,5)) # n mol, M g/mol, m kg
 return(list(v=v,V=V,P=P,T=T,s=s,S=S,W=W,Q=Q,cv=cv,cp=cp,gamma=gamma,WQtot=WQtot,call=call,nM=nM))
}

isobar <- function(x){
 mode <- 'none'
 res <- 0.001
 # arguments: V(l), P(bar), n(mol), T(C), M(g/mol)
 if(is.null(x$V)) {mode='V';Vab=c(NA,NA)} else Vab <- x$V
 if(is.null(x$P)) {mode='P';Pab=c(NA,NA)} else Pab <- x$P
 if(is.null(x$T)) {mode='T';Tab=c(NA,NA)} else Tab <- x$T
 if(is.null(x$S)) Sab=c(1,NA) else Sab <- x$S
 if(is.null(x$n))   n<-1 else n<- x$n
 if(is.null(x$M))   M<-28.97 else M <- x$M
 R <- 8.314 # universal gas constant J/(mol K)

 liter.m3 <- 0.001; bar.pa <- 10^5
 Vab <- Vab*liter.m3; Pab <- Pab*bar.pa; Tab <- Tab +273
 if(mode=='V') Vab <- n*R*Tab/Pab[1]
 if(mode=='P') Pab[1] <- n*R*Tab/Vab[1]
 if(mode=='T') Tab <- Pab[1]*Vab/(n*R)

 dV <- (Vab[2]-Vab[1])*res
 V <- seq(Vab[1],Vab[2],dV); nV <- length(V)
 P <- rep(Pab[1],nV)
 T <- P*V/(n*R); TC <- T-273
 W=array();Q <- array(); S <- array() 
 W[1] <- 0; Q[1] <-0; S[1] <- Sab[1]

 for(i in 2:nV){
  dW <- - P[i]*dV # increasing V makes W -
  dQ <- n*simpson(fcp,TC[i-1],TC[i],pts=2) #dT + makes Q +
  dS <- dQ/T[i]
  W[i] <- W[i-1] + dW
  Q[i] <- Q[i-1] + dQ
  S[i] <- S[i-1] + dS 
 }
 # specific volume; mass kg from mol and molar mass
 v <- round(V/(n*M*0.001),5)	
 V <- round((V/liter.m3),5)
 P <- round((P/bar.pa),5)
 T <- round(TC,5)
 S <- round(S,5)
 # specific entropy
 s <- round(S/(n*M),5)	#(kJ/K)/kg
 W <- round(W/(n*M),5)  #kJ/kg
 Q <- round(Q/(n*M),5)  #kJ/kg
 cv <- round(cp.cv(TC)$cv.kg,3)
 cp <- round(cp.cv(TC)$cp.kg,3)
 gamma <- round(cp.cv(TC)$cp.cv,3)
 WQtot <- round(c(W[length(W)],Q[length(Q)]),1)  #kJ
 call <- noquote(unlist(x))
 nM <- c(n,M,round(0.001*n*M,5)) # n mol, M g/mol, m kg
 return(list(v=v,V=V,P=P,T=T,s=s,S=S,W=W,Q=Q,cv=cv,cp=cp,gamma=gamma,WQtot=WQtot,call=call,nM=nM))
}

isotherm <- function(x){
 mode <- 'none'
 res <- 0.001

 # arguments: V(l), P(bar), n(mol), T(C), M(g/mol)
 if(is.null(x$V)) {mode='V';Vab=c(NA,NA)} else Vab <- x$V
 if(is.null(x$P)) {mode='P';Pab=c(NA,NA)} else Pab <- x$P
 if(is.null(x$T)) {mode='T';Tab=c(NA,NA)} else Tab <- x$T
 if(is.null(x$S)) Sab=c(1,NA) else Sab <- x$S
 if(is.null(x$n))   n<-1 else n<- x$n
 if(is.null(x$M))   M<-28.97 else M <- x$M

 R <- 8.314 # universal gas constant J/(mol K)
 # convert to m3 and pascal
 liter.m3 <- 0.001; bar.pa <- 10^5
 Vab <- Vab*liter.m3; Pab <- Pab*bar.pa; Tab <- x$T+273
 # decide calc mode
 if(mode=='P') Pab <- n*R*Tab[1]/Vab
 if(mode=='V') Vab <- n*R*Tab[1]/Pab
 if(mode=='T') Tab[1] <- Pab[1]*Vab[1]/(n*R)

 dV <- (Vab[2]-Vab[1])*res
 V <- seq(Vab[1],Vab[2],dV); nv <- length(V)
 T <- rep(Tab[1],nv); TC <- T-273
 P <- n*R*T/V
 W <- array(); Q <- array(); S <- array() 
 W[1]<- 0; Q[1] <-0; S[1] <- Sab[1]

 for(i in 2:nv){
  #dW <- - P[i]*dV # increasing V makes W -
  dW <- - simpson.pv(nRT=n*R*T[i],V[i-1],V[i],pts=2)
  dQ <- - dW # opposite sign of W
  dS <- dQ/T[i]
  W[i] <- W[i-1] + dW
  Q[i] <- Q[i-1] + dQ
  S[i] <- S[i-1] + dS 
 }

 # specific volume; mass kg from mol and molar mass
 v <- round(V/(n*M*0.001),5)	
 V <- round((V/liter.m3),5)
 P <- round((P/bar.pa),5)
 T <- round(TC,5)
 S <- round(S,5)
 # specific entropy
 s <- round(S/(n*M),5)	#(kJ/K)/kg
 W <- round(W/(n*M),5)  #kJ/kg
 Q <- round(Q/(n*M),5)  #kJ/kg
 cv <- round(cp.cv(TC)$cv.kg,3)
 cp <- round(cp.cv(TC)$cp.kg,3)
 gamma <- round(cp.cv(TC)$cp.cv,3)
 WQtot <- round(c(W[length(W)],Q[length(Q)]),1)  #kJ
 call <- noquote(unlist(x))
 nM <- c(n,M,round(0.001*n*M,5)) # n mol, M g/mol, m kg
 return(list(v=v,V=V,P=P,T=T,s=s,S=S,W=W,Q=Q,cv=cv,cp=cp,gamma=gamma,WQtot=WQtot,call=call,nM=nM))
}  

adiabat <- function(x){
 mode <- 'none'
 res <- 0.001
 # arguments: V(l), P(bar), n(mol), T(C), M(g/mol)
 if(is.null(x$V)) {mode='V';Vab=c(NA,NA)} else Vab <- x$V
 if(is.null(x$P)) {mode='P';Pab=c(NA,NA)} else Pab <- x$P
 if(is.null(x$T)) {mode='T';Tab=c(NA,NA)} else Tab <- x$T
 if(is.null(x$S)) Sab=c(1,NA) else Sab <- x$S
 if(is.null(x$n))   n<-1 else n<- x$n
 if(is.null(x$M))   M<-28.97 else M <- x$M

 R <- 8.314 # universal gas constant J/(mol K)
 # specific heats and ratio gamma
 # convert to m3 and pascal
 liter.m3 <- 0.001; bar.pa <- 10^5
 Vab <- Vab*liter.m3; Pab <- Pab*bar.pa;Tab <- x$T+273
 # decide calc mode
 
 if(mode=='P')	Pab[1] <- n*R*Tab[1]/Vab[1]
 if(mode=='V')	Vab[1] <- n*R*Tab[1]/Pab[1]
 if(mode=='T') Tab[1] <- Pab[1]*Vab[1]/(n*R)

  if(!anyNA(Vab)){
  dV <- (Vab[2]-Vab[1])*res
  V <- seq(Vab[1],Vab[2],dV); nV <-length(V)
  P <- V; T <- V
  P[1] <- Pab[1];T[1] <- Tab[1]
  for(i in 2:nV){
   gamma.i <- cp.cv(TC=T[i-1]-273)$cp.cv
   T[i] <- T[i-1]*(V[i-1]/V[i])^(gamma.i-1)
   T.i <- (T[i]+T[i-1])/2 
   gamma.i <-  cp.cv(TC=T.i - 273)$cp.cv 
   T[i] <- T[i-1]*(V[i-1]/V[i])^(gamma.i-1)
   P[i] <- P[i-1]*(V[i-1]/V[i])^gamma.i
  }
  npts <- nV
 }

  if(!anyNA(Pab)){
  dP <- (Pab[2]-Pab[1])*res
  P <- seq(Pab[1],Pab[2],dP); nP <-length(P)
  V <- P; T <- P
  V[1] <- Vab[1];T[1] <- Tab[1]
  for(i in 2:nP){
   gamma.i <- cp.cv(TC=T[i-1]-273)$cp.cv
   V[i] <- V[i-1]*(P[i]/P[i-1])^(-1/gamma.i)
   T[i] <- T[i-1]*(V[i-1]/V[i])^(gamma.i-1)
   T.i <- (T[i]+T[i-1])/2 
   gamma.i <-  cp.cv(TC=T.i - 273)$cp.cv 
   V[i] <- V[i-1]*(P[i]/P[i-1])^(-1/gamma.i)
   T[i] <- T[i-1]*(V[i-1]/V[i])^(gamma.i-1)
  }
  npts <- nP
 }

  if(!anyNA(Tab)){
  dT <- (Tab[2]-Tab[1])*res
  T <- seq(Tab[1],Tab[2],dT); nT <-length(T)
  P <- T; V <- T
  P[1] <- Pab[1];V[1] <- Vab[1]
  for(i in 2:nT){
   T.i <- (T[i]+T[i-1])/2 
   gamma.i <-  cp.cv(TC=T.i-273)$cp.cv 
   V[i] <- V[i-1]*(T[i]/T[i-1])^(1/(1-gamma.i))
   P[i] <- P[i-1]*(V[i-1]/V[i])^gamma.i
  }
  npts <- nT
 }

 W <- array(); Q <- array(); S <- array() 
 W[1]<- 0; Q[1] <-0; S[1] <- Sab[1]
 for(i in 2:npts){
  #dW <- - P[i]*dV[i-1] # increasing V makes W -
  dW <- - simpson.pv(P[i]*V[i],V[i-1],V[i],pts=2)
  dQ <- 0 # # zero heat adiabatic
  dS <- dQ/T[i]
  W[i] <- W[i-1] + dW
  Q[i] <- Q[i-1] + dQ
  S[i] <- S[i-1] + dS 
 }

 # specific volume; mass kg from mol and molar mass
 v <- round(V/(n*M*0.001),5)	
 V <- round((V/liter.m3),5)
 P <- round((P/bar.pa),5)
 T <- round(T-273,5)
 S <- round(S,5)
 # specific entropy
 s <- round(S/(n*M),5)	#(kJ/K)/kg
 W <- round(W/(n*M),5)  #kJ/kg
 Q <- round(Q/(n*M),5)  #kJ/kg
 cv <- round(cp.cv(TC=T)$cv.kg,3)
 cp <- round(cp.cv(TC=T)$cp.kg,3)
 gamma <- round(cp.cv(TC=T)$cp.cv,3)
 WQtot <- round(c(W[length(W)],Q[length(Q)]),1)  #kJ
 call <- noquote(unlist(x))
 nM <- c(n,M,round(0.001*n*M,5)) # n mol, M g/mol, m kg
 return(list(v=v,V=V,P=P,T=T,s=s,S=S,W=W,Q=Q,cv=cv,cp=cp,gamma=gamma,WQtot=WQtot,call=call,nM=nM))
}  


# TS plane phase function
phase <- function(){
 s <- seq(0,10,0.1)
 scrit <- 4.5; svar <- 1.8
 y <- 1875*((1/sqrt(2*pi*svar^2))*exp(-(s-scrit)^2/(2*svar^2))) 
 ymax <- max(y)
 plot(s,y,type="l",lty=2,xlab="Entropy s [kJ/kgK]", ylab=expression("Temperature T ["*degree*"C]"),
	  ylim=c(0,1.4*ymax))
 plot(s,y,type="l",lty=2,xlab="Entropy s [kJ/kgK]", ylab=expression("Temperature T ["*degree*"C]"),
	  xaxt="n",yaxt="n",ylim=c(0,1.4*ymax))
 return(list(s=s,y=y))
}
 

