# path.cycles draw cycles in PV and Ts planes with calc and units
 
path.cycles <- function(x,plane='Pv',shade.cycle=FALSE){

if(x$cty=='carnot'){
 TH=x$TH;TL=x$TL;V1=x$V1;V4=x$V4
 if(is.null(x$S1)) S1=1 else S1 <- x$S1
  x1 <- list(V=c(V1,NA),T=c(TL,TH),S=c(S1,NA),path='adiabat',lab=c("1","2")) 
  y1 <- path.calc(x1); V2 <- y1$V[length(y1$V)]; V3=V2*V4/V1
  x2 <- list(V=c(V2,V3),T=c(TH,NA),S=c(S1,NA),path='isotherm',lab=c("2","3"))
  y2 <- path.calc(x2); S2 <- y2$S[length(y2$S)]
  x3 <- list(V=c(V3,NA),T=c(TH,TL),S=c(S2,NA),path='adiabat',lab=c("3","4"))
  y3 <- path.calc(x3)
  x4 <- list(V=c(V4,V1),T=c(TL,NA),S=c(S2,NA),path='isotherm',lab=c("4","1"))
  y4 <- path.calc(x4)
  z <- list(x1,x2,x3,x4)
  path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
  y <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
  return(list(y=y)) 
}

if (x$cty=='brayton'){
 TH=x$TH;TL=x$TL;P1=x$P1;P2=x$P2
 if(is.null(x$S1)) S1=1 else S1 <- x$S1
  x1 <- list(P=c(P1,P2),T=c(TL,NA),S=c(S1,NA),path='adiabat',lab=c("1","2")) 
  y1 <- path.calc(x1); T2 <- y1$T[length(y1$T)]; S2<-S1 
  x2 <- list(P=c(P2,P2),T=c(T2,TH),S=c(S2,NA),path='isobar',lab=c("2","3"))
  y2 <- path.calc(x2); S3 <- y2$S[length(y2$S)]
  x3 <- list(P=c(P2,P1),T=c(TH,NA),S=c(S3,NA),path='adiabat',lab=c("3","4"))
  y3 <- path.calc(x3); T4 <- y3$T[length(y3$T)];S4 <- S3
  x4 <- list(P=c(P1,P1),T=c(T4,TL),S=c(S4,NA),path='isobar',lab=c("4","1"))
  y4 <- path.calc(x4)
  z <- list(x1,x2,x3,x4)
  path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
  y <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
  return(list(y=y)) 
}

if (x$cty=='otto'){
 TH=x$TH;TL=x$TL;V1=x$V1;V2=x$V2
 if(is.null(x$S1)) S1=1 else S1 <- x$S1
  x1 <- list(V=c(V1,V2),T=c(TL,NA),S=c(S1,NA),path='adiabat',lab=c("1","2")) 
  y1 <- path.calc(x1); T2 <- y1$T[length(y1$T)]; S2<-S1 
  x2 <- list(V=c(V2,V2),T=c(T2,TH),S=c(S2,NA),path='isochor',lab=c("2","3"))
  y2 <- path.calc(x2); S3 <- y2$S[length(y2$S)]
  x3 <- list(V=c(V2,V1),T=c(TH,NA),S=c(S3,NA),path='adiabat',lab=c("3","4"))
  y3 <- path.calc(x3); T4 <- y3$T[length(y3$T)];S4 <- S3
  x4 <- list(V=c(V1,V1),T=c(T4,TL),S=c(S4,NA),path='isochor',lab=c("4","1"))
  y4 <- path.calc(x4)
  z <- list(x1,x2,x3,x4)
  path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
  y <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
  return(list(y=y)) 
}

if (x$cty=='diesel'){
 TL=x$TL;V1=x$V1;V2=x$V2;V3=x$V3
 if(is.null(x$S1)) S1=1 else S1 <- x$S1
  x1 <- list(V=c(V1,V2),T=c(TL,NA),S=c(S1,NA),path='adiabat',lab=c("1","2")) 
  y1 <- path.calc(x1); P2 <- y1$P[length(y1$P)]; S2<-S1 
  x2 <- list(V=c(V2,V3),P=c(P2,P2),S=c(S2,NA),path='isobar',lab=c("2","3"))
  y2 <- path.calc(x2); S3 <- y2$S[length(y2$S)]; T3 <- y2$T[length(y2$T)]; V4 <- V1
  x3 <- list(V=c(V3,V4),T=c(T3,NA),S=c(S3,NA),path='adiabat',lab=c("3","4"))
  y3 <- path.calc(x3); T4 <- y3$T[length(y3$T)];S4 <- S3
  x4 <- list(V=c(V4,V4),T=c(T4,TL),S=c(S4,NA),path='isochor',lab=c("4","1"))
  y4 <- path.calc(x4)
  z <- list(x1,x2,x3,x4)
  path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
  y <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
  return(list(y=y)) 
}

if (x$cty=='stirling'){
 TH=x$TH;TL=x$TL;V1=x$V1;V4=x$V4
 if(is.null(x$S1)) S1=1 else S1 <- x$S1
  x1 <- list(V=c(V1,NA),T=c(TL,TH),S=c(S1,NA),path='isochor',lab=c("1","2")) 
  y1 <- path.calc(x1); S2 <- y1$S[length(y1$S)];  V2<-V1; V3 <- V4
  x2 <- list(V=c(V2,V3),T=c(TH,NA),S=c(S2,NA),path='isotherm',lab=c("2","3"))
  y2 <- path.calc(x2);S3 <- y2$S[length(y2$S)]
  x3 <- list(V=c(V3,NA),T=c(TH,TL),S=c(S3,NA),path='isochor',lab=c("3","4"))
  y3 <- path.calc(x3);S4 <- y3$S[length(y3$S)]
  x4 <- list(V=c(V4,V1),T=c(TL,NA),S=c(S4,NA),path='isotherm',lab=c("4","1"))
  y4 <- path.calc(x4)
  z <- list(x1,x2,x3,x4)
  path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
  y.s <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
  return(list(y.summary=y.s)) 
}

if(x$cty=='box'){
 V1=x$V1;V3=x$V3;P1=x$P1;P2=x$P2
 x1 <- list(V=c(V1,NA),P=c(P1,P2),path='isochor',lab=c("1","2")) 
 y1 <- path.calc(x1)
 x2 <- list(P=c(P2,NA),V=c(V1,V3),path='isobar',lab=c("2","3"))
 y2 <- path.calc(x2)
 x3 <- list(V=c(V3,NA),P=c(P2,P1),path='isochor',lab=c("3","4"))
 y3 <- path.calc(x3)
 x4 <- list(P=c(P1,NA),V=c(V3,V1),path='isobar',lab=c("4","1"))
 y4 <- path.calc(x4)
 z <- list(x1,x2,x3,x4)
 path.lines(z,plane=plane,lab.cycle=TRUE,shade.cycle=shade.cycle)
 y <- list(path.summary(y1),path.summary(y2),path.summary(y3),path.summary(y4))
 return(list(y=y)) 
}

}

path.cycles.summary <- function(y){

 nvars <- length(y[[1]][[1]]$start.end[1,])
 pts <- y[[1]][[1]]$pts
 nM <- y[[1]][[1]]$nM
 Y1 <- data.frame(unname(y[[1]][[1]]$start.end)[-1,])
 Y2 <- data.frame(unname(y[[1]][[2]]$start.end)[-1,])
 Y3 <- data.frame(unname(y[[1]][[3]]$start.end)[-1,])
 Y4 <- data.frame(unname(y[[1]][[4]]$start.end)[-1,])
 end.state <- data.frame(rbind(Y1,Y2,Y3,Y4))
 row.names(end.state) <- c("1-2","2-3","3-4","4-1")
 names(end.state) <-c("v(m3/kg)","V(l)","P(b)","T(C)","s(J/Kkg)",
                      "S(J/K)","W(kJ/kg)","Q(kJ/kg)","cv(kJ/Kkg)","cp(kJ/Kkg)","cp/cv")

# WQnet <- data.frame(end.state[2,8],end.state[4,8],sum(end.state[,7]),sum(end.state[,8]),
#                     1-abs(end.state[4,8])/abs(end.state[2,8]))
# names(WQnet) <- c("Qin(kJ/kg)", "Qout(kJ/kg)","Wnet(kJ/kg)","Qnet(kJ/kg)","Eff")

 WQnet <- data.frame(end.state[2,8],end.state[4,8],round(1-abs(end.state[4,8])/abs(end.state[2,8]),2))
 names(WQnet) <- c("Qin(kJ/kg)", "Qout(kJ/kg)","Eff")

 return(list(WQnet=WQnet,end.state=end.state,pts=pts,nM=nM))

}



