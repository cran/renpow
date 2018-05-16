# magnetic elements and circuits

reluctance <- function(x){
 mu0 = 4*pi*10^-7
 nr <- length(x); r <- array()
 for (i in 1:nr)
 r[i] <- x[[i]][2]/(x[[i]][1]*mu0*x[[i]][3]) 
 # r in MA-turn/Wb
 r <- round(r*10^-6,2)
 prnt <- paste("Rel=",r," MA-turn/Wb")
 return(list(rel=r,prnt=prnt))
}

inductor <- function(x){
 # convert to A-turn/Wb
 rel <- x$rel*10^6
 L <- (x$N^2/rel)*10^3 # mH
 L <- round(L,2)
 prnt <- paste("L=",L,"mH")
 return(list(rel=rel,L=L, prnt=prnt))
 } 

flux <- function(x){
 mmf <- x$N*x$i
 rel.M <- x$rel # MA-turn/Wb
 rel <- x$rel*10^6 #A-turn/Wb
 L <- (x$N^2/rel)*10^3 # mH
 L <- round(L,2)

 flux <- round(mmf/rel.M,2)
 prnt <- paste("mmf=",mmf,"A-turn", "Rel=", rel.M,"MA-turn/Wb",
         "flux=",flux,"uWb","L=",L,"mH")
 return(list(mmf=mmf,rel=rel.M,flux=flux,L=L, prnt=prnt))
 } 

