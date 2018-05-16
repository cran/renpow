\name{Wind Power}

\alias{windpower}
\alias{pow.rho.v3.table}
\alias{pow.rho.v3}
\alias{rho.pT.air}
\alias{rho.zT.air}
\alias{pow.v3.plot}
\alias{pow.wind}
\alias{v.H}
\alias{betz}
\alias{cal.vH}
\alias{weibull.plot}
\alias{cdf.plot}
\alias{fit.wind}
\alias{pow.class}
\alias{power.curve}
\alias{prob.power.curve}
\alias{wind.energy}

\title{Wind Power calculations and statistics}

\description{Basic calculations in wind power generation}

\usage{
pow.rho.v3.table(x)
pow.rho.v3(xw)
rho.pT.air(pT)
rho.zT.air(zT)
pow.v3.plot(x)
pow.wind(pw)
v.H(vh)
betz()
cal.vH(calvh)
weibull.plot(xmax,scale,shape)
cdf.plot(rv,xlab,ylab)
fit.wind(xd,vlabel)
pow.class(wc)
power.curve(pc)
prob.power.curve(pc,avg)
wind.energy(pc,Pow,avg)
}

\arguments{
  \item{x}{a list with possible list(rho= density, v=wind speed and A=cross-section area)}
  \item{xw}{list(rho=array of density values,v=array of wind speed values, A=cross-section area)}
  \item{pT}{a list(p,T.C,punit)}
  \item{zT}{list(z=1000, T.C=20, punit="bar", lapse=6)}
  \item{pw}{list with z, array T.C, array v, punit press unit,lapse,yleg,ylabel} 
  \item{vh}{list with array alpha exponent , and array rough of roughness coeffs} 
  \item{calvh}{list with list(v1.v2=v1.v2,H1=2,H2=10)} 
  \item{xmax}{max x value for Weibull}
  \item{scale}{parameter}
  \item{shape}{parameter}  
  \item{rv}{data var for cumulative distribution function}
  \item{xlab}{label x axis}  
  \item{ylab}{label y axis}  
  \item{wc}{wind speed for wind class chart}  
  \item{xd}{wind data object from two anemometers}
  \item{vlabel}{label for graph}  
  \item{pc}{a list cutin, vrated, cutout, A, v }
  \item{avg}{avg wind speed for prop.power.curve and wind.energy}  
  \item{Pow}{result of power curve and used as argument to wind energy}  
  }

\details{Basic calculations of wind power}

\value{
 \item{X1}{result of pow.rho.v3.table is list(X,P)}
 \item{X2}{result of pow.rho.v3 is list(rho=,v=,Pow=)}
 \item{X3}{result of rho.pT.air is array with c("Pressure(kPa)","Temp(C)","Density(kg/m3)")}
 \item{X4}{result of rho.zT.air is list(X,rho)}
 \item{X5}{result of v.H is v.v0.exp and v.v0.log}
 \item{X6}{result of cal.vH alpha exponent, and array rough of roughness coeffs}
 \item{X7}{result of pow.class is list(Pow.x,class.x)}
 \item{X8}{result of pow.curve is list(Pow, z,P.rated)}
 \item{X9}{result of prob.power.curve is hours in various ranges list(h.cutin,h.cutout,h.rated.nstop, h.run.below.rated)}
 \item{X10}{result of wind.energy is list(energy,unit,CF)}
 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 13 of Acevedo (2018)}

\examples{

x <- list(rho=1.225,v=10,A=1); pow.rho.v3.table(x)

x <- list(rho=c(0.9,1,1.1,1.225,1.3),v=seq(0,10),A=1)
X <- pow.rho.v3(x)

x <- list(v=X$v,y=X$rho,Pow=X$Pow,yleg="rho",ylabel="Density(kg/m3)")
pow.v3.plot(x)

x <- list(z=1000, T.C=10, punit="bar")
rho.zT.air(x)
x <- list(z=100, T.C=30, punit="bar")
rho.zT.air(x)

x <- list(z=1000, T.C=c(10,20,30), v=seq(0,30),punit="bar",
          lapse=6,yleg="T.C",ylabel="Temperature (C)")
pow.wind(x)

x <- list(alpha=c(0.1,0.25,0.4),rough=c(0.1,0.4,1.6))
v.H(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
