\name{Hydro Power}

\alias{hydropower}
\alias{P.hA}
\alias{P.Qh}
\alias{Pmax.Qh}
\alias{Pe.Pw}
\alias{Pmax.Qh.plot}
\alias{turbine.regions}
\alias{turbine.regions.all}
\alias{pipe.loss}

\title{Hydroelectric Power Functions}

\description{Basic calculations in hydroelectric power generation}

\usage{
P.hA(x)
P.Qh(x)
Pmax.Qh(x)
Pe.Pw(x)
Pmax.Qh.plot(x)
turbine.regions(type)
turbine.regions.all()
pipe.loss(pipe)
}

\arguments{
  \item{x}{a list with possible h, A, Q, nu, plab: head, cross-area, flow, efficiency,label}
  \item{type}{a turbine type from 'kaplan', 'francis', 'pelton', 'crossflow', 'slh'}
  \item{pipe}{a list Q, d, L, mat: flow, diameter, length, and material. Component mat is one of 'pvc', 'concrete', 'steel', 'galvanized', 'poly'}
  }

\details{Basic calculations of hydropower}

\value{
 \item{X}{result of P.hA(x) is array with Head(m),Vel(m/s),Area(m2),Flow(m3/s), and Power(in kW or MW)}
 \item{X}{result of P.Qh(x) is array with Head(m),Flow(m3/s),and Power(in kW or MW)}
 \item{X}{result of Pmax.Qh(x) is array with Gross head (m),Net head (m),,Flow(m3/s),and Power(in kW or MW)}
 \item{X}{result of Pe.Pw(x) is array with GrossHead(m),NetHead(m),Flow(m3/s),Press(kPa),Eff, PowWater9in kW or MW), and PowGen(in kW or MW)}
 \item{X}{result of pipe.loss(pipe) is array with Head loss(m), and Roughness}

 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 12 of Acevedo (2018)}

\examples{

# head 3 m and cross-sectional area of 2 m2
x <- list(h=3,A=2); P.hA(x) 

x <- list(Q=1000,h=15); Pmax.Qh(x)

x <- list(Q=1000,h=15,nu=0.9); Pe.Pw(x)

x <- list(h=1,Q=100); P.Qh(x)

x <- list(Q=1000,h=15,plab="A"); Pmax.Qh.plot(x)
turbine.regions(type='francis')

x <- list(Q=0.01,d=0.075,L=200,mat='pvc')
hL <- pipe.loss(x)[1,1]
x <- list(h=30-hL,Q=0.01,nu=0.9)
P.Qh(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
