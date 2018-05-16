\name{sequestration}

\alias{Cseq}
\alias{forest.seq}

\title{Carbon Sequestration}

\description{Simulates carbon sequestration}

\usage{
forest.seq(x, y)
}

\arguments{
  \item{x}{parameters for forest growth; a list(t,B0,Bmax,nu,r): t time sequence,B0 Biomass at time zero, Bmax Maximum Biomass (parameter of Richards equation), nu Coefficient exponent to adjust shape (parameter of Richards equation), r Growth rate (parameter of Richards equation)}
  \item{y}{Specifications of power plant emissions a list(kgCO2.kWh,P,C): kgCO2.kWh kg of CO2 per kWh produced, P Power produced, C carbon content}
}

\details{
forest.seq produces rates, area, and a plot of C forest sequestration dynamics 
}

\value{
 \item{yr.max}{time in years to achieve biomass for max rate}
 \item{B.i}{required biomass in t/ha to sequester at max rate}
 \item{max.dCO2}{max rate}
 \item{tCO2.emiss}{carbon emission sequestered}
 \item{area}{area required}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 9 of Acevedo (2018)}

\examples{

forest <- list(t=seq(0,100),B0=1,Bmax=150,nu=0.5,r=0.3)
plant <- list(kgCO2.kWh=0.4,P=10^9,C=1.0) 
forest.seq(x=forest,y=plant)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
