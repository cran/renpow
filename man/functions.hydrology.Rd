\name{Hydro Power Hydrology}

\alias{hydrology}
\alias{area.vol}
\alias{model.flow}
\alias{flow.plot}
\alias{flow.exc.plot}
\alias{exceed}
\alias{annual.avg}

\title{Hydro Power Hydrology functions}

\description{Hydrological calculations in hydroelectric power generation}

\usage{
area.vol(xav)
model.flow(mf)
flow.plot(flow, label)
flow.exc.plot(flow, exc, label)
exceed(flow)
annual.avg(mf, nyrs)
}

\arguments{
  \item{xav}{area.vol is a list(H,B,W,L) H pool elevation (m), B bottom elevation (m), W width (km), L tail length (km)}
  \item{mf}{list(base.flow,peak.flow,day.peak,length.season,variab,coef)}
  \item{flow}{flow time series resulting from model.flow}
  \item{label}{label for flow}
  \item{exc}{ exceedance levels as a result of exceed}
  \item{nyrs}{number of years}
  }

\details{Basic hydrological calculations for hydropower}

\value{
 \item{X}{list(y,proby, prob, Q, Prob.Qmean, prob.Q)}
 \item{flow}{flow time series from model.flow}
 \item{Xtm}{annual avg of flow time series from model.flow}
 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 12 of Acevedo (2018)}

\examples{

x <- list(H=130,B=100,W=10,L=100)
area.vol(x)

x <- list(base.flow=20,peak.flow=100,day.peak=200,length.season=90,
          variab=c(0.01,2),coef=c(0.4,0.3,0.2,0.1))
flow <- model.flow(x)
flow.plot(flow,label="Simulated flow (m3/s)")
exc <- exceed(flow)
exc$prob.Q
flow.exc.plot(flow,exc,label="Simulated flow (m3/s)")

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
