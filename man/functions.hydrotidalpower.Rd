\name{Tidal Power}

\alias{tidalpower}
\alias{read.tide}
\alias{harmonics.tide}
\alias{power.barrage.cycle}
\alias{find.peaks}
\alias{tide.current.abs}
\alias{tidal.power}

\title{Tidal Power Functions}

\description{Calculations in tidal power generation}

\usage{
read.tide(file)
harmonics.tide(x, days, ylabel, plot)
power.barrage.cycle(xba)
find.peaks(tz, band)
tide.current.abs(tz, ylabel, plot)
tidal.power(tz, Aflow)
}

\arguments{
  \item{file}{filename for tide files in extdata}
  \item{x}{a tidal dataset read from file}
  \item{days}{number of days}
  \item{ylabel}{label for y axis of plot default for harmonics.tide is "Tide wrt MSL (m)" default for tide.current.abs is "Current abs (m/s)" }
  \item{plot}{logical to decide to plot default is TRUE}
  \item{xba}{list(a,Abasin,z,nu): a center of mass at half the tidal range, Abasin tidal basin area, z tidal cycle range, nu efficiency}
  \item{tz}{time and tide height in tidal signal produced by harmonics.tide or rectified by tide.current.abs}
  \item{band}{threshold band to find peaks}
  \item{Aflow}{cross section area for tidal power generation}
  }

\details{Basic calculations for tidal power}

\value{
 \item{X}{list(t,z,x) tide time series}
 \item{X}{list(xp,tp,range)}
 \item{X}{list(pow.tide.MW,pow.gen.MW,gen.MWh)}
 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 12 of Acevedo (2018)}

\examples{

# using extdata
x <- read.tide(system.file("extdata","AnchorageTide.csv",package = "renpow"))
harmonics.tide(x,days=29)

x <- read.tide(system.file("extdata","VelocityTide.csv",package="renpow"))
y <- harmonics.tide(x, days=15,ylabel="Velocity m/s)")
tide.current.abs(y, ylabel="Velocity (abs val) (m/s)", plot=TRUE)
y <- harmonics.tide(x, days=365,ylabel="Velocity m/s)",plot=FALSE)
z <- tide.current.abs(y, ylabel="Current abs (m/s)", plot=TRUE)
tidal.power(z,Aflow=1)

x <- read.tide(system.file("extdata","ElevationTide.csv",package="renpow"))
y <- harmonics.tide(x,days=29)
y <- harmonics.tide(x,days=365)
z <- find.peaks(y, band=c(0,1))

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
