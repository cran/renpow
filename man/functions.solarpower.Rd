\name{Solar Power}

\alias{solarpower}
\alias{I0.orbit}
\alias{days.mo}
\alias{declination}
\alias{sun.elev}
\alias{read.tau}
\alias{beam.diffuse}
\alias{I0.blackbody}
\alias{spectral}
\alias{useful.waste}
\alias{sun.path}
\alias{sun.diagram}
\alias{collector}
\alias{month.prod}
\alias{tilt.adj}
\alias{one.axis.tracking}
\alias{two.axis.tracking}

\title{Solar Power Functions}

\description{Basic calculations in solar power generation including irradiance and sun path}

\usage{
I0.orbit(leap = FALSE, plot = TRUE)
days.mo(day, leap = FALSE)
declination(leap=FALSE,plot=TRUE)
sun.elev(xdec,lat,plot=TRUE)
read.tau(file)
beam.diffuse(dat,plot=TRUE)
I0.blackbody(T.sun, wl.nm, plot = TRUE)
spectral(X, label, wl.lim.nm, T.sun, plot.surf = FALSE)
useful.waste(I0.bb)
sun.path(lat, nday, plot = TRUE)
sun.diagram(lat)
collector(Ibd, sunpath, tilt, azi.c, fr, label = "")
month.prod(dat)
tilt.adj(lat, days, labels)
one.axis.tracking(dat, mode = "PNS")
two.axis.tracking(dat)
}

\arguments{
  \item{leap}{logical to use a leap year}
  \item{plot}{logical to plot}
  \item{day}{day of the year for a given day in the month}
  \item{xdec}{result from declination}
  \item{lat}{latitude}
  \item{file}{filename}  
  \item{T.sun}{temp sun}
  \item{wl.nm}{wavelength in nm}
  \item{X}{data read from file}
  \item{label}{label for plot}
  \item{wl.lim.nm}{limits for wavelength}
  \item{plot.surf}{logical plot together with surface data}
  \item{dat}{data from file}
  \item{I0.bb}{output from I0.blackbody}
  \item{nday}{day number}
  \item{Ibd}{Irradiance: direct beam and diffuse}
  \item{sunpath}{elev and azimuth from the output of sun.path}
  \item{tilt}{collector tilt angle}
  \item{azi.c}{orientation of collector}
  \item{fr}{fraction reflected}
  \item{days}{days}
  \item{labels}{labels}
  \item{mode}{tracking mode}
  }

\details{Basic calculations of solar power}

\value{
 \item{I0}{result of I0.orbit: values of daily ET solar irradinace}
 \item{day.mo}{result of day.month: day numbers beginning at day}
 \item{dec}{result of declination: ne equinox day and values of daily declination}
 \item{elev}{result of sun elev: ne equinox day, latitude, and and values of daily sun elev}
 \item{tau}{list(loc, lat.long.elev, tau)} 
 \item{Ibd}{list(day21.mo,tau,air.mass,lat,Ib,Id,Id.Ib)} 
 \item{I0.bb}{list lambda,I.sun.nm,I0}
 \item{sunpath}{list nday, hr.noon, azi, elev}
 \item{Ic}{list(Ib,Id,Ibc,Idc,Irc,Ic,I.h)}
 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 14 of Acevedo (2018)}

\seealso{DC circuits \code{\link{PVcell}}, \code{\link{PVcell.plot}}}

\examples{

I0.orbit()
declination()

x <- declination(plot=FALSE)
sun.elev(x,lat=32.9)

tau <- tauGolden
# or alternatively
# tau <- read.tau(system.file("extdata","tauGolden.csv",package="renpow"))
Ibd <- beam.diffuse(tau)

I0.blackbody(T.sun=5800,wl.nm=seq(150,2500))






}

\keyword{ hplot }
\keyword{ dplot }
\keyword{ math }
