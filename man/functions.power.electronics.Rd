\name{AC-DC power electronics}

\alias{AC.DC.AC}
\alias{rectifier}
\alias{ac.plot.rect}
\alias{inverter}

\title{AC-DC and DC-AC in the time domain}

\description{Calculates and plots AC sinusoidal waves in the time domain and phasors. }

\usage{
rectifier(v.t, full = FALSE)
ac.plot.rect(V.t, v.lab = "v(t)", v.units = "V", y.lab = "v(t)[V]", rms = FALSE)
inverter(x)
}

\arguments{
  \item{v.t}{Values of AC variable at time intervals calculated using waves}
  \item{full}{Logical to use full wave rectifier}
  \item{V.t}{output of rectifier}
  \item{v.lab}{Label for variable in time domain plot or phasor plot; time domain would typically include (t) and phasor would be upper case} 
  \item{v.units}{Units for variable in time domain plot or phasor plot; time domain would typically include (t) and phasor would be upper case} 
  \item{y.lab}{Label for y axis composed of variable label and units}
  \item{rms}{Logical for whether the RMS value is addedd to the plot} 
  \item{x}{list(f,vin,nc): frequency, input voltage, number of cycles}
}
\details{Calculations of basic power electronics devices}
\value{
  \item{V.t}{output of rectifier list(w,t,nw,ym,ang,y,yrms,yavg)}
  \item{vtr}{Values of rectified AC variable at time intervals calculated using rectifier}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)

Irwin, J.D. and R.M. Nelms. 2011. Basic Engineering Circuit Analysis. 11th edition. 2011: Wiley.
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Used in Chapter 5 of Acevedo 2018}

\seealso{AC functions \code{\link{ac.plot}}, \code{\link{phasor.plot}}, \code{\link{waves}}}

\examples{

# rectified wave
v.AC <- list(c(170,0))
v.t <- waves(v.AC,nc=4) 
V.t <- rectifier(v.t)
ac.plot.rect(V.t)

# inverter
x <- list(f=60,vin=170,nc=2) 
inverter(x)

}
\keyword{hplot}
\keyword{dplot}
\keyword{math}

