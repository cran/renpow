\name{AC power}

\alias{ACpower}
\alias{inst.pow.calc}
\alias{inst.pow.plot}
\alias{inst.pow.leg}
\alias{complex.pow.calc}
\alias{complex.pow.plot}
\alias{complex.pow.tri}
\alias{pf.corr}
\alias{pf.corr.tri}

\title{AC power in the time and frequency domain}

\description{Calculates and plots AC power in the domain and phasors. Includes complex number calculations.}

\usage{
inst.pow.calc(x, freq = 60, nc = 2)
inst.pow.plot(x, rms = FALSE, freq = 60, nc = 2)
inst.pow.leg(ang, lab, ym, w, units, pf)
complex.pow.calc(xc, dig = 2, res = TRUE)
complex.pow.plot(cp)
complex.pow.tri(cp)
pf.corr(P, V, pf, pfc, w = 377, dig = 2)
pf.corr.tri(xpfc)
}

\arguments{
  \item{x}{list of voltage and current phasors}
  \item{freq}{frequency in Hz}
  \item{nc}{number of cycles to calculate and/or plot}
  \item{rms}{logical to decide to plot the rms value}
  \item{ang}{angles for legend}
  \item{lab}{label for legend}
  \item{ym}{magnitude for legend}
  \item{w}{angular frequency for legend and for pf correction}
  \item{units}{units for legend}
  \item{pf}{power factor for legend of inst.pow.leg or to be corrected by pf.corr}  
  \item{xc}{list with voltage, current, and phase angle}
  \item{dig}{number of digits to round complex power}
  \item{res}{logical to print results of complex power} 
  \item{cp}{output of complex power; a list}
  \item{P}{real power for pf correction}
  \item{V}{voltage for pf correction}
  \item{pfc}{target power factor}
  \item{xpfc}{output of pf.corr; a list}  
  }

\details{
inst.pow.calc calculates time domain isntantaneous power given current and voltage phasors. Its output is passed to inst.pow.plot.
inst.pow.plot produces time domain plots given current and voltage phasors. It uses inst.pow.leg to display legend
 
}

\value{
\item{w}{angular frequency in rad/s}
\item{vm}{array of magnitude}
\item{ang}{array of angle}
\item{pf}{power factor}
\item{pavg}{average power}
\item{t}{time sequence}
\item{v.i}{matrix with voltage and current as a function of time}
\item{p}{instantaneous power or power as a function of time}
\item{units}{units for complex power result}
\item{S}{apparent power}
\item{theta}{phase angle}
\item{P}{real power}
\item{Q}{reactive power}
\item{pf}{power factor}
\item{prnt}{string to print results}
\item{cp}{complex power output of pf.corr}
\item{cpc}{complex power after correction}
\item{y}{Voltage,current, and corrected params}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)

Irwin, J.D. and R.M. Nelms. 2011. Basic Engineering Circuit Analysis. 11th edition. 2011: Wiley.
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter Chapter 8 of Acevedo (2018)}

\seealso{Power electronics functions \code{\link{ac.plot.rect}}, \code{\link{rectifier}}, \code{\link{inverter}}}

\examples{

# power resistor
vm=170; R=10
x <- list(c(vm,0),c(vm/R,0))
inst.pow.plot(x)

# power capacitor
w <- 377; v.s <- c(170,0)
C=1000*10^-6
# current response
i.res <- c(v.s[1]*(w*C),v.s[2]+90)
x <- list(v.s,i.res)
inst.pow.plot(x)

# calc complex power
V.s=c(170,10); Z.p=c(10,20)
I.p <- div.polar(V.s,Z.p)
V <- V.s[1]/sqrt(2); I <- I.p[1]/sqrt(2)
theta <- V.s[2]-I.p[2] 
cp <- complex.pow.calc(list(V,I,theta))

# pf correction
P=5; V=240; I=40; pfc=0.9
pf <- P*1000/(V*I)
# call pf correction function
pfcorr <- pf.corr(P,V,pf,pfc)
# visualize
pf.corr.tri(pfcorr)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
