\name{AC waves and plots}

\alias{ACcircuits}
\alias{waves}
\alias{ac.plot}
\alias{phasor.plot}
\alias{polar}
\alias{recta}
\alias{mult.polar}
\alias{div.polar}
\alias{horiz.lab}
\alias{wave.leg}
\alias{phas.leg}
\alias{gridcir}
\alias{sinplot}
\alias{arc}
\alias{rot.fig}
\alias{admit}
\alias{vector.phasor}

\title{AC in the time and frequency domain}

\description{Calculates and plots AC sinusoidal waves in the time domain and phasors. Includes complex number calculations.}

\usage{
waves(x, f = 60, nc = 2)
ac.plot(v.t, v.lab = "v(t)", v.units = "V", y.lab = "v(t)[V]", rms = FALSE)
phasor.plot(v.p, v.lab = "V", v.units = "V", lty.p = NULL)
polar(rec)
recta(pol)
mult.polar(x1,x2)
div.polar(x1,x2)
horiz.lab(nw, ym, tmax, ymax, units, yrms, rms)
wave.leg(nw, ang, lab, ym, w, units)
phas.leg(np, mag, ang, lab, units, lty.p)
gridcir(rmax)
sinplot(xlab, ylab)
arc(mag, ang)
rot.fig(vp,v.lab="wt")
admit(Z.r)
vector.phasor(V, I)
}

\arguments{
  \item{x}{AC variable given as a list of arrays. Each array contains two entries: magnitude and phase}
  \item{f}{Frequency in Hz, default 60 Hz}
  \item{nc}{Number of cycles to calculate and plot, default 2 cycles}
  \item{v.t}{Values of AC variable at time intervals calculated using waves}
  \item{v.lab}{Label for variable in time domain plot or phasor plot; time domain would typically include (t) and phasor would be upper case; also specifies label for angle in rot.fig} 
  \item{v.units}{Units for variable in time domain plot or phasor plot; time domain would typically include (t) and phasor would be upper case} 
  \item{y.lab}{Label for y axis composed of variable label and units}
  \item{rms}{Logical for whether the RMS value is added to the plot} 
  \item{v.p}{Phasors to plot}
  \item{lty.p}{set of line types for plot and type of line for legend} 
  \item{rec}{argument to polar array of rectangular coord}
  \item{pol}{argument to rect: array of magnitude and angle}
  \item{x1,x2}{complex numbers to multiply or divide by mult.polar and div.polar}
  \item{nw}{number of waves to assign horizontal lines}
  \item{ym}{magnitude array}
  \item{tmax}{max x axis value}
  \item{ymax}{max y axis value}
  \item{units}{units for lines and for legend}
  \item{yrms}{rms value to be used in y axis}
  \item{ang}{angle of waves or phasors}
  \item{lab}{label for legend}
  \item{w}{angular frequency for legend}
  \item{np}{number of phasors for legend}
  \item{mag}{magnitude of waves or phasors}
  \item{rmax}{max extent of polar grid for gridcir} 
  \item{xlab}{xaxis label for sinplot}
  \item{ylab}{yaxis label for sinplot}
  \item{vp}{voltage for rot fig}
  \item{Z.r}{impedance in rectangular form}
  \item{V}{Voltage}
  \item{I}{Current}
  }

\details{
Waves is first used to calculate values for time and the variable. Then the object creted by wave is used by ac.plot.
Functions polar() and recta() allow polar and rect conversion specifying phasors as arrays.
Function mult.polar and div.polar are used for multiplication and division of phasors.
Function horiz.lab plots horizontal lines with labels for magnitude and rms. Used mostly by other renpow functions.
Function wave.leg and phas.leg write out waves and phasors for legend. Used mostly by other renpow functions.
Function gridcir draws a polar grid. Used mostly by other renpow functions. 
Function arc draws an arc from 0 degrees to a phasor line specified by mag and ang. Used mostly by other renpow functions. 
Function admit calculates admittance given the impedance in rectagular form.
Function vector.phasor specifies volatge and current phasors for plotting.
}

\value{
Function waves:
\item{w}{angular frequency}
\item{t}{time values}
\item{nw}{number of waves}
\item{ym}{magnitude values}
\item{ang}{ang}
\item{y}{values}
\item{yrms}{rms values}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)

Irwin, J.D. and R.M. Nelms. 2011. Basic Engineering Circuit Analysis. 11th edition. 2011: Wiley.
}

\author{Miguel F. Acevedo \email{Acevedo@unt.edu}}

\note{Functions used in Chapter 5, 8, and 10 of Acevedo (2018)}

\seealso{Power electronics functions \code{\link{ac.plot.rect}}, \code{\link{rectifier}}, \code{\link{inverter}}}

\examples{

# from Chapter 5
# one wave show phase angle
x <- list(c(170,30)); v.t <- waves(x); ac.plot(v.t)

# two waves different magnitude and phase
x <- list(c(170,0),c(160,30)); v.t <- waves(x);
v.lab <- c("v1(t)","v2(t)"); v.units <- rep("V",2)
ac.plot(v.t,v.lab,v.units)

# one phasor
v.p <- list(c(170,10)); phasor.plot(v.p)

# phasors phase difference
v.units <- rep("V",2)
v.lab <- c("V1","V2")
# V1 leads V2
v.p <- list(c(170,70),c(170,50))
phasor.plot(v.p,v.lab,v.units)

# rect to polar
polar(c(2,1))
# polar to rect
recta(c(2,45))
# multiplication
x <- polar(c(1,2))
y <- polar(c(2,3))
mult.polar(x,y)

# from Chapter 8
# nodal analysis
Y1 <- 1/(5+5i); Y2 <- 1/(5+5i); Y3 <- 1/(10+10i) 
Y <- matrix(c(Y1+Y2,-Y2,-Y2,Y3+Y2),ncol=2,byrow=TRUE)
Is <- c(1+0i,0+0i)
Vn <- solve(Y,Is)
VpIp <- vector.phasor(Vn,Is)
phasor.plot(VpIp$VI, c("V1","V2","Is1","Is2"),
           c("V","V","A","A"),lty.p=c(2,2,1,1))

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
