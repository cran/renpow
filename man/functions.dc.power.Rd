\name{DC circuits}

\alias{DCcircuits}
\alias{resistor}
\alias{ivplane}
\alias{diode}
\alias{vsource}
\alias{isource}
\alias{PVcell}
\alias{PVcell.plot}
\alias{eff.pow}
\alias{transient}
\alias{fuel.cell}

\title{DC circuits and power calculations}

\description{Calculates and plots DC V-I response as well as transient response }

\usage{
resistor(V,R)
ivplane(x, x0 = FALSE, y0 = FALSE)
diode(V)
vsource(Voc,Rs)
isource(Isc,Rp)
PVcell(x.PVcell)
PVcell.plot(y.PVcell)
eff.pow(x.eff.pow)
transient(ys, tau, ylabel, yslabel)
fuel.cell(x.fcell)
}

\arguments{
  \item{V}{Voltage}
  \item{R}{Resistance}
  \item{x}{Result of function resistor, diode, vsource, or isource}
  \item{x0}{Logical argument to plot2axis}
  \item{y0}{Logical argument to plot2axis}
  \item{Voc}{Open circuit voltage}
  \item{Rs}{Series resistance} 
  \item{Isc}{Short circuit current}
  \item{Rp}{Parallel resistance}
  \item{x.PVcell}{list(I0.A, Isc.A, Area, Rs, Rp, Light)}
  \item{y.PVcell}{output of PVcell}  
  \item{x.eff.pow}{list(Rth,Voc)}
  \item{ys}{Source for transient calculation}
  \item{tau}{Time constant}
  \item{ylabel}{y axis label for transient plots}
  \item{yslabel}{Source label for transient plots}
  \item{x.fcell}{list(area.cm2,Rload.ohm)}
  
  }

\details{Object x contains all the calculated values which can be passed to function ivplane to plot the graphs. 
}

\value{
\item{V}{Voltage}
\item{I}{Current}
\item{P}{Power}
\item{txt}{labels}
\item{Light}{Light levels from PVcell}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)

Irwin, J.D. and R.M. Nelms. 2011. Basic Engineering Circuit Analysis. 11th edition. 2011: Wiley.
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 3 and Chapter 5 of Acevedo (2018)}

\examples{

# from Chapter 3
# resistors
V <- seq(-0.2,1,0.01) # volts
x <- resistor(V,R=1)
ivplane(x)

# diode
V <- seq(-0.1,0.6,0.01) # volts
x <- diode(V)
ivplane(x)

# voltage source
x <- vsource(Voc=24,Rs=1)
ivplane(x)
# current source
x <- isource(Isc=20,Rp=1)
ivplane(x)

# PV cell
x <- list(I0.A=1, Isc.A=40, Area=100, Rs=0.05, Rp=1, Light=1) 
# units: I0.A pA/cm2 Isc.A mA/cm2 Area cm2 Rs  ohm Rp  ohm
X <- PVcell(x)
ivplane(X,x0=TRUE,y0=TRUE)

# efficiency of power transfer
x <- list(Rth=1,Voc=24)
eff.pow(x)

# from Chapter 5
# RC circuit transient
R=0.2;C=1 # Mohm and uF
transient(ys=12,tau=R*C,ylabel="Vc(t) [V]",yslabel="Vs [V]") 

# RL circuit transient
R=20;L=1; # kohm and mH
transient(ys=12/R,tau=L/R,ylabel="iL(t) [A]",yslabel="Vs/R [V]") 

# Chapter 6 fuel cells
x <- list(area.cm2=15,Rload.ohm=NA)
fuel.cell(x)
x <- list(area.cm2=15,Rload.ohm=0.5)
fuel.cell(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
