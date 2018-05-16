\name{panels plots}

\alias{panels.plots}
\alias{panels}
\alias{plot2yaxis}

\title{Plot utilities}

\description{Uses layout to plot several graphs on a page and plot on two vertical axis}

\usage{
panels(wd, ht, rows, cols, pty, int = "r")
plot2yaxis(x, x0 = FALSE, y0 = FALSE)
}

\arguments{
  \item{wd}{Width}
  \item{ht}{Height}
  \item{rows}{number of rows}
  \item{cols}{number of columns}
  \item{pty}{type}
  \item{int}{axis intersection}
  \item{x}{Variable x for horizontal axis and variables y1, y2 for vertical axes}
  \item{x0}{Logical to select axis limit}
  \item{y0}{Logical to select axis limit}
}

\details{Function panels makes more efficient use of space than par(mfrow)}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{Acevedo@unt.edu}}

\note{Util functions used in most chapters of Acevedo (2018)}

\examples{

wd=6; ht=3
panels(wd,ht,1,1,pty="m")
t <- seq(0,60,0.1)
x <- pow.work(t,pow="const",p=2)
pow.work.plot(x)

wd=7; ht=3.5
panels(wd,ht,1,2,pty="m")
# diode
V <- seq(-0.1,0.05,0.01) # volts
x <- diode(V)
ivplane(x)
V <- seq(-0.1,0.6,0.01) # volts
x <- diode(V)
ivplane(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}

