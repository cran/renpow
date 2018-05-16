\name{tau and bd}

\alias{tau.bd}
\alias{tauGolden}
\alias{tauAtlanta}

\docType{data}

\title{Dataset: monthly tau and bd for a location}

\description{Optical depth for a location}

\usage{
data("tauGolden")
data("tauAtlanta")
}

\format{
 list with location, lat-long-elev, and two columns for tau (b, d)
}
\details{Optical depth}

\source{
ASHRAE. ASHRAE Climate Data Center. 2017. Accessed November 2017.
Available from: https://www.ashrae.org/resources--publications/bookstore/climate-data-center#std169.
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}
\examples{
tau <- tauGolden$tau
}
\keyword{datasets}
