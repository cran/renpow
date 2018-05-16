\name{CO2}

\alias{CO2data}
\alias{CO2monthly}
\alias{CO2annual}

\docType{data}

\title{Dataset: CO2 monthly and annualy}

\description{Atmospheric CO2 monthly and annualy}

\usage{
data("CO2monthly")
data("CO2annual")
}
\format{
CO2 monthly is CO2 expressed as a mole fraction in dry air, micromol/mol, abbreviated as ppm
a record has: year, month,decimal-date, average, interpolated, trend-season-corr, days
For instance: 1958   3    1958.208      315.71      315.71      314.62     -1
CO2 annual has three values
Year, Avg, Uncorrected
}
\details{CO2 monthly and annual}

\source{NOAA. Trends in Atmospheric Carbon Dioxide. 2017. NOAA, Earth System Research Laboratory, Global Monitoring Division. Accessed June 2017.
URL: http://www.esrl.noaa.gov/gmd/ccgg/trends/}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}
\examples{
data(CO2monthly)
}
\keyword{datasets}
