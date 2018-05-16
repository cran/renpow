\name{ERCOT}

\alias{ERCOT}
\alias{ERCOT2010}
\alias{ERCOT2010.01}

\docType{data}

\title{Dataset: ERCOT demand for year 2010}

\description{ERCOT load year 2010 by region}

\usage{
data("ERCOT2010")
data("ERCOT2010.01")
}
\format{
Day and Hour, follwoed by demand for 8 regions, and demand total ERCOT
Header:   
Hour_End	COAST	EAST	FAR_WEST	NORTH	NORTH_C	SOUTHERN	SOUTH_C	WEST	ERCOT
Example of record: 
1/1/2010 1:00	7775.679969	1238.21679	1237.781607	950.698047	12406.20808	2467.652377	5032.076127	1059.772842	32168.08584
}
\details{Example demand grid}
\source{ERCOT. Hourly Load Data Archives. 2017. Accessed October 2017. Available from: http://www.ercot.com/gridinfo/load/load_hist}
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}
\examples{
data(ERCOT2010.01)
}
\keyword{datasets}
