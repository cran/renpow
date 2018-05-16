\name{TEO at DP 2017}

\alias{TEODP2017}

\docType{data}

\title{Dataset: Example of wind speed at two heights}

\description{Example summer 2017 TEO station at DP campus. Includes wind speed at two heights, air temp two heigths, BP, RH, Wind dir}

\usage{data("TEODP2017")}

\format{
 Text file with header x specifying the variable name 
}
\details{Example for wind speed height calibration}

\source{http://teo.unt.edu}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}
\examples{
v1.v2 <- TEODP2017[,c(1,8,4)]
}
\keyword{datasets}
