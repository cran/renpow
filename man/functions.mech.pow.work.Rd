\name{pow.work}

\alias{pow.work}
\alias{pow.work.plot}

\title{Mechanical power and work}

\description{Calculates and plots mechanical power and work.}

\usage{
pow.work(t, pow = "const", p)
pow.work.plot(x)
}

\arguments{
  \item{t}{Values of a time sequence for calculation of power}
  \item{pow}{Mode of power calculation, constant or linear}
  \item{p}{Value of power for constant or slope for linear}
  \item{x}{power object calculated from pow.work}
}

\details{
A call to pow.work requires three arguments: a time sequence, a type of power function such as constant or linear, and a value of power for the parameter of the function.
}

\value{
The output of pow.work contains the values of power and work as a function of time. Then we call the plot function using its output as argument to the pow.work.plot function
\item{t}{time sequence}
\item{p.w}{matrix containing with power and work}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 1 of Acevedo (2018)}

\examples{
t <- seq(0,60,0.1)
x <- pow.work(t,pow="const",p=2)
pow.work.plot(x)

x <- pow.work(t,pow="linear",p=0.025)
pow.work.plot(x)
}

\keyword{hplot}
\keyword{dplot}
\keyword{math} 
