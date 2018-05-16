\name{transformer}

\alias{transformer}
\alias{xformer.ckt}

\title{Transformer Circuits}

\description{Calculates tranformer circuits}

\usage{
xformer.ckt(x, dig = 2)
}

\arguments{
  \item{x}{list(N,Vs.p,Zs.r,Zo.r,Zl.r): these are N turns ratio, Vs.p Voltage source in polar form, Zs.r Source impedance in rectangular form, Zo.r Output impedance in rectangular form, Zl.r Load impedance in rectangular form.}
  \item{dig}{number of digits for results}
  }

\details{
xformer.ckt transformer circuit calculations using source, impedances, and load,
}

\value{
 \item{ys}{Source results}
 \item{yx}{Transformer results}
 \item{yl}{Load results}
 \item{prnt}{Results of transformer function for printout form}
 }

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 10 of Acevedo (2018)}

\seealso{AC functions \code{\link{ac.plot}}, \code{\link{phasor.plot}}, \code{\link{waves}}}

\examples{

x <- list(N=c(2,1),Vs.p=c(12,30),Zs.r=c(1,0),Zo.r=c(0,0),Zl.r=c(1,0))
xf <- xformer.ckt(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
