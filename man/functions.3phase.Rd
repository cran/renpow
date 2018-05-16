\name{three-phase}

\alias{generator}

\title{Three-phase generator}

\description{Calculates three-phase generator}

\usage{
generator(x)
}

\arguments{
  \item{x}{parameters for generator: list(S3p, Vl.rms, pf, lead.lag, Zs.r) where: S3p rated three-phase apparent power, V1.rms line voltage as rms, pf load power factor, lead.lag is -1 or +1, and Zs.r armature impedance in rectangular form}
  }

\details{
Function generator produces E, Vp, and Il in polar form and a phasor diagram.
}

\value{
 \item{E.p}{electromotive force in polar form}
 \item{Vp.p}{Phase voltage generated in polar form}
 \item{Il.p}{Line current generated in polar form}
 \item{ys}{Source results}
 }

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 10 of Acevedo (2018)}

\seealso{AC functions \code{\link{ac.plot}}, \code{\link{phasor.plot}}, \code{\link{waves}}}

\examples{

x <- list(S3p = 15*10^6, Vl.rms = 13.8*10^3, pf=0.85, lead.lag=-1, Zs.r = c(0.1,2))
generator(x)

x <- list(S3p = 15*10^6, Vl.rms = 13.8*10^3, pf=0.85, lead.lag=1, Zs.r = c(0.1,2))
generator(x)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
