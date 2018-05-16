\name{harmonics}

\alias{ACharmonics}
\alias{harmonic}

\title{AC Harmonic Distortion}

\description{Calculates harmonic distortion}

\usage{
harmonic(x, harm.odd, lab.units)
}

\arguments{
  \item{x}{list of arrays c(mag,phase) with fundamental magnitude and phase}
  \item{harm.odd}{fraction of odd harmonics with respect to the fundamental}
  \item{lab.units}{label for units}
  }

\details{
Function harmonic performs calculations and plots, including odd harmonics and THD.
}

\value{
 \item{t}{Time sequence}
 \item{Itot}{Total currents}
 \item{Isum}{Sum of currents in neutral}
 \item{THD}{Total Harmonic Distortion} 
 }

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 10 of Acevedo (2018)}

\seealso{AC waves and plots \code{\link{ac.plot}}, \code{\link{phasor.plot}}, \code{\link{waves}} three-phase \code{\link{generator}}}

\examples{

# single phase harmonics
x <- list(c(10,0)); harm.odd <- list(c(0.2,0.1,0.05)); lab.units <- "I [A]"
y <- harmonic(x,harm.odd,lab.units)

# three-phase harmonics
x <- list(c(10,0),c(10,-120),c(10,120))
harm.odd <- list(c(0.2,0.1,0.05),c(0.2,0.1,0.05),c(0.2,0.1,0.05))
lab.units <- "I [A]"
y <- harmonic(x,harm.odd,lab.units)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
