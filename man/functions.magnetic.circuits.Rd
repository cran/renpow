\name{magnetic circuits}

\alias{magnetic.circuits}
\alias{inductor}
\alias{reluctance}
\alias{flux}

\title{Magnetic circuit calculations}

\description{Calculates magnetic elements and circuits}

\usage{
reluctance(x)
inductor(x)
flux(x)
}

\arguments{
  \item{x}{argument to reluctance or to inductance.
  For reluctance it is a list of one or more arrays; for instance, one array c(mu,l,A) where mu is relative permeability, l is path length, and A is cross-area.
  For inductance it is a list(N,rel) where N is number of turns and rel is reluctance
  For flux it is a list(N,i,rel)where N is number of turns, i is current, and rel is reluctance 
  }
}

\details{
Function reluctance calculates reluctance in MA-turn/Wb from mu relative permeability, l path length, and A cross-area.
Function inductor calculates inductance in mH from number of turns N and reluctance
Function flux calculates magnetic flux from number of turns, current, and reluctance 
}

\value{
\item{rel}{reluctance value from reluctance, inductor, or flux functions}
\item{prnt}{printout of reluctance,inductance, or flux including value and units}
\item{L}{inductance value from inductor or flux functions}
\item{mmf}{magnetomotive force or mmf}
\item{flux}{magnetic flux value}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)

Irwin, J.D. and R.M. Nelms. 2011. Basic Engineering Circuit Analysis. 11th edition. 2011: Wiley.
}

\author{Miguel F. Acevedo \email{Acevedo@unt.edu}}

\note{Functions used in Chapters 5 and 10 of Acevedo (2018)}

\examples{

# Chapter 5 simple inductor
mu=1000; l=1*10^-2; A=pi*(1*10^-3)^2
reluc <- reluctance(x=list(c(mu,l,A)))
reluc$prnt
ind <- inductor(x=list(N=159,rel=reluc$rel))
ind$prnt

# Chapter 10 air gaped inductor
mucore=1000; lcore=6*10^-2; Acore=1*10^-4
mugap=1; lgap=1*10^-3; Agap <- Acore
reluc <- reluctance(x=list(c(mucore,lcore,Acore),c(mugap,lgap,Agap)))
reluc$prnt

rel.eq <- reluc$rel[1]+reluc$rel[2] 
ind <- inductor(x=list(N=100,rel=rel.eq))
ind$prnt

magckt <- flux(x=list(N=10,i=1,rel=rel.eq))
magckt$prnt

magckt <- flux(x=list(N=10,i=1,r=reluc$r[2]))
magckt$prnt

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
