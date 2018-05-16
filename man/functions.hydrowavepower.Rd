\name{Wave Power}

\alias{wavepower}
\alias{generate.duration}
\alias{powflux.wave}
\alias{wave.contour}
\alias{duration.wave}
\alias{energy.wave}
\alias{energy.gen}
\alias{rotate}

\title{Wave Power Functions}

\description{Calculations in wave power generation}

\usage{

generate.duration(file)
powflux.wave(Hs, Tp)
wave.contour(X,label,sum, sumlabel)
duration.wave(datafile, file=TRUE)
energy.wave(Pflux, D)
energy.gen(Ew, L, nu)
}

\arguments{
  \item{datafile}{either a data file name or a renpow dataset}
  \item{file}{a logical to decide between a data file or a dataset}
  \item{Hs}{height in HsTp matrix}
  \item{Tp}{duration in HsTp matrix}
  \item{X}{Object to plot as contour; can be Pflux, Ew, or Eg}
  \item{Pflux}{Power flux produced by powflux.wave function}
  \item{Ew}{Energy produced by energy.wave function}
  \item{label}{label for plot; the default is ""}
  \item{sum}{logical default TRUE}
  \item{sumlabel}{label for sum; the default is ""}
  \item{D}{duration produced by duration.wave function}
  \item{L}{length of structure or energy conversion device}
  \item{nu}{efficiency}
 }

\details{Basic calculations for wave power}

\value{
 \item{Eg}{Energy produced by energy.gen function}
 }
 
\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 12 of Acevedo (2018)}

\examples{

 Pflux <- powflux.wave(Hs=seq(1,10),Tp=seq(5,20)) 
 wave.contour(X=Pflux,label="Power flux (kW/m)")
 D <-  duration.wave(WaveHsTp,file=FALSE)
 # alternatively using file in extdata 
 #  D <- duration.wave(system.file("extdata","WaveHsTp.csv",package="renpow"))
 wave.contour(X=D,label="Duration (hrs/yr)")

 Ew <- energy.wave(Pflux,D)
 wave.contour(X=Ew,label="Energy flux (MWh/yr/m)",sum=TRUE,sumlabel="(MWh/yr/m)")
 Eg <- energy.gen(Ew,L=30,nu=0.4)
 wave.contour(X=Eg,label="Energy produced (MWh/yr)",sum=TRUE,sumlabel="(MWh/yr)")

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
