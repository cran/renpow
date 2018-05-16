\name{grid}

\alias{grid}
\alias{infinite.bus}
\alias{seg.ts}
\alias{series.plot}
\alias{daily.min.max}
\alias{load.duration}

\title{Electrical Power Grid}

\description{Electric power systems calculations including buses and loads}

\usage{
infinite.bus(x)
seg.ts(Xts,dh0,dhf,var) 
series.plot(x.t)
daily.min.max(x.t.y, ylabel, inst.cap, inst.lab)
load.duration(x.t.y, inst.cap, inst.lab)
}

\arguments{
  \item{x}{parameters for infinite bus a list(V,X,E,delta):V Voltage, X Reactance, E array of electromotive force, delta array of angle delta}  
  \item{Xts}{data frame with ts values}
  \item{dh0}{intial time stamp}
  \item{dhf}{final time stamp}
  \item{var}{variable number to select}
  \item{x.t}{time series to plot}
  \item{x.t.y}{load time series for daily.min.max plots and load-duration curve}
  \item{ylabel}{label for min max} 
  \item{inst.cap}{array of installed capacities}
  \item{inst.lab}{array of labels for installed capacities}
  }

\details{
infinite.bus uses a fixed V and X, and plots the effect of a set of E and angle delta
daily.min.max draws a plot showing daily maximum and minimum along with horizontal lines for installed capacities 
load.duration builds a load-duration curve and calculate CF for daily load data and a set of generation capacities   
}

\value{
 \item{Hours}{number of hours for each class}
 \item{energy}{total energy under curve}
 \item{possible}{possible utilization for each class}
 \item{utilization}{actual utilization for each class}
 \item{CF}{capacity factor for each class}
 }

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 11 of Acevedo (2018)}

\examples{

x <- list(V=13.8,X=1,E=seq(14,16,1),delta=seq(6,30,1))
infinite.bus(x)

X <- ERCOT2010
x.t <- seg.ts(X,dh0="1/1/2010 1:00",dhf="1/1/2011 0:00",c(1:9))$x.t

# week day Monday Jan 04 2010
x.t.wd.wt <- seg.ts(X,dh0="1/4/2010 0:00",dhf="1/4/2010 23:00",c(1,5,7))$x.t.seg
series.plot(x.t.wd.wt) 

# week Sunday Jan 24 to Saturday Jan 30 2010
x.t.wd.wt <- seg.ts(X,dh0="1/24/2010 0:00",dhf="1/30/2010 23:00",c(1,5,7))$x.t.seg
series.plot(x.t.wd.wt)

# NORTH_C
x.t.y <- x.t[,5]; ylabel <- colnames(x.t)[5]; inst.cap <- c(10000,8000,8000) #MW
inst.lab <- c("Baseload","Intermediate", "Peaking")
daily.min.max(x.t.y,ylabel,inst.cap,inst.lab)

x.t.y <- x.t[,5]
inst.cap <- c(10000,8000,8000) # MW
inst.lab <- c("Baseload","Intermediate", "Peaking")
load.duration(x.t.y, inst.cap,inst.lab)

}

\keyword{hplot}
\keyword{dplot}
\keyword{math}
