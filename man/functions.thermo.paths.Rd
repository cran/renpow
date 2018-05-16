\name{thermo paths}

\alias{thermodynamics}
\alias{cpcv.cal}
\alias{RefCoefAir}
\alias{cp.cv}
\alias{fcp}
\alias{fcpcv}
\alias{fcv}
\alias{fpv}
\alias{simpson}
\alias{simpson.pv}
\alias{path.calc}
\alias{path.summary}
\alias{isochor}
\alias{isobar}
\alias{adiabat}
\alias{isotherm}
\alias{phase}
\alias{path.lines}
\alias{path.cycles}
\alias{path.cycles.summary}

\title{Thermodynamic paths and cycles}

\description{Functions to calculate thermodynamic paths and cycles}

\usage{
cpcv.cal(datafile, plots = FALSE)
RefCoefAir
cp.cv(TC,ref=RefCoefAir)
fcp(TC)
fcpcv(TC)
fcv(TC)
fpv(V,p)
simpson(fun, a, b, pts = 100)
simpson.pv(nRT, a, b, pts = 100)
path.calc(x)
path.summary(y)
isochor(x)
isobar(x)
adiabat(x)
isotherm(x)
phase()
path.lines(x, plane = "Pv", shade.between = FALSE, lab.cycle = FALSE, shade.cycle = FALSE)
path.cycles(x, plane = "Pv", shade.cycle = FALSE)
path.cycles.summary(y)
}

\arguments{
	\item{datafile}{name of file for calibration}
	\item{plots}{Logical TRUE for plot}
	\item{TC}{Temp in deg C}
	\item{ref}{Ref and Coef argument to cp.cv; default is RefCoefAir}
	\item{V}{Volume}
	\item{p}{Pressure}
	\item{fun}{function to integrate}
	\item{a}{lower bound of interval}
	\item{b}{upper bound of interval}
	\item{pts}{number of points within the interval}
	\item{nRT}{gas law to calculate pV}
	\item{x}{argument specifying the states for any of isochor, isobar, adiabat, and isotherm ptah functions. These fucntions are used internally by path.calc. Argument x of path.calc specifies states for each one of isotherm, adiabat, isochor, isobar: V, P, n, T, M, S or to path.lines shade.under is a logical and part of the x argument specification; it is used by path.lines to shade area under the curve for a given path; default=FALSE or no shade. Argument x is also an argument to path.cycles as list of states in the cycle and cycle type; "carnot", "brayton","otto","diesel","stirling","box"}
	\item{y}{argument to path.summary or to path.cycles.summary: return from any path function as list(v,V,P,T,s,S,W,Q,cv,cp,gamma,WQtot,call,nM)}
	\item{shade.between}{logical used by path.lines to shade area between to paths; default=FALSE or no shade}
	\item{lab.cycle}{logical used by path.lines to label paths as a thermodynamic cycle; default=FALSE or no cycle}
	\item{shade.cycle}{logical used by path.lines and path.cycles to shade area within a thermodynamic cycle; default=FALSE or no shade}
	\item{plane}{Thermodynamic plane "Ts" or "Pv" to draw cycle.  Default is Pv}
}

\details{
Function cpcv.cal calibrates heat capacities as a function of temperature from a datafile.
Results for data file AirCvCpTK.csv (same as dataset AirCvCpTK) are stored in RefCoefAir for easy reference.
The result of cpcv.cal for other files can be used for the argument ref of cp.cv.
Function cp.cv calculates heat capacity at a given temperature.
Simpson integrates a function and simpson.pv integrates pv work using gas law.
Function path.calc arguments, V is an array with values of V1 and V2, P is array of corresponding pressure values,
T can be left unspecified and calculated from P and V. Argument path declares the type of process. argument lab are labels for the initial and final states of the path. Default values are one mole of dry air with M =28.97
Other calculation modes such as providing T instead of P, to calculate temperature from pressure and volume.
The call y<-path.calc(x) produces y that contains values of specific volume (m3/kg), volume (liters), pressure (bar), temperature (degC), specific heat W and work Q in kJ/kg, heat capacities in kJ/kgK, and gamma.
By default there are 1001 points calculated. 
Function path.summary(y) provides a view of the path.calc results.
Function path.lines(x)uses path.calc and helps visualize a paths x; or a list of paths, e.g., list(x1,x2) where x1 and x2 are paths.
These paths could be more than 2 and form a thermodynamic cycle.
Function path.cycles calculates and plots a thermodynamic cycle
Function path.cycles.summary(y) provides a view of the path.cycles results}

\value{
	\item{TK.ref}{T ref in K}
	\item{cv.ref}{cv ref}
	\item{cpcv.ref}{cp to cv ratio ref}
	\item{cv.coef}{cv coefficents}
	\item{cpcv.coef}{cp to cv ratio coefficents}
	\item{cv}{cv per mole}
	\item{cp}{cp per mole}
	\item{cp.cv}{cp to cv ratio}
	\item{cv.kg}{cv per kg}
	\item{cp.kg}{cp per kg}
	\item{s}{integration result by Simpson method}
	\item{y}{return list from path.calc or path.cycles. Used by path.summary and by path.lines as well as by path.cycles.summary}
	\item{start.end}{start and end of path}
	\item{WQtot}{Total work and heat}
	\item{pts}{mumber of points calculated}
	\item{call=y$call}{function called}
	\item{nM}{Number of moles}
}

\references{
Acevedo, M.F. 2018. Introduction to Renewable Electric Power Systems and the Environment with R. Boca Raton, FL: CRC Press. (ISBN 9781138197343)
}

\author{Miguel F. Acevedo \email{acevedo@unt.edu}}

\note{Functions used in Chapter 4 of Acevedo (2018)}

\examples{

# heat capacities at a given temperature, say 100C
cp.cv(TC=100)

# cv at the middle of a range of temperature 
cp.cv(TC=(323+25)/2)$cv

# arguments: V(l), P(bar), T(C), n(mol), M(g/mol)
# default n=1,M=28.9
# example specify V and P
x <- list(V=c(24.78,NA),P=c(1,2),path='isochor',lab=c("1","2"))
y <- path.calc(x)
path.summary(y)

# specify volume and calculate pressure
x <- list(V=c(24.78,34.78),P=c(1,NA),path='isobar',lab=c("3","4"))
y <- path.calc(x)
path.summary(y)

# example specify V and T visualize the path
x <- list(V=c(10,30),T=c(30,NA),path='isotherm',lab=c("1","2"))
path.lines(x)

# visualize two paths
x1 <- list(V=c(24.78,NA),P=c(1,2),path='isochor',lab=c("1","2"))
x2 <- list(V=c(24.78,34.78),P=c(1,NA),path='isobar',lab=c("3","4"))
path.lines(list(x1,x2))

# visualize two paths and shade in between curves
x1 <- list(V=c(10,30),T=c(200,200),path='isotherm',lab=c("3","4"))
x2 <- list(V=c(10,30),T=c(30,30),path='isotherm',lab=c("1","2"))
x <- list(x1,x2)
path.lines(x,shade.between=TRUE)

# carnot cycle
x <- list(TH=200,TL=30,V1=20,V4=40,cty='carnot')
y <- path.cycles(x,shade.cycle=TRUE)

}

\keyword{ hplot }
\keyword{ dplot }
\keyword{ math }

