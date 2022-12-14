% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlyWeaGen.R
\name{dlyWeaGenPoints}
\alias{dlyWeaGenPoints}
\title{Daily Weather Generator}
\usage{
dlyWeaGenPoints(
  temp,
  prec,
  bsdf,
  year = 2000,
  aprchTEMP = c("hip", "tsi", "const"),
  aprchPREC = c("tsi", "hip", "const"),
  aprchBSDF = c("hip", "const"),
  dvTEMP = rep(0.7, 12),
  dvPREC = rep(0.7, 12),
  argCkd = FALSE
)
}
\arguments{
\item{temp}{'numeric' R object with one-year time series of monthly mean air temperature (in °C)}

\item{prec}{'numeric' R object with one-year time series of monthly precipitation sum (in mm)}

\item{bsdf}{'numeric' R object with one-year time series of monthly mean relative sunshine duration (dimensionless)}

\item{year}{'numeric' vector with values of the year (using astronomical year numbering)}

\item{aprchTEMP}{'character' vector of length 1 that indicates the scheme used to generate daily values of the
daily mean air temperature for a specific year. Valid values are as follows: \cr
(a) \code{'hip'} -
this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
monthly mean air temperature in order to generate daily values; \cr
(b) \code{'tsi'} -
this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
air temperature, in order to generate a synthetic time series of the selected meteorological variable at a
temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled to a
daily resolution; \cr
(c) \code{'const'} -
this scheme is assumed that values of the daily mean air temperature are constant within each month.}

\item{aprchPREC}{'character' vector of length 1 that indicates the scheme to generate daily values of the
daily precipitation sum. Valid values are as follows: \cr
(a) \code{'tsi'} -
this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
precipitation intensity, in order to generate a synthetic time series of the selected meteorological variable
at a temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled
to a daily resolution; \cr
(b) \code{'hip'} -
this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
monthly mean precipitation intensity in order to generate daily values; \cr
(c) \code{'const'} -
this scheme is assumed that values of the daily precipitation sum are constant within each month (the monthly
precipitation sum is divided equally across each day of the month).}

\item{aprchBSDF}{'character' vector of length 1 that indicates the scheme used to generate daily values of the
daily fractional sunshine duration for a specific year. Valid values are as follows: \cr
(a) \code{'hip'} -
this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
monthly mean relative sunshine duration in order to generate daily values; \cr
(b) \code{'const'} -
this scheme is assumed that values of the daily relative sunshine duration are constant within each month.}

\item{dvTEMP}{'numeric' vector of length 12 with monthly values of the damping variable for the air temperature
data.}

\item{dvPREC}{'numeric' vector of length 12 with monthly values of the damping variable for the precipitation data.}

\item{argCkd}{'logical' scalar that indicates whether or not the checking and correction of arguments can be
omitted.}
}
\value{
A list with three 365- or 366-column matrices that contain quasi-daily time series for the three basic
    climate variables:

    \itemize{
      \item{\code{TEMP}: daily mean air temperature (in °C)}
      \item{\code{PREC}: daily precipitation sum (in mm)}
      \item{\code{BSDF}: daily fractional sunshine duration (dimensionless)}
    }
}
\description{
Generates quasi-daily time series from the monthly mean values of temperature, precipitation and
    sunshine data.
}
\details{
For relative sunshine duration and air temperature, it is recommended the 'harmonic' interpolation
    technique (\code{'hip'}) described by Epstein (1991), with a correction of physically impossible values. This
    technique can also be used to values of the mean precipitation intensity, however, in the case of
    precipitation, the temporal scaling (\code{'tsi'}) using an iterative interpolation technique described by
    Lüdeke et al. (1994) is recommended, with a damping variable of 0.7 for each month. The damping variable can
    be set separately for each month and must be between 0 and 1. This iterative scheme can also be applied to
    monthly mean data of air temperature. Furthermore, for all three climate variables, it is possible to ignore
    intra-month variability. For this, the setting \code{'const'} must be used.
}
\note{
As with any point function, a set of basic input data is defined here. In this case, they are as follows:
   \code{'temp'} (one-year time series of monthly mean air temperature), \code{'prec'} (one-year time series
   of monthly precipitation sum), and \code{'bsdf'} (one-year time series of monthly mean relative sunshine
   duration.) The objects \code{'temp'}, \code{'prec'} and \code{'bsdf'} must be either vectors of length 12 or
   12-column matrices. The first dimensions of these matrices have to be the same length. The function
   automatically converts vectors into single-row matrices during the error handling, and then uses these
   matrices. The first dimensions of these matrices determines the number of rows in the result matrix. In the
   case of arguments that do not affect the course of the calculation procedure or the structure of the return
   object, scalar values (i.e., 'numeric' vector of length 1) may also be allowed. In this case, it is as
   follow: \code{'year'} (year using astronomical year numbering). This scalar is converted to a vector by the
   function during the error handling, and this vector is applied in the further calculations. If these data are
   stored in vectors of length at least 2, their length must be the same size of first dimension of the matrices
   containing the basic data.
}
\examples{
library(graphics)

# Loading mandatory data for the Example 'Points'
data(inp_exPoints)

with(inp_exPoints, {
# Generate quasi-daily time series for basic climate variables with default settings,
# at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
year <- trunc(mean(seq(1981, 2010)))
wea01 <- dlyWeaGenPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), year = year)

# Modify the daily weather data generation techniques
# To temperature data, apply the iterative interpolation technique with basic settings
# To precipitation data, change the value of the damping variable, over the whole year
# To sunshine data, assume that its values are constant within each month
wea02 <- dlyWeaGenPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), aprchTEMP = "tsi",
    aprchBSDF = "const", dvPREC = rep(0.6, 12), year = year)

# Check the differences
vars <- c("TEMP", "PREC", "BSDF")
lbls <- list(expression(italic(T[a])~(~degree*C)), expression(italic(P[n])~(mm)),
    expression(italic(S[f])~(unitless)))
ys <- c(20, 2.5, 0.6)
ats <- list(seq(-4, 24, 4), seq(0, 3, 0.5), seq(0., 0.8, 0.2))
cols <- c("black", "green")
opar <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
for (i in 1 : length(vars)) {
  par(mar = c(2, 5, 1, 1))
  matplot(t(rbind(wea01[[vars[i]]], wea02[[vars[i]]])), type = "l", lwd = 2, col = cols,
      xaxt = "n", xlab = NA, ylab = NA, axes = FALSE)
  axis(side = 1, las = 1, tck = -0.03, labels = NA, at = seq(-60, 720, 30))
  axis(side = 2, las = 1, tck = -0.03, labels = NA, at = ats[[i]])
  axis(side = 2, las = 1, lwd = 0, line = -0.4, cex.axis = 1.6, at = ats[[i]])
  if (i == length(vars)) {
    axis(side = 1, las = 1, lwd = 0, line = -0.4, at = seq(-60, 720, 30), cex.axis = 1.6)
  }
  mtext(side = 2, lbls[[i]], line = 3, cex = 1.1)
  legend(1, ys[i], legend = c("default", "modified"), col = cols, lty = 1 : 2, lwd = 2, xpd = TRUE)
}
par(opar)
})

}
\references{
\cite{Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
    \doi{10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2}}

\cite{Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
    Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented model of
    seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I. Model description
    and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166. \doi{10.3354/cr004143}}
}
