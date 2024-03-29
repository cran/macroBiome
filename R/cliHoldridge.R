#' Vegetation Classifier Using the HLZ System
#'
#' @description Calculates the values of bioclimatic indices used in the Holdridge life zone (HLZ) system (Holdridge
#'     1947, 1967), and designates the HLZ type using these values, by using the monthly time series of temperature
#'     and precipitation.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#'
#' @details To classify vegetation, the HLZ system developed by Holdridge (1947, 1967) uses the values of the
#'     following 3 bioclimatic indices:
#'
#'     \itemize{
#'       \item{\code{abt}: Mean Annual Biotemperature (Eq 1 in Szelepcsényi et al. (2014); in °C)}
#'       \item{\code{tap}: Total Annual Precipitation (in mm)}
#'       \item{\code{per}: Potential Evapotranspiration Ratio (Eq 4 in Szelepcsényi et al. (2014); dimensionless)}
#'     }
#'
#'     For details about calculating bioclimatic indices, see the function
#'     \code{\link[macroBiome]{cliBioCliIdxPoints}}. \cr
#'     The HLZ system classifies the vegetation type based on the distance from the ideal (theoretical) point in the
#'     3-dimensional space of bioclimatic indices. Numerous variants of the HLZ system are known (e.g.,
#'     Henderson-Sellers 1994; Yates et al. 2000). Here, one of its most widely used versions ('version with no
#'     altitudinal belts') is implemented, in accordance with works of Szelepcsényi et al. (2014, 2018). In this
#'     version, a total of 39 HLZ types are distinguished (see \code{\link[macroBiome]{vegClsNumCodes}}).
#'
#' @return Depending on the setting, a data frame with one or more columns where the HLZ types are stored in the last
#'     (character) column, while the additional columns contain the values of bioclimatic indices used. The
#'     abbreviations of HLZ types can be found in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a one-column data frame with the HLZ types.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'     follows: \code{'temp'} (one-year time series of monthly mean air temperature), and \code{'prec'} (one-year
#'     time series of monthly precipitation sum). The objects \code{'temp'} and \code{'pre'} must be either vectors
#'     of length 12 or 12-column matrices. The first dimensions of these matrices have to be the same length. The
#'     function automatically converts vectors into single-row matrices during the error handling, and then uses
#'     these matrices. The first dimensions of these matrices determines the number of rows in the result matrix.
#'
#' @references
#'
#' \cite{Henderson-Sellers A (1994) Global terrestrial vegetation ‘prediction’: the use and abuse of climate and
#'     application models. Prog Phys Geogr 18(2):209–246. \doi{10.1177/030913339401800203}}
#'
#' \cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' \cite{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
#'     original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
#'     \doi{10.2478/s13533-012-0189-5}}
#'
#' \cite{Szelepcsényi Z, Breuer H, Kis A, Pongrácz R, Sümegi P (2018) Assessment of projected climate change in the
#'     Carpathian Region using the Holdridge life zone system. Theor Appl Climatol 131(1–2):593–610.
#'     \doi{10.1007/s00704-016-1987-3}}
#'
#' \cite{Yates DN, Kittel TGF, Cannon RF (2000) Comparing the Correlative Holdridge Model to Mechanistic
#'     Biogeographical Models for Assessing Vegetation Distribution Response to Climatic Change. Clim Chang
#'     44(1–2):59–87. \doi{10.1023/A:1005495908758}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#' data(vegClsNumCodes)
#'
#' # Designate the HLZ type (using the related bioclimatic indices),
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' HLZ <- cliHoldridgePoints(colMeans(temp), colMeans(prec), verbose = TRUE)
#' numCode <- which(sapply(vegClsNumCodes$Code.HLZ, identical, HLZ[, "vegCls"]))
#' cbind(HLZ[,-c(4)], vegClsNumCodes[numCode, c("Name.HLZ", "Code.HLZ")])
#' })
#'
#' @importFrom stats complete.cases setNames
#'
#' @export
#'
cliHoldridgePoints <- function(temp, prec, verbose = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(temp = temp, prec = prec)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.bci <- c("abt", "tap", "per")

  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of each bioclimatic index required to classify vegetation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bioCliIdx <- cliBioCliIdxPoints(temp, prec, bciOpVar = cv.bci, argCkd = T)
  list2env(unclass(as.data.frame(bioCliIdx)), envir = environment())

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object containing vegetation class
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vegCls <- rep(NA_character_, length = lgth)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Determine the vegetation class by using the Holdridge life zone system
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set a board containing the distance from the optimal value for each classes
  distBds <- matrix(nrow = lgth, ncol = nrow(hlzDefSubset), dimnames = list(NULL, hlzDefSubset$ID))
  for (i_vcl in 1 : nrow(hlzDefSubset)) {
    optVabt <- log2(hlzDefSubset$abt[i_vcl]) + 0.5
    tmpVabt <- ifelse(abt == 0., NA, (log2(abt) - optVabt) ** 2.)
    optVtap <- log2(hlzDefSubset$tap[i_vcl]) + 0.5
    tmpVtap <- ifelse(tap == 0., NA, (log2(tap) - optVtap) ** 2.)
    optVper <- log2(hlzDefSubset$per[i_vcl]) + 0.5
    tmpVper <- ifelse(per == 0., NA, (log2(per) - optVper) ** 2.)
    distBds[, i_vcl] <- sqrt(tmpVabt + tmpVtap + tmpVper)
  }

  # Calculate minimum distances in a three-dimensional space of bioclimatic indices
  minVal <- do.call(pmin, as.data.frame(distBds))

  # Select all possible classes using minimum distances
  psblVegCls <- t(sapply(1 : nrow(distBds), function(i) { distBds[i, ] == minVal[i] }))

  # Number of possible classes
  n_pvc <- rowSums(psblVegCls)

  # Set some magic numbers
  # Polar temperature line
  plTL <- 1.5

  # Minimum value of the total annual precipitation
  mn_tap <- 62.5

  # Minimum value of the potential evapotranspiration ratio
  mn_per <- 0.125

  # Frost or critical temperature line
  frTL <- 2. ** (log2(12.) + 0.5)

  # Vegetation classes in warm temperate and subtropical regions
  WtCls <- seq(17, 23)
  names(WtCls) <- c("WtD", "WtDs", "WtTs", "WtDf", "WtMf", "WtWf", "WtRf")
  StCls <- seq(24, 30)
  names(StCls) <- c("StD", "StDs", "StTw", "StDf", "StMf", "StWf", "StRf")

  # Classify the vegetation type under certain boundary conditions
  gr        <- !is.na(n_pvc)
  grI       <- gr & (tap < mn_tap | per < mn_per)
  grII      <- gr & (tap >= mn_tap & per >= mn_per)
  grIIA     <- grII & abt < plTL
  grIIB     <- grII & abt >= plTL
  grIIB1    <- grIIB & n_pvc == 1
  grIIB2    <- grIIB & n_pvc != 1
  grIIB2a   <- grIIB2 & apply(psblVegCls[, c(WtCls, StCls), drop = FALSE], 1, function(x) !any(x))
  grIIB2b   <- grIIB2 & apply(psblVegCls[, c(WtCls, StCls), drop = FALSE], 1, any)
  grIIB2bwt <- grIIB2b & abt < frTL
  grIIB2bst <- grIIB2b & abt >= frTL

  selVegCls <- function(x) { names(which(x))[1] }

  # I. Bare soil and no vegetation
  vegCls[grI]       <- "BaSl"

  # II.A Polar desert
  vegCls[grIIA]     <- "PD"

  # II.B.1 One possible valid vegetation type
  vegCls[grIIB1]    <- apply(psblVegCls[grIIB1, , drop = FALSE], 1, function(x) names(which(x)))

  # II.B.2.a Two or more possible valid vegetation types (without warm temperate and subtropical regions)
  vegCls[grIIB2a]   <- apply(psblVegCls[grIIB2a, , drop = FALSE], 1, selVegCls)

  # II.B.2.bwt Two or more possible valid vegetation types (with warm temperate region)
  vegCls[grIIB2bwt] <- apply(psblVegCls[grIIB2bwt, -c(StCls), drop = FALSE], 1, selVegCls)

  # II.B.2.bst Two or more possible valid vegetation types (with subtropical region)
  vegCls[grIIB2bst] <- apply(psblVegCls[grIIB2bst, -c(WtCls), drop = FALSE], 1, selVegCls)

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (verbose) {
    rslt <- data.frame(bioCliIdx, vegCls = vegCls)
  } else {
    rslt <- data.frame(vegCls = vegCls)
  }
  return(rslt)

}


#' Vegetation Classifier Using the HLZ System
#'
#' @description Calculates the values of bioclimatic indices used in the Holdridge life zone (HLZ) system
#'     (Holdridge 1947, 1967), and designates the HLZ type using these values, for a given region, by using the
#'     monthly time series of temperature and precipitation.
#'
#' @param rs.temp multi-layer Raster*/SpatRaster object with one-year time series of monthly mean air temperature
#'     (in °C)
#' @param rs.prec multi-layer Raster*/SpatRaster object with one-year time series of monthly precipitation sum
#'     (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[terra]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliHoldridgePoints}}.
#'
#' @return Depending on the setting, a SpatRaster object with one or more layers where the numeric integers encoding
#'     the HLZ type are stored at the last layer, while the additional layers contain the values of bioclimatic
#'     indices used. The meaning of integers is given in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a single-layer SpatRaster obejct with numeric integers encoding
#'     the HLZ type.
#'
#' @note The objects \code{'rs.temp'} and \code{'rs.prec'} must be 12-layer Raster*/SpatRaster objects. These
#'     Raster*/SpatRaster objects must have the same bounding box, projection, and resolution.
#'
#' @references
#'
#' \cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' @examples
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#'
#' # Designate the HLZ types (using the related bioclimatic indices)
#' # for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' rs.HLZ <- cliHoldridgeGrid(temp, prec, verbose = TRUE)
#' rs.HLZ
#' })
#'
#' @importFrom methods as
#' @import terra
#'
#' @export
#'
cliHoldridgeGrid <- function(rs.temp, rs.prec, verbose = FALSE, filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  cv.arg <- c("rs.temp", "rs.prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  err_han <- errorHandlingGrid(rs.temp = rs.temp, rs.prec = rs.prec)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  rs.aux <- terra::subset(rs.temp, 1)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  n_lyr <- ifelse(verbose, 4, 1)

  rs.rslt <- terra::rast(rs.aux, nlyrs = n_lyr)

  cv.mly_var <- c("rs.temp", "rs.prec")
  cv.arg <- cv.mly_var
  for (i_arg in 1 : length(cv.arg)) {
    if (!is.null(get(cv.arg[i_arg]))) {
      x <- get(cv.arg[i_arg])
      terra::readStop(get(cv.arg[i_arg]))
      if (!terra::readStart(get(cv.arg[i_arg]))) { stop(x@ptr$messages$getError()) }
      on.exit(terra::readStop(get(cv.arg[i_arg])))
      rm(x)
    }
  }

  overwrite <- list(...)$overwrite
  if (is.null(overwrite)) overwrite <- FALSE
  wopt <- list(...)$wopt
  if (is.null(wopt)) wopt <- list()

  b <- terra::writeStart(rs.rslt, filename, overwrite, wopt = wopt)

  for (i in 1 : b$n) {
    for (i_arg in 1 : length(cv.arg)) {
      if (!is.null(get(cv.arg[i_arg]))) {
        x <- get(cv.arg[i_arg])
        assign(substring(cv.arg[i_arg], 4), terra::readValues(x, row = b$row[i], nrows = b$nrows[i], col = 1,
                                                              ncols = ncol(x), mat = TRUE))
        rm(x)
      } else {
        assign(substring(cv.arg[i_arg], 4), NULL)
      }
    }
    df.rslt <- cliHoldridgePoints(temp, prec, verbose = verbose)
    numCode <- hlzDefSubset$Numeric.code[match(df.rslt[["vegCls"]], hlzDefSubset$ID)]
    df.rslt[["vegCls"]] <- numCode

    terra::writeValues(rs.rslt, as.matrix(df.rslt), b$row[i], b$nrows[i])
  }
  terra::writeStop(rs.rslt)

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  names(rs.rslt) <- colnames(df.rslt)
  return(rs.rslt)

}
