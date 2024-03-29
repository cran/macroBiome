# macroBiome 0.4.0

Released 2024-01-15

## Changed

-   for the functions with a grid mode, some function arguments have changed, now these functions can handle SpatRaster objects in addition to Raster objects
-   the 'Grid' functions have changed considering the class of the return value, now these functions return SpatRaster objects
-   under the hood, some remarkable changes have been made in the function `cliBrtSunDurFrcPoints()`, (e.g., the map used to classify continents and regions has been replaced).


# macroBiome 0.3.0

Released 2023-05-29

## Fixed

-   fixed some bugs to comply with changes in package dependencies


# macroBiome 0.2.0

Released 2023-02-12

## Added

-   function `cliKoppenPoints()` designate the Köppen-Geiger classification (KGC) type and calculate the associated bioclimiatic indicies using the monthly time series of temperature and precipitation
-   function `cliKoppenGrid()` to apply the above algorithm to the appropriate raster datasets
-   in the functions `cliBioCliIdxPoints()` and `cliBioCliIdxGrid()`, the range of selectable bioclimatic indices has been expanded by those used in the KGC system

## Changed

-   the classification schemes based on the Holdridge Life Zone (HLZ) system are now vectorized

## Fixed

-   function `cliHoldridgeGrid()` ignored the class "BaSl" (Bare soil and no vegetation) marked with a value of 39
