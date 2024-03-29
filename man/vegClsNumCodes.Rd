% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vegClsNumCodes}
\alias{vegClsNumCodes}
\title{Supplemental Data Frame for Decoding Outputs of Vegetation Classifiers}
\format{
A data frame that allows for interpreting the return objects provided by climate-based vegetation
    classifiers implemented here. Two columns belong to each vegetation classification approach. Columns whose
    names begin with the string \code{'Name.'} contain the full names of the vegetation classes. While columns
    whose names begin with the string \code{'Code.'} summarize the abbreviations used by functions with a point
    mode. Row numbers of the data frame have a special role because they are the same as the numbers returned by
    the functions with a grid mode.
}
\description{
The key to the classes used by climate-based vegetation classifiers implemented here. Currently, three
    bioclimatic vegetation classification approaches are implemented:
    \itemize{
      \item{\code{HLZ}: a version with no altitudinal belts of the Holdridge life zone (HLZ) system (Holdridge
      1947, 1967), in accordance with works of Szelepcsényi et al. (2014, 2018)}
      \item{\code{KGC}: the Köppen-Geiger classification (KGC) system (Köppen 1936) with some modifications
      suggested by Peel et al. (2007)}
      \item{\code{BIOME}: the initial version of the BIOME model developed by Prentice et al. (1992)}
    }
}
\examples{
data(vegClsNumCodes)
str(vegClsNumCodes)

}
\references{
\emph{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
    105(2727):367–368. \doi{10.1126/science.105.2727.367}}

\emph{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}

\cite{Köppen W (1936) Das geographische System der Klimate. In: Köppen W, Geiger R (eds) Handbuch der
    Klimatologie. Verlag von Gebrüder Borntraeger, Berlin, Germany, pp 1–44}

\cite{Peel MC, Finlayson BL, McMahon TA (2007) Updated world map of the Köppen-Geiger climate classification.
    Hydrol Earth Syst Sci 11(5):1633–1644. \doi{10.5194/hess-11-1633-2007}}

\emph{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A global biome model based on
    plant physiology and dominance, soil properties and climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}

\emph{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
    original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
    \doi{10.2478/s13533-012-0189-5}}

\emph{Szelepcsényi Z, Breuer H, Kis A, Pongrácz R, Sümegi P (2018) Assessment of projected climate change in the
    Carpathian Region using the Holdridge life zone system. Theor Appl Climatol 131(1–2):593–610.
    \doi{10.1007/s00704-016-1987-3}}
}
