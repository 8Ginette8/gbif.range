% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluateRM.R
\name{evaluateRM}
\alias{evaluateRM}
\title{Evaluates the sensitivity & precision of range maps based on validation
data, such as predictions of species distributions (SDMs) or IUCN expert
range maps.}
\usage{
evaluateRM(
  root.dir = NULL,
  valData.dir = NULL,
  ecoRM.dir = NULL,
  valData.type = NULL,
  verbose = TRUE,
  print.map = TRUE,
  mask = NULL,
  res.fact = NULL
)
}
\arguments{
\item{root.dir}{Character. Root directory to files}

\item{valData.dir}{Numeric. Buffer width parameter}

\item{ecoRM.dir}{Numeric. Number of observation points}

\item{valData.type}{Character. Type of valData - either "SHP" or "TIFF"}

\item{verbose}{Logical. Optional - report details while running}

\item{print.map}{Logical. Optional - if verbose=TRUE should a overlap map be printed}

\item{mask}{?}

\item{res.fact}{Integer. Factor for coarsening the original resolution}
}
\description{
Evaluates the sensitivity & precision of range maps based on validation
data, such as predictions of species distributions (SDMs) or IUCN expert
range maps.
}
\references{
Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential of 
ecoregional range maps for boosting taxonomic coverage in ecology and 
conservation. Ecography, 12, e06794.
}
