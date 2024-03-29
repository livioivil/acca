% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccbiplot.R
\name{ccbiplot}
\alias{ccbiplot}
\title{ccbiplot}
\usage{
ccbiplot(
  pcobj,
  choices = 1:2,
  scale = 1,
  pc.biplot = TRUE,
  obs.scale = 1 - scale,
  var.scale = scale,
  groups = NULL,
  shapes = NULL,
  ellipse = FALSE,
  ellipse.prob = 0.68,
  labels = NULL,
  labels.size = 3,
  alpha = 1,
  var.axes = TRUE,
  circle = FALSE,
  circle.prob = 0.69,
  varname.size = 3,
  varname.adjust = 1.5,
  varname.abbrev = FALSE,
  arrows.color = NULL,
  ...
)
}
\arguments{
\item{pcobj}{an object returned by acca:cc()}

\item{choices}{which PCs to plot}

\item{scale}{covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.}

\item{pc.biplot}{for compatibility with biplot.princomp()}

\item{obs.scale}{scale factor to apply to observations}

\item{var.scale}{scale factor to apply to variables}

\item{groups}{optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups}

\item{shapes}{which shapes to use when plotting observations}

\item{ellipse}{draw a normal data ellipse for each group?}

\item{ellipse.prob}{size of the ellipse in Normal probability}

\item{labels}{optional vector of labels for the observations}

\item{labels.size}{size of the text used for the labels}

\item{alpha}{alpha transparency value for the points (0 = transparent, 1 = opaque)}

\item{var.axes}{draw arrows for the variables?}

\item{circle}{draw a correlation circle?}

\item{circle.prob}{size of the circle in qchisq(circle.prob, df = 2)}

\item{varname.size}{size of the text for variable names}

\item{varname.adjust}{adjustment factor the placement of the variable names, >= 1 means farther from the arrow}

\item{varname.abbrev}{whether or not to abbreviate the variable names}

\item{arrows.color}{color of the arrows}

\item{...}{not used yet}
}
\value{
It returns a \code{\link[ggplot2]{ggplot}} object.
}
\description{
The function adapts the \code{ggbiplot::ggbiplot()} function to an \code{acca}-object.
}
