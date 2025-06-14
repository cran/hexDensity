% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHexDensity.R
\name{plotHexDensity}
\alias{plotHexDensity}
\title{Plotting method for hexagonal Kernel Density Estimation}
\usage{
plotHexDensity(
  hexDensity,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xaxt = TRUE,
  yaxt = TRUE,
  lcex = 1,
  colramp = colorRampPalette(col.viridis),
  colorcut = 1024,
  legend = TRUE,
  legendWidth = 0.05,
  legendDistance = 0.15,
  aspectRatio = diff(hexDensity@xbnds)/diff(hexDensity@ybnds),
  margin = 0.18,
  newpage = TRUE
)
}
\arguments{
\item{hexDensity}{hexbin object returned by hexDensity}

\item{main}{Main title}

\item{xlab, ylab}{x-axis and y-axis label}

\item{xaxt, yaxt}{Logical. Whether to plot x,y axes}

\item{lcex}{Expansion factor for all letters.}

\item{colramp}{Color function that accept an integer n and return n colors.}

\item{colorcut}{An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1].}

\item{legend}{Legend is currently non-functional and should be ignored.}

\item{legendWidth}{Expansion factor for legend width.}

\item{legendDistance}{Expansion factor for the space between the plot and the legend.s}

\item{aspectRatio}{width to height ratio of the plot. Default is the (inverse of) shape value of hexDensity.}

\item{margin}{Minimum guaranteed margin for the plot. Different aspect ratio between the screen and the plot means that margin can be larger on certain sides.}

\item{newpage}{logical for whether to plot on a new page.}
}
\value{
No return value
}
\description{
Adapted the plotting function from \link[hexbin]{hexbin}. X and Y axes now
have the same scale with option for different aspect ratio. Ribbon legend 
for continuous data.
}
\section{SIDE EFFECTS}{
 Create kernel density estimate plot with hexagons
}

\examples{
set.seed(133)
d = hexDensity(x=rnorm(200),y=rnorm(200),bandwidth=0.15)
plotHexDensity(d)
}
\author{
Dan Carr <dcarr@voxel.galaxy.gmu.edu>; ported and extended by 
Nicholas Lewin-Koh nikko@hailmail.net. Modified by Quoc Hoang Nguyen 
<nguyen.q@wehi.edu.au> for hexDensity.
}
