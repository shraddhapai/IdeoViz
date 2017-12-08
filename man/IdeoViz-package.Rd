\name{IdeoViz-package}
\alias{IdeoViz-package}
\docType{package}
\title{Plots continuous or discrete data along chromosomal ideogram}
\description{Plotting discrete or continuous dataseries in the context of chromosomal location has several useful applications in genomic analysis. Examples of possible metrics include RNA expression levels, densities of epigenetic marks or genomic variation, while applications could range from the analysis of a single variable in a single context, to multiple measurements in several biological contexts (e.g. age/sex/tissue/disease context). Visualization of metrics against the chromosome could identify:
\enumerate{
\item{could identify distinctive spatial distribution that could further hypotheses about the functional role of the metric (e.g. telocentric or pericentromeric enrichment)}
\item{could highlight distribution differences between different groups of samples, suggesting different regulatory mechanisms; in extreme cases, visualization may identify large genomic foci of differences}
\item{could confirm that a quantitative difference measured between groups of interest is consistent throughout the genome (i.e. that there are no foci, and that the change is global).}
}
   This package provides a method to plot one or several dataseries against the chromosomal ideogram. It provides some simple options (vertical/horizontal orientation, display in bars or linegraphs). Data are expected to be binned; IdeoViz provides a function for user-specified bin widths. Ideograms for the genome of choice can also be automatically downloaded from UCSC using the getIdeo() function. 
}
\details{
\tabular{ll}{Package: \tab IdeoViz\cr
Type: \tab Package\cr
Title: \tab Plots data (continuous/discrete) along chromosomal ideogram\cr
Version: \tab 0.99.1\cr
Date: \tab 2013-06-26\cr
Author: \tab Shraddha Pai <Shraddha.Pai@camh.ca>, Jingliang Ren\cr
Maintainer: \tab Shraddha Pai <Shraddha.Pai@camh.ca>\cr
Depends: \tab Biobase, IRanges, GenomicRanges, RColorBrewer, rtracklayer\cr
biocViews: \tab Visualization,Microarray\cr
License: \tab GPL-2\cr}
}
\author{Shraddha Pai <Shraddha.Pai@camh.ca>, Jingliang Ren}

\keyword{ package }


