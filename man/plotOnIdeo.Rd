\name{plotOnIdeo}
\alias{plotOnIdeo}
\title{Plot data superimposed on chromosomal ideogram}
\description{Main function to plot binned data alongside chromosomal ideogram.}
\usage{plotOnIdeo(chrom = stop("enter chromosome(s) to plot"), ideoTable, 
    values_GR, value_cols = "values", plotType = "lines", col = "orange", 
    bpLim = NULL, val_range = NULL, addScale = TRUE, scaleChrom = TRUE, 
    vertical = FALSE, addOnetoStart = TRUE, smoothVals = FALSE, cex.axis = 1,     plot_title = NULL, ablines_y = NULL, cex.main=1, ...)}
\arguments{
  \item{chrom}{(character) chromosome(s) to create ideograms for}
  \item{ideoTable}{(data.frame) ideogram table. See \code{getIdeo()}}
  \item{values_GR}{(GenomicRanges) data to be plotted must be in metadata columns}
  \item{value_cols}{(character) which series to plot. Should be column names of the \code{mcols()} slot of \code{values_GR}}
  \item{plotType}{(character) Plot type for each series. Values can be "\code{lines}" or "\code{rect}" to plot lines or barplots respectively. The latter is not recommended when several series are to be plotted on the same axis.)}
  \item{col}{(character) vector of colors for data series}
  \item{bpLim}{(numeric) (xlim); display only a section of the chromosome and the corresponding values}
  \item{val_range}{(numeric) (ylim); y-axis scale for data series}
  \item{addScale}{(logical) if TRUE, bp positions will be shown along the chromosomes. This feature should be turned off if numerous chromosomes' worth of data are being plotted and all objects don't fit on the final graphics device.}
  \item{scaleChrom}{(logical) if FALSE, all chroms will display as the same size. scaleChrom will be ignored if bpLim is not NULL}
  \item{vertical}{(logical) if TRUE, chromosomes will be plotted vertically}
  \item{addOnetoStart}{(logical) if TRUE, adds 1 to chromStart. Useful to convert data in half-open coordinates - which is all data from the UCSC genome browser, including cytoBandIdeo, into 1-base.}
  \item{smoothVals}{(logical) if T, smoothes each trendline. Currently hard-coded to lowess smoothing with span=0.03}
  \item{cex.axis}{(integer) axis font size}
  \item{plot_title}{(character) title for overall graph}
  \item{ablines_y}{(numeric) when supplied, draws reference lines on the y-axis}
  \item{cex.main}{(numeric) font size for plot title}
  \item{\dots}{other graphing options for barplot (i.e. main="Values", to title bar plot "Values")}
}
\details{
\code{plotOnIdeo()} is the main function of this package. It is the one the end-user is expected to call to generate plots. 
Input is provided as a \emph{GRanges} object (\code{values_GR}), with data to be plotted contained in its metadata slot. The user is responsible for providing pre-binned data, if binning is required. Data can also be binned using the avgByBin() function in this package. 
The ideogram table (\code{ideoTable}) is the same as the \emph{cytoBandIdeo} table available from the UCSC genome browser database for a given genome is a  can be either automatically downloaded from UCSC (see \code{getIdeo()}) or read in from a local-file and passed to this function. 

There are numerous arguments which control the appearance of the plot. The main decision points are:
  \enumerate{
		\item{\code{vertical}: Whether the entire plot should have a horizontal or vertical orientation}
		\item{\code{plotType}: Whether each dataseries should be shown as a trendline or as a barplot. The latter is not recommended for cases where multiple series are to be plotted on the same axis.}
	}
Other considerations:
\itemize{
  \item{The size of the graphics device limits the number of chromosomes that can be plotted. A simple solution may be to set \code{addScale=FALSE}. However, it is recommended to call \code{plotOnIdeo()} multiple times, and plotting a fewer number of chromosomes on each page. }
  \item{The code expects coordinates of \code{values_GR} to be in 1-base. Set \code{addOneToStart=TRUE} if supplied coordinates are in 0-base.}
}
  }
\examples{
data(binned_multiSeries)
data(hg18_ideo)
plotOnIdeo(chrom=seqlevels(binned_multiSeries),ideoTable=hg18_ideo, values_GR=binned_multiSeries, value_cols=colnames(mcols(binned_multiSeries)), col=1:5)

}



\author{Shraddha Pai <Shraddha.Pai@camh.ca>, Jingliang Ren}




