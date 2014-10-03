\name{plotChromValuePair}
\alias{plotChromValuePair}
\title{Plot a chromosome-value pair}
\description{Base function which plots the ideogram and superimposed data for a single chromosome. \code{plotOnIdeo()} calls this function and stacks the resulting output.}
\usage{plotChromValuePair(chrom, cytoTable, bpLim, vertical, values_GR, 
    val_range, col, value_cols = "values", default_margins, addScale, 
    ablines_y, smoothVals, span=0.03, verbose = FALSE, ...)}
\arguments{
  \item{chrom}{(character) chromosome(s) to create ideograms for }
  \item{cytoTable}{(data.frame) loaded ideogram table. (see \emph{ideoTable} argument to \code{plotOnIdeo()})}
  \item{bpLim}{(numeric) (aka xlim); display only a section of the chromosome and the corresponding values}
  \item{vertical}{(logical) if TRUE, chromosomes will be plotted vertically}
  \item{values_GR}{(GenomicRanges) data to be plotted must be in metadata columns}
  \item{val_range}{(numeric) (aka ylim); y-axis scale for data series}
  \item{col}{(character) colour for series}
  \item{value_cols}{(character) column name for series to plot}
  \item{default_margins}{(numeric) page inner margins (in inches)}
  \item{addScale}{(logical) if FALSE, bp positions will be hidden}
  \item{ablines_y}{(numeric) when specified, will draw reference lines on the y-axis}
  \item{smoothVals}{(logical) when T applies loess() to each series}
  \item{span}{(numeric)\code{span} argument for \code{loess} function}
  \item{verbose}{(logical) print messages}
  \item{\dots}{arguments to \code{axis()}, \code{line()}, and \code{rect()}}
}
\details{Plots one unit of chromosome ideogram with dataseries superimposed. Usually, the user can avoid this function and directly call \code{plotOnIdeo()}. However, this function may be used in cases where further plot customization is required.}
\examples{
data(hg18_ideo)
data(binned_multiSeries)
layout(matrix(1:2, byrow=TRUE,ncol=1),heights=c(2.5,1))
plotChromValuePair("chr1",hg18_ideo, 
	values_GR=binned_multiSeries,value_cols=colnames(mcols(binned_multiSeries)),plotType='lines',
	col=1:5,val_range=c(0,10),bpLim=NULL,vertical=FALSE,addScale=TRUE,ablines_y=NULL,
	smoothVals=FALSE,default_margins=c(0.5,.5,.1,.1))
}
\seealso{\code{plotOnIdeo()}}
\author{Shraddha Pai <Shraddha.Pai@camh.ca>, Jingliang Ren}





