#' Plot data superimposed on chromosomal ideogram
#' @usage
#' plotOnIdeo(chrom = stop("enter chromosome(s) to plot"), ideoTable, 
#'     values_GR, value_cols = "values", plotType = "lines", 
#' 	 col = "orange", 
#'     bpLim = NULL, val_range = NULL, addScale = TRUE, scaleChrom = TRUE, 
#'     vertical = FALSE, addOnetoStart = TRUE, smoothVals = FALSE, 
#' 	cex.axis = 1,plot_title = NULL, ablines_y = NULL, cex.main=1, ...)
#' 
#' @param chrom (character) chromosome(s) to create ideograms for
#' @param ideoTable (data.frame) ideogram table. See getIdeo()
#' @param values_GR (GenomicRanges) data to be plotted must be in metadata 
#' columns
#' @param value_cols (character) which series to plot. Should be column names 
#' of the mcols() slot of values_GR
#' @param plotType (character) Plot type for each series. Values can be "lines"
#'  or "rect" to plot lines or barplots respectively. The latter is not recommended when several series are to be plotted on the same axis.)
#' @param col (character) vector of colors for data series
#' @param bpLim (numeric) (xlim); display only a section of the chromosome and the corresponding values
#' @param val_range (numeric) (ylim); y-axis scale for data series
#' @param addScale (logical) if TRUE, bp positions will be shown along the chromosomes. This feature should be turned off if numerous chromosomes' worth of data are being plotted and all objects don't fit on the final graphics device.
#' @param scaleChrom (logical) if FALSE, all chroms will display as the same size. scaleChrom will be ignored if bpLim is not NULL
#' @param vertical (logical) if TRUE, chromosomes will be plotted vertically
#' @param addOnetoStart (logical) if TRUE, adds 1 to chromStart. Useful to convert data in half-open coordinates - which is all data from the UCSC genome browser, including cytoBandIdeo, into 1-base.
#' @param smoothVals (logical) if T, smoothes each trendline. Currently hard-coded to lowess smoothing with span=0.03
#' @param cex.axis (integer) axis font size
#' @param plot_title (character) title for overall graph
#' @param ablines_y (numeric) when supplied, draws reference lines on the y-axis
#' @param cex.main (numeric) font size for plot title. 
#' @param \dots other graphing options for barplot (i.e. main="Values", to title bar plot "Values")
#' @details Main function to plot binned data alongside chromosomal ideogram.
#' plotOnIdeo() is the main function of this package. It is the one the end-user is expected to call to generate plots. 
#' Input is provided as a GRanges object (values_GR), with data to be plotted contained in its metadata slot. The user is responsible for providing pre-binned data, if binning is required. Data can also be binned using the avgByBin() function in this package. 
#' The ideogram table (ideoTable) is the same as the cytoBandIdeo table available from the UCSC genome browser database for a given genome is a  can be either automatically downloaded from UCSC (see getIdeo()) or read in from a local-file and passed to this function. 
#' 
#' There are numerous arguments which control the appearance of the plot. 
#' The main decision points are:
#'   \enumerate{
#' 		\item{vertical: Whether the entire plot should have a horizontal or vertical orientation}
#' 		\item{plotType: One of [rect|lines|seg_tracks]. 
#' Type of plot, trendline ("lines"), barplot ("rect") or tracks of 
#' GenomicRanges (seg_tracks). "rect" only works when there is a single
#' data series (single set of values) to be plotted on the same axis. 
#' }}
#' 	
#' Other considerations:
#' \itemize{
#'   \item{The size of the graphics device limits the number of chromosomes that can be plotted. A simple solution may be to set addScale=FALSE. However, it is recommended to call plotOnIdeo() multiple times, and plotting a fewer number of chromosomes on each page.}
#'   \item{The code expects coordinates of values_GR to be in 1-base. Set addOneToStart=TRUE if supplied coordinates are in 0-base.}
#' }
#' 
#'   
#' @examples
#' data(binned_multiSeries)
#' data(hg18_ideo)
#' plotOnIdeo(chrom=seqlevels(binned_multiSeries),
#'	ideoTable=hg18_ideo, values_GR=binned_multiSeries, 
#'	value_cols=colnames(mcols(binned_multiSeries)), col=1:5)
#' 
#' @export
plotOnIdeo <- function(chrom=stop("enter chromosome(s) to plot"),
ideoTable,values_GR,value_cols='values',plotType='lines',col='orange', 
bpLim=NULL,val_range=NULL, addScale=TRUE,scaleChrom=TRUE,vertical=FALSE,
addOnetoStart=TRUE,smoothVals=FALSE,cex.axis=1.0,plot_title=NULL,
ablines_y=NULL,cex.main=1,... ) 
{ 
    def_las <- par("las")
    def_cexaxis <- par("cex.axis")
    def_fontaxis <- par("font.axis")
    if (vertical) {
        default_margins <- c(0,0.01,1,0)
        omi <- c(0.3,0,1,0)
    } else {
        default_margins <- c(0.01,0.4,0,0)
        omi <- c(0.1,0.3,0.3,0.3)
    }
    
    if(!all(chrom %in% ideoTable$chrom)) 
        stop("Not all chroms specified found in ideoTable.")
    ideoTable <- ideoTable[ideoTable$chrom %in% chrom, ]
    
    if (addOnetoStart) ideoTable$chromStart <- ideoTable$chromStart+1
    if (is.null(bpLim)){
        if (scaleChrom) { 
         bpLim <- c(min(ideoTable$chromStart),
					max(ideoTable$chromEnd)); 
        }
    }
	#cat("this point\n")

    # don't convert large numbers to scientific notation
    options(scipen=10) 

    #layout settings 
    par(las=1, cex.axis=cex.axis, font.axis=2,omi=omi,bty="n")
    panel_id <-1:(2*length(chrom))
    if(vertical) {
        layout(matrix(panel_id, byrow=FALSE,ncol=length(panel_id)), 
        widths=rep(c(1,2.5),length(chrom)) )
    } else {
        layout(matrix(panel_id, byrow=TRUE,ncol=1),heights=rep(c(2.5,1),
        length(chrom)))
    #par(mar=c(0,0,0,0))
    }
    	
    #title
    ctr <- 1
    suppressWarnings(chrom_num <- as.integer(sub("chr","",chrom)))
    # set artificially high ordering to force sex chromosomes to sort to
    # end
    if (any(grep("chrX",chrom))) { chrom_num[grep("chrX",chrom)] <- 1000 }
    if (any(grep("chrY",chrom))) { chrom_num[grep("chrY",chrom)] <- 1002 }
    if (any(grep("chrM",chrom))) { chrom_num[grep("chrM",chrom)] <- 1003 }
    idx <- order(chrom_num)
    
	#cat("that point\n")
    isFirst <- TRUE
    for (chr in chrom[idx]) {
    #cat(sprintf("%s\n",chr))
    plotChromValuePair(chr,ideoTable, bpLim=bpLim, 
        vertical=vertical, values_GR=values_GR,val_range=val_range,
        value_cols=value_cols,addScale=addScale, 
        default_margins=default_margins, plotType=plotType,col=col, 
        smoothVals=smoothVals,ablines_y=ablines_y,    
        ...)
    ctr <- ctr+1
    
    if (isFirst) 
        mtext(plot_title,side=3,outer=TRUE,line=1,cex=cex.main)
    isFirst <- FALSE; plot_title <- NULL
    
}
}
