#' Aggregates data by genomic bins
#'
#' @details Computed mean value of binned data. This function assumes 
#' that all elements in \code{featureData} have identical width. 
#' If provided with elements of disparate widths, the respective widths
#' are not weighted averaging. This behaviour may change in future 
#' versions of IdeoViz.
#' This function allows the user to bin data if this hasn't already 
#' been done, and is a step involved in preparing the data for 
#' \code{plotOnIdeo()}. This function computes binned within-sample 
#' average of probes overlapping the same range.  Where a range 
#' overlaps multiple bins, it gets counted in all. 
#' @usage avgByBin(xpr, featureData, target_GR, justReturnBins = FALSE, 
#'    getBinCountOnly = FALSE, FUN = mean, doSampleCor = FALSE, 
#'    verbose = FALSE)
#' @param xpr(data.frame or matrix) Locus-wise values. Rows correspond to genomic intervals (probes, genes, etc.,) while columns correspond to individual samples
#' @param featureData (data.frame or GRanges) Locus coordinates. Row order must match \code{xpr}. Column order should be: 1. chrom, 2. locus start, 3. locus end. All elements are assumed to be of identical width. Coordinates must be zero-based or one-based, but not half-open. Coordinate system must match that of \code{target_GR}.
#' @param target_GR (GRanges) Target intervals, with coordinate system matching that of featureData.
#' @param justReturnBins (logical) when TRUE, returns the coordinates of the bin to which each row belongs. Does not aggregate data in any way. This output can be used as input for more complex functions with data from each bin.
#' @param getBinCountOnly (logical) when TRUE, does not aggregate or expect xpr. Only returns number of overlapping subject ranges per bin. Speeds up computation.
#' @param FUN (function) function to aggregate data in bin
#' @param doSampleCor (logical) set to TRUE to compute mean pairwise 
#' sample correlation (Pearson correlation) for each bin; when TRUE, 
#' this function overrides \code{FUN}.
#' @param verbose (logical) print status messages
#' @return (GRanges) Binned data or binning statistics; information 
#' returned for non-empty bins only.
#' The default for this function is to return binned data; alternately,
#' if \code{justReturnBins=TRUE} or \code{getBinCountOnly=TRUE} the 
#' function will return statistics on bin counts. The latter may be
#' useful to plot spatial density of the input metric. \cr The flags 
#' and output types are presented in order of evaluation precedence:
#' \enumerate{
#'  \item{If \code{getBinCountOnly=TRUE}, returns a list with a single 
#' entry: \emph{bin_ID}: (data.frame) bin information: chrom, start, 
#' end, width, strand, index, and count. "index" is the row number of 
#' \emph{target_GR} to which this bin corresponds}
#'  \item{If \code{justReturnBins=TRUE} and 
#' \code{getBinCountOnly=FALSE}, returns a list with three entries:
#'    \enumerate{
#'      \item{\emph{bin_ID}: same as \emph{bin_ID} in output 1 above}
#'      \item{\emph{xpr}:(data.frame) \emph{B}-by-\emph{n} columns 
#' where \emph{B} is total number of [\emph{target_GR},
#' \emph{featureData}] overlaps (see next entry, \emph{binmap_idx}) 
#' and \emph{n} is number of columns in \emph{xpr}; column order 
#' matches \emph{xpr}. Contains sample-wise data "flattened" so that 
#' each [target,subject] pair is presented. More formally, entry 
#' [\emph{i},\emph{j}] contains expression for overlap of row 
#' \emph{i} from \emph{binmap_idx} for sample \emph{j} (where 1 <= 
#' \emph{i} <= \emph{B}, 1 <= \emph{j} <= \emph{n})}
#'       \item{\emph{binmap_idx}:(matrix) two-column matrix:  
#' 1) target_GR row, 2) row of featureData which overlaps with index 
#' in column 1.  (matrix output of 
#' \code{GenomicRanges::findOverlaps())})}
#'      }}
#'  \item{Default: If \code{justReturnBins=FALSE} and \code{getBinCountOnly=FALSE}, returns a GRanges object. Results are contained in the \code{elementMetadata} slot. For a dataset with \emph{n} samples, the table would have \emph{(n+1)} columns; the first column is \emph{bin_count}, and indicates number of units contained in that bin. Columns (2:(\emph{n}+1)) contain binned values for each sample in column order corresponding to that of \emph{xpr}. \cr
#'  For \code{doSampleCor=TRUE}, result is in a metadata column with name "mean_pairwise"cor". Bins with a single datapoint per sample get a value of NA.}
#'  }
#' @examples
#' ideo_hg19 <- getIdeo("hg19")
#' data(GSM733664_broadPeaks)
#' chrom_bins <- getBins(c("chr1","chr2","chrX"), 
#'		ideo_hg19,stepSize=5*100*1000)
#' # default binning
#' mean_peak <- avgByBin(data.frame(value=GSM733664_broadPeaks[,7]),  
#'	GSM733664_broadPeaks[,1:3], chrom_bins)
#' # custom function
#' median_peak <- avgByBin(
#'	data.frame(value=GSM733664_broadPeaks[,7]), 
#'	GSM733664_broadPeaks[,1:3], chrom_bins, FUN=median)
#' # mean pairwise sample correlation
#' data(binned_multiSeries)
#' bins2 <- getBins(c("chr1"), ideo_hg19, stepSize=5e6)
#' samplecor <- avgByBin(mcols(binned_multiSeries)[,1:3], binned_multiSeries, bins2, doSampleCor=TRUE)
#' # just get bin count
#' binstats <- avgByBin(data.frame(value=GSM733664_broadPeaks[,7]), 
#' GSM733664_broadPeaks[,1:3], chrom_bins, getBinCountOnly=TRUE)
#' @seealso \code{getIdeo()}, \code{getBins()}
#' @export
avgByBin <- function(xpr,featureData,target_GR,justReturnBins=FALSE,
getBinCountOnly=FALSE,FUN=mean,doSampleCor=FALSE,verbose=FALSE) {
    if (verbose) cat('* make GRanges - get overlap')
    if (class(featureData)!="GRanges") 
        probe_GR <- GRanges(featureData[,1],
            IRanges(start=featureData[,2], end=featureData[,3]))
     else 
        probe_GR <- featureData
    
    target_GR <- target_GR[as.character(seqnames(target_GR)) %in% 
						   as.character(seqlevels(probe_GR))]
    seqlevels(target_GR) <- seqlevels(probe_GR)
    
    ol <- as.matrix(findOverlaps(target_GR, probe_GR))
    
    if (verbose) cat('\n* Taking averages of windows')
    avg_DF <- as.data.frame(target_GR[ol[,1]])
    
    col2use <- "index"
    avg_DF <- cbind(avg_DF, index=ol[,1])
    
    tmp <- ave(rep(1,nrow(avg_DF)), avg_DF[,col2use],FUN=sum)
    avg_DF <- cbind(avg_DF, bin_count=tmp)
    rm(tmp)
    
    if (getBinCountOnly)  return(list(bin_ID=avg_DF)) 
    
    # this line is key to matching intervals contained in the 
    # respective target ranges
    xpr2 <- data.frame(xpr[ol[,2],])
    colnames(xpr2) <- colnames(xpr)
    
    if (justReturnBins & !getBinCountOnly) 
        return(list(bin_ID=avg_DF,xpr=xpr2,binmap_idx=ol))
    
    if (!doSampleCor){
    for (k in 1:ncol(xpr2)) { 
        tmp <- ave(xpr2[,k], avg_DF[,col2use], FUN=FUN)
        avg_DF <- cbind(avg_DF,tmp)
    }
    colnames(avg_DF)[which(colnames(avg_DF)=="tmp")] <- colnames(xpr2)
    
    } else {
        if (class(xpr2) =="matrix") xpr2 <- as.data.frame(xpr2)
        corFunc <- function(x) {
        if (nrow(x) ==1) return(NA)# cannot do intersample correlation 
                                   # with one datapoint in the bin
        x <- as.matrix(x)
        y <- matrix(x, ncol=ncol(xpr2), byrow=FALSE)
        z <- cor(y,method='p')
        diag(z) <- NA
        cors <- z[upper.tri(z,diag=FALSE)]
        return(mean(cors))
        } 
        tmp <- ave(xpr2, avg_DF[,col2use], FUN=corFunc )
        avg_DF <- cbind(avg_DF,tmp[,1]); 
        colnames(avg_DF)[ncol(avg_DF)] <- "mean_pairwise_cor"
    }
    
    avg_DF <- avg_DF[!duplicated(avg_DF[,col2use]),]
    rm(xpr2,target_GR)
    
    avg_GR <- GRanges(avg_DF[,1],IRanges(avg_DF[,2],avg_DF[,3]))
    mcols(avg_GR) <- avg_DF[,c(7:ncol(avg_DF))]
    
    return(avg_GR)
}

