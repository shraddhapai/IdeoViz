avgByBin <- function(xpr,featureData,target_GR,justReturnBins=FALSE,
getBinCountOnly=FALSE,FUN=mean,doSampleCor=FALSE,verbose=FALSE)
{
    if (verbose) cat('* make GRanges - get overlap')
    if (class(featureData)!="GRanges") 
        probe_GR <- GRanges(featureData[,1],
            IRanges(start=featureData[,2], end=featureData[,3]))
     else 
        probe_GR <- featureData
    
    target_GR <- target_GR[seqnames(target_GR) %in% seqlevels(probe_GR)]
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

