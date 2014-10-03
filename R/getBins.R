getBins <-function(chroms, ideo, binLim=NULL,stepSize) 
{
    # make intervals for averaging
    if (missing(chroms)) chroms <- unique(ideo[,1])
    chromSteps <- cbind(CHROM=NA, START=NA,END=NA)
    
    for (chr in chroms) {
        csize <- max(ideo$chromEnd[ideo$chrom==chr])
        bstart <- seq(1,csize, stepSize)
        bend <- bstart + (stepSize-1); bend[length(bend)] <- csize
        df <- cbind(chr, bstart,bend)
        chromSteps <- rbind(chromSteps,df)
    }
    chromSteps <- chromSteps[-1,]
    full_chromInt <- GRanges(seqnames=chromSteps[,1], 
        ranges=IRanges(start=as.numeric(chromSteps[,2]),
        end=as.numeric(chromSteps[,3])))
    
    return(full_chromInt)
}
