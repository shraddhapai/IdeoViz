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
         bpLim <- c(min(ideoTable$chromStart),max(ideoTable$chromEnd)); 
        }
    }
    
    options(scipen=10) # don't convert large numbers to scientific notation

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
