#' Base function which plots the ideogram and superimposed data for a single chromosome. plotOnIdeo() calls this function and stacks the resulting output.
#'
#' @details Plots one unit of chromosome ideogram with dataseries 
#' superimposed. Usually, the user can avoid this function and directly
#' call plotOnIdeo(). However, this function may be used in 
#' cases where further plot customization is required.
#' @usage plotChromValuePair(chrom, cytoTable, bpLim, vertical, 
#' values_GR,val_range, col, value_cols = "values", default_margins, 
#' addScale, ablines_y, smoothVals, span=0.03, verbose = FALSE, ...)
#' @param chrom(character) chromosome(s) to create ideograms for 
#' @param cytoTable(data.frame) loaded ideogram table. 
#' (see ideoTable argument to plotOnIdeo())
#' @param bpLim(numeric) (aka xlim); display only a section of the 
#' chromosome and the corresponding values
#' @param vertical(logical) if TRUE, chromosomes will be plotted 
#' vertically
#' @param values_GR(list or GenomicRanges) If plotType is "lines" or
#' "rect", the function expects this to be a GRanges() object with 
#' data series in metadata columns. 
#' If plotType is "seg_tracks" this is a list of GRanges(), each 
#' entry of which represents a track. 
#' @param val_range(numeric) (aka ylim); y-axis scale for data series
#' @param col(character) colour for series
#' @param value_cols(character) column name for series to plot
#' @param default_margins(numeric) page inner margins (in inches)
#' @param addScale(logical) if FALSE, bp positions will be hidden
#' @param ablines_y(numeric) when specified, will draw reference lines 
#' on the y-axis
#' @param smoothVals(logical) when T applies loess() to each series
#' @param span(numeric) loess::span()
#' @param verbose(logical) print messages
#' @param ... arguments to \code{axis()}, \code{line()}, and 
#' \code{rect()}
#' @examples
#' data(hg18_ideo)
#' data(binned_multiSeries)
#' layout(matrix(1:2, byrow=TRUE,ncol=1),heights=c(2.5,1))
#' plotChromValuePair("chr1",hg18_ideo, 
#' 	values_GR=binned_multiSeries,
#'	value_cols=colnames(mcols(binned_multiSeries)),plotType='lines',
#' 	col=1:5,val_range=c(0,10),bpLim=NULL,vertical=FALSE,addScale=TRUE,
#'	ablines_y=NULL,smoothVals=FALSE,default_margins=c(0.5,.5,.1,.1))
#' @seealso plotOnIdeo()
#' @export
plotChromValuePair <- function(chrom,cytoTable,bpLim,vertical,
	values_GR,val_range,col,value_cols='values',default_margins,
	addScale,ablines_y,smoothVals,span=0.03,verbose=FALSE,...) {
  ideo <- cytoTable[cytoTable$chrom==chrom,]
  chromStart <- min(ideo$chromStart)
  chromEnd <- max(ideo$chromEnd)
  if(is.null(bpLim))  bpLim <- c(chromStart, chromEnd)	
  args <- list(ideo=ideo, 
               chromStart=chromStart, chromEnd=chromEnd, 
               bpLim=bpLim, chrom=chrom,vertical=vertical,
               addScale=addScale, ...)
  plotType <- args$plotType
#print(plotType)
  if (any(grep("plotType",names(args)))) args[["plotType"]] <- NULL
  if (any(grep("xlab",names(args)))) args[["xlab"]] <- NULL
  if (any(grep("ylab",names(args)))) args[["ylab"]] <- NULL
  
  # limit to chroms being plotted
  if (class(values_GR)=="list") {
	values_GR <- lapply(values_GR, function(x) {
		x[as.character(seqnames(x)) %in% chrom]
	})
  } else {
	values_GR <- values_GR[as.character(seqnames(values_GR)) 
	  %in% chrom]
}
  plotValue_args <- list(
	bpLim=bpLim, chrom=chrom, 
  	values_GR=values_GR,
    vertical=vertical, bty='n',val_range=val_range,
	ablines_y=ablines_y,smoothVals=smoothVals,...) 
  plotValue_args[["span"]] <- NULL

  newPlot <- TRUE
  if (vertical) {
	curr_mai <- c(0,0,default_margins[3],0)
	#curr_mai <- c(0,0,0,0)
    if (addScale) curr_mai <- replace(curr_mai,2, 0.3)
        
        par(mai=curr_mai)
        do.call(".plotChromosome", args)
cat("Plot chromosome done\n")
        ctr <- 1

		if (plotType %in% "seg_tracks") {
			plotValue_args <- modifyList(plotValue_args, 
	            	list(value_cols="", newPlot=TRUE,col=col,
						 yaxt='n'))
	            	par(mai=default_margins)
	        # cat("\tValue: "); .printMargins()
	        suppressWarnings(do.call('.plot_values',plotValue_args))
		} else {
			# plot each series in turn
	        for (k in value_cols) {
	            if (verbose) cat(sprintf('\t%s\n',k))
	            
	            col_idx <- ctr %% length(col); 
	            if (col_idx==0) col_idx=length(col)

	            plotValue_args <- modifyList(plotValue_args, 
	            	list(value_cols=k, newPlot=newPlot,
	            	col=col[col_idx]))
	            	par(mai=default_margins)
	            # cat("\tValue: "); .printMargins()
	            suppressWarnings(
					do.call('.plot_values',plotValue_args)
				)
	            newPlot <- FALSE; ctr <- ctr+1
	        }
		}
 } else { # horizontal
    ctr <- 1

	if (plotType %in% "seg_tracks") {
		plotValue_args <- modifyList(plotValue_args, 
	           	list(value_cols="", newPlot=TRUE,yaxt='n',col=col))
	           	par(mai=default_margins)
	      #cat("\tValue: "); .printMargins()
	    suppressWarnings(
			do.call('.plot_values',plotValue_args)
		)
	} else { # lines,rect
		par(mai=default_margins)
     	#cat("\tValue: "); .printMargins()
	    for (k in value_cols) {
	        if (verbose) cat(sprintf('\t%s\n',k))
	        col_idx <- ctr %% length(col); 
			if (col_idx==0) col_idx=length(col)
	
	        plotValue_args <- modifyList(plotValue_args, 
	            list(value_cols=k, newPlot=newPlot,
	                col=col[col_idx]))
	        
	        suppressWarnings(do.call('.plot_values',plotValue_args))
	        newPlot <- FALSE; ctr <- ctr+1
	    }
	}
    #abline(h=0.9,lty=3,col='black')
    
    curr_mai <- c(0,default_margins[2],0,0)
    if (addScale) { curr_mai <- replace(curr_mai,1, 0.3) }
    par(mai=curr_mai)
    # cat("\tIdeo after call: "); .printMargins()
    suppressWarnings(do.call(".plotChromosome", args))
}; 
}

.plotChromosome <- function(chrom,ideo, chromStart, chromEnd,bpLim=NULL,
titleVertical=1,titleHorizontal=4,vertical=TRUE,addScale=TRUE,
chromName_cex=0.8,verbose=FALSE,...)
{    
    chrom_width <- c(-0.1, 1.1) 
    if (verbose) cat("\n plotting...",  chrom, "\n")
    args <- list(type='n', 
        yaxt='n',xaxt='n',xlab=' ', ylab=' ', 
        bty='n',x=0,y=0,
        cex.axis=1,font.axis=2,cex.lab=1.3,font=2)

    if (vertical) {
        args <- modifyList(args, list(yaxs='i', ylim=rev(bpLim), 
            xlim=chrom_width))
            #mar=replace(def_mar,1,0.1)))
    } else {
        args <- modifyList(args, list(xaxs='i', xlim=bpLim, 
            ylim=chrom_width))
            #mar=replace(def_mar,1,0.1)))
    }
    #browser()
    suppressWarnings(do.call('plot',args))
    if (addScale) {
        if(vertical) {
            at <-2 
        } else {
            at <- 1 
        }
        axt <- axTicks(side=at)
        suppressWarnings(axis(at, xpd=FALSE,at=axt, labels=axt/1000000, ...))
    }
    
    #drawing centromere
    centro_idx <- which(ideo$gieStain=="acen") # gieStain values of "acen" 
                                               # represent the centromere
    if(any(centro_idx)){
        centro_coord <- range(c(ideo[centro_idx, "chromStart"], 
            ideo[centro_idx, "chromEnd"]))
        ideo <- ideo[-centro_idx, ]
    
        if(vertical) polygon(c(0, 1, 0, 1), c(centro_coord[1], 
            centro_coord[2], centro_coord[2], centro_coord[1]), 
            col="darkred", border="black", lwd=2) 
        else polygon(c(centro_coord[1], centro_coord[2], centro_coord[2], 
            centro_coord[1]), c(0, 1, 0, 1), col="darkred", 
            border="black", lwd=2) 
        rm(centro_coord)
        
    }; 
    rm(centro_idx)
    
    ##setting band colors 
    ## gneg & gposN(where N is an integer [1,100] ) prepresent densities. 
    ## gneg has density 0
    col <- character(nrow(ideo))
    col[1:length(col)] <- "#000000" # gvar(constitutive heterochromatins) 
                                    # and gpos100
    # colors will become darker as density increases from 0 -> 100
    col[grep("gneg", ideo$gieStain)] <- "gray90"
    col[grep("gpos25", ideo$gieStain)] <- "gray75"
    col[grep("gpos33", ideo$gieStain)] <- "gray66"
    col[grep("gpos50", ideo$gieStain)] <- "gray50"
    col[grep("stalk", ideo$gieStain)] <- "darkred" # repetitive areas
    col[grep("gpos66", ideo$gieStain)] <- "gray33"
    col[grep("gpos75", ideo$gieStain)] <- "gray25"
    ##--------------------------------------------------    
    #adding bands
    for (k in 1:nrow(ideo)) {
        if(vertical) {
            rect(xleft=0, ybottom=ideo[k,2],xright=1, ytop=ideo[k,3], 
            col=col[k], border=col[k])
        }
        else {
            rect(xleft=ideo[k,2],ybottom=0,xright=ideo[k,3],ytop=1,
            col=col[k], border=col[k])
        }
        
    }; rm(col)
    # add chrom name
    if (vertical) {
        mtext(chrom, side=titleVertical,line=2, outer=FALSE,
			col='gray50',adj=c(1,1),cex=chromName_cex,font=2)
    } else {
        mtext(chrom, side=titleHorizontal,line=0.5,outer=FALSE, 
        col='gray50', cex=chromName_cex,font=2,las=3)
    }
    
    #add border
    p_arm <- grep("p", ideo$name)
    if(any(p_arm)){
        p_rect_coord <- range(c(ideo[p_arm, "chromStart"], 
        ideo[p_arm, "chromEnd"]))
        if(vertical) rect(0, p_rect_coord[1], 1, p_rect_coord[2], 
            col=NA, border="black", lwd=2)
        else rect(p_rect_coord[1], 0, p_rect_coord[2], 1, 
            col=NA, border="black", lwd=2)
        rm(p_rect_coord)
    }
    q_arm <- grep("q", ideo$name)
    if(any(q_arm)){
        q_rect_coord <- range(c(ideo[q_arm, "chromStart"], 
            ideo[q_arm, "chromEnd"]))
        if(vertical) rect(0, q_rect_coord[1], 1, q_rect_coord[2], 
            col=NA, border="black", lwd=2)
        else rect(q_rect_coord[1], 0, q_rect_coord[2], 1, col=NA , 
            border="black", lwd=2)
        rm(q_rect_coord)
    }
    
}    

.plot_values <- function(chrom, bpLim,values_GR, vertical=TRUE, 
plotType='rect', col='gray50',val_range=NULL,  value_cols="values", 
 ablines_y=NULL, newPlot=TRUE,smoothVals=FALSE,verbose=FALSE,...) {

	# initialize plot limits
	if (plotType %in% "seg_tracks") {
		padding <- max(7-length(values_GR),1)
		val_range 	<- c(0,length(values_GR)+padding)
		rng_start 	<- start(ranges(values_GR[[1]]))
		rng_end		<- end(ranges(values_GR[[1]])) 
	} else {
		idx 	<- which(as.character(seqnames(values_GR)) == chrom) 
  		values_GR	<- values_GR[idx]
    	cvalues 	<- mcols(values_GR)[,value_cols]
    	if (is.null(val_range)) { 
			val_range <- c(min(cvalues,na.rm=TRUE),
						   max(cvalues,na.rm=TRUE)) 
		}
    	rng_start <- start(ranges(values_GR))
    	rng_end <- end(ranges(values_GR))
	} 
    
    idx <- order(rng_start)
    rng_start <- rng_start[idx]
    rng_end <- rng_end[idx]

	if (!plotType %in% "seg_tracks") {
    	cvalues <- cvalues[idx]
	}
    
	# create a new plot
    if (verbose) cat("\n plotting...", chrom, "values\n")
    if (newPlot) {
      plotDefaults <- list(
        axes=TRUE, ann=FALSE,
        xpd=FALSE,          # clip to plotting region
        xaxs="i", yaxs="i", # let R decide what works
        main="",
        xlab="", ylab="",x=0,y=0,type='n', mgp=c(3,1,0)
      )
      args <- modifyList(plotDefaults,list(...)); 
      rm(plotDefaults) #overrides default values with user inputs
    }

    if (smoothVals) {
      tmp <- loess(y~x,data=list(x=rng_start, y=cvalues), span=span)
      cvalues <- predict(tmp)
      rm(tmp)
    }

	# add data series: vertical
    if (vertical) {
        if (newPlot) {
        	args <- modifyList(args, 
				list(ylim=rev(bpLim), xlim=val_range,
            		yaxt='n',xaxt="n",xpd=NA))
        	if (any(grep("ylab",names(args)))) { 
            	ax_lab <- args$ylab; args$ylab=NULL 
        	} else {
           		ax_lab <- NULL
        	}
        
			suppressWarnings(do.call('plot', args))
        	suppressWarnings(axis(side=3,...))
        	if (!is.null(ax_lab)) mtext(ax_lab, side=3, line=2,
            	outer=FALSE,cex=0.8)
        }
        
        if (plotType=="rect") {
        	na_idx <- which(is.na(cvalues))
        	if (any(na_idx)) {
         		suppressWarnings(rect(ybottom=rng_start[-na_idx],
           			ytop=rng_end[-na_idx],xleft=0, 
            		xright=cvalues[-na_idx], border=NA,col=col,...))
        	} else {
         		suppressWarnings(rect(ybottom=rng_start,ytop=rng_end,
					xleft=0,xright=cvalues, border=NA,col=col,...))
         	}
        } else if (plotType=="lines"){
            suppressWarnings(lines(y=rng_start, 
								   x=cvalues, col=col,...))
        } else if (plotType == "seg_tracks") {
			cat("at seg_tracks\n")
						browser()	
		} else {
			stop('invalid plotType')
		}
        
        #rect(xleft=cvalues,xright=0,ybottom=rng_start, ytop=rng_end, 
        #    border=NA,...)
        #lines(y=bpLim,x=c(0,0),col='black',lty=1)
        
    # add data series: horizontal
    } else{
    	if (newPlot) {
        	args <- modifyList(args, list(xlim=bpLim,ylim=val_range,
            	xaxt='n',...))
        suppressWarnings(do.call('plot', args))
    	}
    
	    if (plotType=="rect") {
	        suppressWarnings(rect(xleft=rng_start,xright=rng_end,
	            ybottom=0, ytop=cvalues, border=NA,col=col, ...))
	    } else if (plotType=="lines"){
	        suppressWarnings(lines(x=rng_start, y=cvalues, 
								   col=col, ...))
	    } else if (plotType == "seg_tracks") {
			#cat("horiz:at seg_tracks\n")
		# plot each track
			for (k in 1:length(values_GR)) {
				idx <- k %% length(col)
				if (idx == 0) idx <- length(col)

				cury <- c(k-0.9,k-0.5)
				df <- as.data.frame(values_GR[[k]])
				rect(ybottom=rep(cury[1],nrow(df)),
					 ytop=rep(cury[2],nrow(df)),
					 xleft=df$start, xright=df$end,
					 border=NA,col=col[idx])
			}
		} else {
			stop('invalid plotType') 
		}
        lines(x=bpLim,y=c(0,0),col='black',lty=1)
    }
    
    if (!is.null(ablines_y)) {
        if (vertical) abline(v=ablines_y,lty=3,col='red')
        else abline(h=ablines_y,lty=3,col='red')
    }
#if(addBorder) box(col="black") #barplot border
}

.printMargins <- function() {
    str <- (
        sprintf("\tmar=(%s),mai=(%s),oma=(%s),omi=(%s)\n", 
        paste(round(par("mar"),digits=2),sep=",",collapse=","),
        paste(round(par("mai"),digits=2),sep=",",collapse=","),
        paste(round(par("oma"),digits=2),sep=",",collapse=","),
        paste(round(par("omi"),digits=2),sep=",",collapse=","))
    )
    print(str)
}
