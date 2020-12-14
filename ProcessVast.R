#!/usr/bin/env Rscript
### Work up vast-tools output, save tables and do a few analysis of changing events
###
### U. Braunschweig 2019-2020
###
### Changes:  Option 'depth' has been replaced with 'filter', with the default settings now allowing
###           a quality score of OK/LOW/VLOW for CE/MIC depth and OK/B1/B2/Bl/Bn for balance, 
###           and IR minimum reads of 10. This rescues ~10-20% of events compared to the previous default.

cArgs <- commandArgs(TRUE)

## Check dependencies
libMissing <- !require("optparse", quietly=T) && stop("Failed to load R package 'optparse'")


### Capture input
opt.list <- list(
    make_option(c("-s", "--sampleTab"),     action="store", type="character",      metavar="CSVFILE",
                help="Sample table in CSV format with, at least, columns Sample (an arbitrary, unique name),
                and Type (treatment etc.)"),
    make_option(c("-c", "--contrTab"),      action="store", type="character",      metavar="CSVFILE",
                help="Contrasts table in CSV format with, at least, columns Experimental, Control, and File."),
    make_option(c("-p", "--specEventTable"),action="store", type="character",      metavar="CSVFILE",
                help="Optional table indicating specific events to be used for an analysis of changing
                events. Required are a colsumn EVENTS and one or more other columns containing a string
                annotation such as 'neural-HIGH', 'neural-LOW' etc."),
    make_option(c("-v", "--vastDir"),       action="store", default="./vast", metavar="DIR",
                help="Directory with vast-tools output [%default]"),
    make_option(c("-o", "--outDir"),        action="store", default="./splicing", metavar="DIR",
                help="Output directory. Will be created if it does not exist [%default]"),
    make_option("--continue",               action="store_true", default="FALSE",
                help="Try to locate saved output tables in --outDir and load them rather then doing full analysis"),
    make_option(c("-E", "--vast2EnsemblGene"),action="store", type="character",      metavar="FILE",
                help="Optional tab-delimited file with required for outputting files for GO analysis,
                with at least columns EVENT and EnsemblGeneID"),
    make_option("--noGO",                   action="store_true", default="FALSE",
                help="Do not save files for GO analysis."),
    make_option(c("-m", "--minRepFrac"), action="store", default=0.5, metavar="REAL",
                help="Minimum fraction of replicates of each type that needs to survive filtering [%default]"),
    make_option(c("-f", "--filter"),        action="store", default="DEFAULT", metavar="DEFAULT|STRICT|LEGACY",
                help="Read depth filtering and balance filtering for individual events and replicates.
		DEFAULT: CE/MIC/Alt3/Alt5 depth must be SOK/OK/LOW/VLOW; CE balance OK/B1/B2/Bl/Bn; IR depth >= 10; 
                         IR balance pval >= 0.05.
                STRICT:  CE/MIC/Alt3/Alt5 depth must be SOK/OK/LOW; CE balance OK/B1/B2/Bl/Bn; IR depth >= 15; 
                         IR balance pval >= 0.05.
                LEGACY:  CE/MIC/Alt3/Alt5 depth must be SOK/OK/LOW; CE balance OK/B1; IR depth >= 15; IR balance pval >= 0.05.
                         This is the same as the previous default setting.
                [%default]"),
    make_option("--strictDiff",             action="store_true", default="FALSE",
                help="If set, significant events will be those where |dPSI| (p>0.95) is greater than 10.
                By default, significant events will be those where |dPSI| (p>0.95) is greater
                than 0 and the |dPSI| point estimate is greater than 5/10/15 (10 for plots)"),
    make_option("--keepFullAlt",            action="store_true", default="FALSE",
                help="By default, Alt5 and Alt3 events will be reduced to one score for each event,
                with 0 indicating full usage of the closest alternative splice site and 1 only the furthest.
                When this option is set, tables will be stored as an .Robj that still contain the original
                separate sub-events."),
    make_option("--scatterForce",           action="store_true", default="FALSE",
                help="By default, if there are > 9 contrasts, mutual pairwise scatterplots of dPSI by event
                type are suppressed. Set this flag to force plotting."),
    make_option("--colors",                 action="store",
                default="dodgerblue,grey50,darkmagenta,darkorange1,darkolivegreen2,darkslategray4,burlywood3,cyan3,darkgoldenrod3,firebrick3,navy,seagreen4",
                metavar="NUM",
                help="Colors per sample type for plotting [%default]"),   
    make_option(c("-C", "--cores"),         action="store", default=1, metavar="INT",
                help="CPUs to use [%default]")
)

opt <- parse_args(OptionParser(
    usage = "Filter vast-tools results, average replicates, generate tables, analyze differential events",
    option_list=opt.list),
    args=cArgs)

libMissing <- !require("gplots", quietly=T)   && stop("Failed to load R package 'gplots'")
libMissing <- !require("parallel", quietly=T) && {
    warning("Failed to load R package 'parallel', running with a single core.")
    opt$cores <- 1
}


### Constants ###

maxPlotContr <- 9

### Function defs ###
saveLog <- function(logName, opt, sampleTab, contrTab) {
### Save a log file
    titl <- "*** Processing vast-tools results ***\n\n== Options =="
    if (opt$continue) {
        selOpt <- which(names(opt) %in% c("sampleTab","contrTab","vast2EnsemblGene","specEventTable","continue",
                                          "scatterForce","outDir","noGO")
                        )
    } else {
        selOpt <- which(names(opt) %in% c("sampleTab","contrTab","vast2EnsemblGene","specEventTable","continue",
                                          "strictDiff","keepFullAlt","minRepFrac","filter",
					  "scatterForce","vastDir","outDir","noGO")
                        )
    }
    optTable <- data.frame(option = names(opt)[selOpt],
                           value  = unlist(opt[selOpt]),
                           stringsAsFactors=F
                           )
    optTable$option <- paste0(optTable$option,
                              sapply(max(nchar(optTable$option)) + 1 - nchar(optTable$option), FUN=function(x) {
                                  paste(rep(" ", x), collapse="")
                              })
                              )
    suppressWarnings({
        write.table(titl,
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName)
        write.table(paste(optTable$option, optTable$value, sep=": "),
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Samples ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table(sampleTab[,names(sampleTab) %in% c("Sample","Type")],
                    row.names=F, col.names=T, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Contrasts ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table(contrTab[,names(contrTab) %in% c("Experimental","Control","File")],
                    row.names=F, col.names=T, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Time stamp ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
    })
}

cleanAS <- function(x,                                        # Quality column of one sample
                    complexCol,                               # Colum 'COMPLEX' of vast-tools output
                    covThresCE    = c("SOK","OK","LOW"),      # Allowable coverage score (score 3) for CE
                    covThresOther = c("SOK","OK","LOW"),      # Allowable coverage score (score 3) for ME, Alt3, Alt5
                    covThresIR    = 15,                       # Minimum coverage (score 3) for IR
                    balThresCE    = c("OK","B1"),             # Allowable balance score (score 4) for CE
                    balThresIR    = 0.05                      # Maximum balance p-value (score 4, bigger is better) for IR
                    ) {
### Filter PSI output of vast-tools pipeline using the quality column
### Value: logical vector indicating which events passed the criteria
    
    is.ce <- complexCol %in% c("S","C1","C2","C3","ANN")
    is.ir <- grepl("IR", complexCol)
    is.o  <- complexCol %in% c("Alt5", "Alt3", "MIC")

    score1 <- sub("([^,]+),.+", "\\1", x)
    score3 <- sub("[^,]+,[^,]+,([^,]+),.+", "\\1", x)
    score4 <- sub("[^,]+,[^,]+,[^,]+,([^,]+),.+", "\\1", x)
    score5 <- sub("[^,]*,[^,]*,[^,]*,[^,]*,([^,@]*)[@,].*", "\\1", x)

    oldVast <- ifelse(any(grepl("=", score3)), FALSE, TRUE)
    cat("Vast-tools version", ifelse(oldVast, "<", ">="), "2.2.2 detected", "\n")

    if (oldVast) {  # The content of score4 has changed in vast-tools 2.2.2
        ok.ce <- is.ce & score3 %in% covThresCE & score4 %in% balThresCE
        ok.o  <- is.o  & score3 %in% covThresOther
    } else {
        ok.ce <- is.ce & score1 %in% covThresCE & score4 %in% balThresCE
        ok.o  <- is.o  & score1 %in% covThresOther
    }

    suppressWarnings(
        nreads <- as.integer(sub("([^=]*)=.+", "\\1", score4)) +
            as.integer(sub("[^=]*=([^=]*)=.*", "\\1", score4)) +
            as.integer(sub("[^=]*=[^=]*=(.*)", "\\1", score4))
    )
    suppressWarnings(
        ok.ir <- is.ir & nreads >= covThresIR & as.numeric(score5) > balThresIR
    )

    ok.ce | ok.ir | ok.o
}

removeVastDupEvents <- function(vast) {
### For vast-tools version with duplicate event IDs, remove duplicates, prioritizing MIC events
    cat("Duplicated vast-tools EVENT IDs found... removing.\n")
    vast[!(vast$EVENT %in% vast$EVENT[duplicated(vast$EVENT)]),]
}

reshapeVast <- function(vast) {
### Re-introduce dummy 'quality' columns into filtered INCLUSION... table
    tmp <- matrix(nrow=nrow(vast), ncol=(ncol(vast) - 6) * 2)
    colnames(tmp) <- paste0(rep(names(vast)[-c(1:6)], each=2), c("","-Q"))
    tmp[,seq(1, ncol(tmp) - 1, 2)] <- as.matrix(vast[,7:ncol(vast)])
    data.frame(vast[,1:6], tmp)
}

mergeToVast <- function(x, events) {
### Merge DIFF tables to conform to the vast-tools main table
    x <- merge(x, data.frame(EVENT=vast$EVENT, vastInd=1:nrow(vast)), by.x=2, by.y=1, all.y=T)
    x <- x[order(x$vastInd), c(2,1,3:6)]
}


diffPointEst <- function(x, contrTab, diffres, sampleTab, vast, minFrac=0.5) {
### Extract PSI point estimates from DIFF tables, combine if multiple,
### and filter by minimum number of sample that survived quality filtering.
### x: sample type
    cat(x, "\n")
    tmp1 <- tmp2 <- diffres[[1]][,c()]
    if (length(grep(x, contrTab$Experimental)) > 0) {
        tmp1 <- sapply(which(contrTab$Experimental == x), FUN=function(y) {sapply(diffres[y], "[[", 3)})
    }
    if (length(grep(x, contrTab$Control)) > 0) {
        tmp2 <- sapply(which(contrTab$Control == x), FUN=function(y) {sapply(diffres[y], "[[", 4)})
    }

    if (!is.numeric(tmp1)) {tmp1 <- as.numeric(rep(NA, nrow(tmp2)))}
    if (!is.numeric(tmp2)) {tmp2 <- as.numeric(rep(NA, nrow(tmp1)))}
    tmp <- 100 * rowMeans(cbind(tmp1, tmp2), na.rm=T)

    filt <- vast[,sampleTab$vastCol[sampleTab$Type == x]]
    if (any(is.null(dim(filt)))) {filt <- as.matrix(filt)}
    minNotNa <- ceiling(ncol(filt) * minFrac)
    tmp[apply(filt, MAR=1, FUN=function(y) {length(which(!is.na(y)))}) < minNotNa] <- NA
    tmp
}

dpsiMake <- function(i, contrTab, uqContr, means) {
### Generate dPSI from DIFF point estimates
    means[,which(uqContr == contrTab$Experimental[i])] - means[,which(uqContr == contrTab$Control[i])]
}

altCombinePSI <- function(vast, dat) {
### Combine mean PSI/PIR values of Alt5/Alt3 into a compound score
    info <- data.frame(vast[,1:6], vastInd=1:nrow(vast))[order(vast$EVENT),]
    dat  <- as.matrix(dat[info$vastInd,])
    info$EVENT <- sub("(.+ALT[AD][0-9]+)-.*", "\\1", info$EVENT)
    
    dat.a5 <- dat[info$COMPLEX == "Alt5",]
    info.a5 <- info[info$COMPLEX == "Alt5",]
    posInd.a5 <- altIndex(dat.a5, uniqEvent=info.a5$EVENT, aggreg="altInd")
    info.a5 <- info.a5[nrow(info.a5):1,]
    info.a5 <- info.a5[!duplicated(info.a5$EVENT),]
    info.a5 <- info.a5[!is.na(info.a5$EVENT),]
    info.a5 <- info.a5[nrow(info.a5):1, 1:6]

    dat.a3 <- dat[info$COMPLEX == "Alt3",]
    info.a3 <- info[info$COMPLEX == "Alt3",]
    posInd.a3 <- altIndex(dat.a3, uniqEvent=info.a3$EVENT, aggreg="altInd")
    info.a3 <- info.a3[nrow(info.a3):1,]
    info.a3 <- info.a3[!duplicated(info.a3$EVENT),]
    info.a3 <- info.a3[!is.na(info.a3$EVENT),]
    info.a3 <- info.a3[nrow(info.a3):1, 1:6]

    out <- rbind(data.frame(info[,-7], dat)[!(info$COMPLEX %in% c("Alt5","Alt3")),],
                 data.frame(info.a5, posInd.a5[,-1]),
                 data.frame(info.a3, posInd.a3[,-1])
                 )

    type <- character(nrow(out))
    type[out$COMPLEX %in% c("C1","C2","C3","S","ANN")] <- "CE"
    type[grep("IR", out$COMPLEX)] <- "IR"
    type[out$COMPLEX == "MIC"]  <- "MIC"
    type[out$COMPLEX == "Alt5"] <- "Alt5"
    type[out$COMPLEX == "Alt3"] <- "Alt3"
    type[is.na(out$COMPLEX) & grepl("EX[0-9]+", out$EVENT) & out$LENGTH >  27] <- "CE"  # some events have no COMPLEX
    type[is.na(out$COMPLEX) & grepl("EX[0-9]+", out$EVENT) & out$LENGTH <= 27] <- "MIC" #
    type[type == "" & grepl("A_", out$COMPLEX)] <- "CE"
    type[type == "CE" & out$LENGTH <= 27]  <- "MIC"

    
    data.frame(out[,1:6],
               TYPE = type,
               out[,7:ncol(out)]
               )
}

altCombine.sig <- function(vast, diffsig) {
### Convert DIFF results such that for Alt5/Alt3, only the biggest MV.dPSI is taken
    info <- data.frame(vast[,1:6], vastInd=1:nrow(vast))[order(vast$EVENT),]
    dat <- as.matrix(diffsig[info$vastInd,])
    info$EVENT <- sub("(.+ALT[AD][0-9]+)-.*", "\\1", info$EVENT)

    dat.a5 <- as.matrix(dat[info$COMPLEX == "Alt5",])
    info.a5 <- info[info$COMPLEX == "Alt5",]
    anySig.a5 <- altIndex(dat.a5, uniqEvent=info.a5$EVENT, aggreg="max")
    rownames(anySig.a5) <- anySig.a5[,1]

    dat.a3 <- as.matrix(dat[info$COMPLEX == "Alt3",])
    info.a3 <- info[info$COMPLEX == "Alt3",]
    anySig.a3 <- altIndex(dat.a3, uniqEvent=info.a3$EVENT, aggreg="max")
    rownames(anySig.a3) <- anySig.a3[,1]

    rbind(as.matrix(dat[!(info$COMPLEX %in% c("Alt5","Alt3")),]),
          as.matrix(anySig.a5[,-1]),
          as.matrix(anySig.a3[,-1])
	    ) 
}

altIndex <- function(dat, uniqEvent, aggreg=c("altInd","max","any")[1]) {
    posInd <- aggregate(dat,
                        by  = list(uniqEvent),
                        FUN = switch(aggreg,
                                     altInd = function(x) {mean(x[-1] * 1:(length(x) - 1))},
                                     max    = max,
                                     any    = any)
                        )
    posInd[!(posInd[,1] %in% uniqEvent[duplicated(uniqEvent)]), -1] <- NA # only one subevent
    posInd
}

orderEvents <- function(info) {
### Generate a sorting for events
    ind <- numeric(nrow(info))
    ind[info$TYPE == "CE"]   <- 1
    ind[info$TYPE == "MIC"]  <- 2
    ind[info$TYPE == "IR"]   <- 3
    ind[info$TYPE == "Alt5"] <- 4
    ind[info$TYPE == "Alt3"] <- 5
    order(ind, info$GENE, info$EVENT)    
}

setRowCol <- function(x) {
### Select a layout for x panels
    if (x == 1)                ncol=1
    if (x %in% c(2,4))         ncol=2
    if (x %in% c(3,5,6,9))     ncol=3
    if (x %in% c(7,8) | x > 9) ncol=4
    nrow <- ceiling(x / ncol)
    c(nrow, ncol)
}

changingEventTypes <- function(dpsi, sig, info, includeANN=T, main=NA, percentage=T,
                               cols=c("gold","dodgerblue"),
                               ylim=NA, cex.main=1.2) {
### Barplot of up- and downregulated events, by event type
###   dpsi       :  vector of dPSI
###   sig        :  logical vector indicating which dpsi are signficant
###   eventCol   :  column EVENT from vast-tools splicing table
###   complexCol :  column COMPLEX from vast-tools splicing table
###   lengthCol  :  column LENGTH -"-
###   MICsize    :  maxium size for cassette exons to be classified as microexons
###   includeANN :  Include 'ANN' events, which are mostly constitutive exons, in counts for CE and MIC?
###   main       :  (optional) title
### Changes:     New option 'includeANN'

    if (any(is.na(main) | length(main) != ncol(dpsi))) {main <- rep("", ncol(dpsi))}
    if (!includeANN) {info$TYPE[info$COMPLEX == "ANN"] <- "ANN"}

    ## Count numbers of events
    change <- array(dim=c(2, length(setdiff(unique(info$TYPE), NA)), ncol(dpsi)),
                    dimnames=list(c("up","down"), c("CE","MIC","Alt3","Alt5","IR"), c()))
    allev <- matrix(nrow=ncol(dpsi), ncol=length(setdiff(unique(info$TYPE), NA)))
    for (i in 1:ncol(dpsi)) {
        change[1,,i] <- sapply(dimnames(change)[[2]],
                               FUN=function(x) {length(which(dpsi[sig[,i] & info$TYPE == x, i] > 0))})
        change[2,,i] <- sapply(dimnames(change)[[2]],
                               FUN=function(x) {length(which(dpsi[sig[,i] & info$TYPE == x, i] < 0))})
        allev[i,] <- sapply(dimnames(change)[[2]],
                            FUN=function(x) {length(which(!is.na(dpsi[info$TYPE == x, i])))})
    }
    if (percentage) {
        ylab <- "% of changing events"
        plotchange <- change
        for (i in 1:ncol(dpsi)) {
            plotchange[,,i] <- 100 * t(t(plotchange[,,i]) / allev[i,])
        }
    } else {
        ylab <- "Number of changing events"
        plotchange <- change
    }
    if (any(is.na(ylim))) {
        ylim <- c(-max(plotchange[2,,], na.rm=T), max(plotchange[1,,], na.rm=T))
	ylim[1] <- min(ylim[1], -0.1 * (ylim[2] - ylim[1]))  # make sure the labels fit
	ylim[2] <- max(ylim[2],  0.1 * (ylim[2] - ylim[1]))
    }


    ## Plotting
    rowCol <- setRowCol(ncol(dpsi))
    par(mfrow=rowCol, mar=c(3,4,3,1))
    for (i in 1:ncol(dpsi)) {
        plot(1,1, type="n", bty="n", xaxt="n", xlim=c(0.3,7.7), ylim=ylim,
             xlab="", ylab=ylab, main=main[i], cex.main=cex.main)
        axis(1, at=seq(1, 7 ,1.5), tick=F, line=F, labels=colnames(change), cex.axis=1.2)
        rect(seq(0.5, 6.5, 1.5),          0, seq(1.5, 7.5, 1.5), plotchange[1,,i], col=cols[1], border=cols[1])
        rect(seq(0.5, 6.5, 1.5), -plotchange[2,,i], seq(1.5, 7.5, 1.5),         0, col=cols[2], border=cols[2])
        abline(h=0)
        text(seq(1, 7 ,1.5),  0.05 * (ylim[2] - ylim[1]), change[1,,i])
        text(seq(1, 7 ,1.5), -0.05 * (ylim[2] - ylim[1]), change[2,,i])
        par(xpd=NA)
        text(4, ylim[2] + 0.06 * (ylim[2] - ylim[1]) ,
             paste(format(sum(allev[i,]), big.mark=","), "events with coverage"))
        par(xpd=F)
    }
}

changingSpecEvents <- function(specEv, dpsi.i, sig.i, info, includeANN=T, percentage=T, main=NA,
                               ylim=NA, cex.main=1.2) {
### Plot the % of events among ones submitted in an extra file as bar graphs,
### with one panel per grouping (e.g. neural-differential)
    specEv <- merge(data.frame(EVENT = info$EVENT, ind=1:nrow(info)),
                    specEv,
                    by.x=1, by.y=which(names(specEv) == "EVENT"), all.x=T)
    specEv <- as.data.frame(specEv[order(specEv$ind),-c(1,2)])

    ## Count numbers of events
    if (!includeANN) {info$TYPE[info$COMPLEX == "ANN"] <- "ANN"}

    specKinds <- apply(specEv, MAR=2, getSpecKinds)
    nSpecKinds <- sapply(specKinds, nrow)

    dat <- lapply(1:ncol(specEv), FUN=function(x) {
        change <- array(dim=c(2, length(setdiff(unique(info$TYPE), NA)), nSpecKinds[x]),
                              dimnames=list(c("up","down"), c("CE","MIC","Alt3","Alt5","IR"), specKinds[[x]]$kind))
        allev <- matrix(nrow=nSpecKinds[x], ncol=length(setdiff(unique(info$TYPE), NA)))
        for (i in 1:nSpecKinds[x]) {
            change[1,,i] <- sapply(dimnames(change)[[2]],
                                   FUN=function(y) {length(which(dpsi.i[sig.i & info$TYPE == y &
                                                                      specEv[,x] == specKinds[[x]]$kind[i]] > 0))})
            change[2,,i] <- sapply(dimnames(change)[[2]],
                                   FUN=function(y) {length(which(dpsi.i[sig.i & info$TYPE == y &
                                                                      specEv[,x] == specKinds[[x]]$kind[i]] < 0))})
            allev[i,] <- sapply(dimnames(change)[[2]],
                            FUN=function(y) {length(which(!is.na(dpsi.i[info$TYPE == y &
                                                                      specEv[,x] == specKinds[[x]]$kind[i]])))})
        }
        if (percentage) {
            ylab <- "% of changing events"
            plotchange <- change
            for (i in 1:nSpecKinds[x]) {
                plotchange[,,i] <- 100 * t(t(plotchange[,,i]) / allev[i,])
            }
        } else {
            ylab <- "Number of changing events"
            plotchange <- change
        }
        
        list(ylab       = ylab,
             change     = change,
             allev      = allev,
             plotchange = plotchange
             )        
    })
        
    if (any(is.na(ylim))) {
        ylim <- c(-max(sapply(dat, FUN=function(x) {max(x$plotchange[2,,], na.rm=T)})),
                   max(sapply(dat, FUN=function(x) {max(x$plotchange[1,,], na.rm=T)})))
		   ylim[1] <- min(ylim[1], -0.1 * (ylim[2] - ylim[1]))  # make sure the labels fit
		   ylim[2] <- max(ylim[2],  0.2 * (ylim[2] - ylim[1]))
    }
    
    specRowCol <- setRowCol(ncol(specEv))
        par(mfrow=specRowCol, mar=c(3,4,3,1))
    for (i in 1:ncol(specEv)) {
        xpos <- seq(1, 7, 1.5)
        xfine <- seq(0, 1, length.out=1 + nSpecKinds[i]) - 0.5
        xnumb <- (xfine[-1] + xfine[-length(xfine)]) / 2
        plot(1,1, type="n", bty="n", xaxt="n", xlim=c(0.3,7.7), ylim=ylim,
             xlab="", ylab=dat[[i]]$ylab, cex.main=cex.main)
        if (i == 1) {title(main=main, adj=0, cex.main=1.5)}
        axis(1, at=xpos, tick=F, line=F, labels=colnames(dat[[1]]$change), cex.axis=1.2)
        for (j in 1:nSpecKinds[i]) {
            rect(xpos + xfine[j],                           0, xpos + xfine[j + 1], dat[[i]]$plotchange[1,,j],
                 col=specKinds[[i]]$col[j], border=specKinds[[i]]$col[j])
            rect(xpos + xfine[j],  -dat[[i]]$plotchange[2,,j], xpos + xfine[j + 1],                         0,
                 col=specKinds[[i]]$col[j], border=specKinds[[i]]$col[j])
            text(xpos + xnumb[j],  0.03 * (ylim[2] - ylim[1]), dat[[i]]$change[1,,j], srt=90, adj=c(0,0.5))
            text(xpos + xnumb[j], -0.03 * (ylim[2] - ylim[1]), dat[[i]]$change[2,,j], srt=90, adj=c(1,0.5))
            text(xpos + xnumb[j], ylim[1], dat[[i]]$allev[j,], srt=90, adj=c(0,0.5))
        }
        abline(h=0)

        legend("topright", legend=specKinds[[i]]$kind, fill=specKinds[[i]]$col, border="white", bty="n")
    }   
}

getSpecKinds <- function(specEv1) {
### Called by changingSpecEvents(), extracts the levels of a grouping and assigns colors
    reds  <- c("brown1","brown4")
    blues <- c("dodgerblue","darkblue")
    cols  <- c("dodgerblue","darkmagenta","darkorange1","darkolivegreen2","darkslategray4","burlywood3",
               "cyan3","darkgoldenrod3","firebrick3","navy")
    kinds <- sort(unique(specEv1))
    kinds <- data.frame(kind=kinds, col=NA)
    kinds$col[grep("HIGH",kinds$kind)] <- colorRampPalette(reds)(length(grep("HIGH",kinds$kind)))
    kinds$col[grep("LOW",kinds$kind)]  <- colorRampPalette(blues)(length(grep("LOW",kinds$kind)))
    if (any(is.na(kinds$col))) {
        kinds$col[is.na(kinds$col)] <- rep(cols, ceiling(length(which(is.na(kinds$col))) / length(cols)))[
            1:length(which(is.na(kinds$col)))]
    }

    kinds
}

scatterXY <- function(x,y, cols=c("dodgerblue","brown1","darkmagenta","darkolivegreen3")) {
    x <- which(colnames(dpsi) == x)
    y <- which(colnames(dpsi) == y)
    notna <- !is.na(dpsi[,x]) & !is.na(dpsi[,y])
    change <- which((diff10[,x] | diff10[,y]) & notna)

    par(mfrow=c(1,5), mar=c(5,4,3,1))
    for (i in unique(info$TYPE)) {
        evType <- info$TYPE %in% i
        use <- which(notna & evType)
        lims <- range(as.numeric(c(dpsi[use, x], dpsi[use, y])), na.rm=T)
        sigcols <- rep(NA, length(use))
        sigcols[diff10[use, x]]  <- cols[1]
        sigcols[diff10[use, y]] <- cols[2]
        sigcols[diff10[use, x] & diff10[use, y]] <- cols[3]
        sigcols[diff10[use, x] & diff10[use, y] & sign(dpsi[use, x]) != sign(dpsi[use, y])] <- cols[4]
        

        smoothScatter(dpsi[use[is.na(sigcols)],x], dpsi[use[is.na(sigcols)],y], pch=".", xlim=lims, ylim=lims,
                      colramp=colorRampPalette(c("white","grey60","grey10","grey3","black")),
                      bandwidth=c(0.05,0.05), nbin=500,
                      nrpoints=0, transformation=function(x) {x^0.15},
                      xlab=colnames(dpsi)[x], ylab=colnames(dpsi)[y])
        points(dpsi[use,x], dpsi[use,y], pch=20, col=sigcols)
 
        title(main=paste(format(length(use), big.mark=","), i, "events"), adj=1)

        abline(a=0, b=1, lty=2)
        
        legend("topleft",
               legend=c("Significant changes (dPSI>10):",
                        paste(c(colnames(dpsi)[x], colnames(dpsi)[y],"both","divergent"), " (",
                              sapply(cols, FUN=function(x) {length(which(sigcols == x))}), ")", sep="")),
               pch=20, col=c(NA,"dodgerblue","brown1","darkmagenta","darkolivegreen3"), bty="n")

        legend("bottomright",
               legend=paste("r =", round(cor(dpsi[use,x], dpsi[use,y], use="c"), 2)),
               bty="n")
    }    
}

cleanV2G <- function(v2g, info) {
    if (any(grepl("ALT[AD][0-9]+-[0-9]+/", v2g$EVENT))) {
        v2g$EVENT <- sub("(.+ALT[AD][0-9]+)-.*", "\\1", v2g$EVENT)
        v2g <- v2g[!duplicated(v2g$EVENT),]
    }
    v2g <- merge(data.frame(EVENT=info$EVENT, ind=1:nrow(info)),
                 v2g,
                 by.x=1, by.y=grep("EVENT", names(v2g)), all.x=T
                 )
    v2g <- v2g[order(v2g$ind), -2]
    if (length(which(sub("-[0-9]+/[0-9]+", "", info$EVENT[info$COMPLEX != "ANN"]) %in%
                     sub("-[0-9]+/[0-9]+", "", v2g$EVENT))) / length(which(info$COMPLEX != "ANN")) < 0.9) {
        stop("Less than 90% of vast-tools events found in ", basename(opt$specEventTable))
    }
    v2g
}

saveGOfiles <- function(v2g, info, opt, contrTab, diff, dpsi) {
    GOoutPath <- file.path(opt$outDir, "GO", "input")
    for (i in 1:nrow(contrTab)) {
        write.table(getBG(i,c("MIC","CE"), diff), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("BG_", contrTab$Contrast[i], ".CE.MIC.txt")))
        write.table(getBG(i,"MIC", diff), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("BG_", contrTab$Contrast[i], ".MIC.txt")))
        write.table(getBG(i,"IR", diff), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("BG_", contrTab$Contrast[i], ".IR.txt")))
        write.table(getBG(i,"Alt5", diff), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("BG_", contrTab$Contrast[i], ".Alt5.txt")))
        write.table(getBG(i,"Alt3", diff), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("BG_", contrTab$Contrast[i], ".Alt3.txt")))

        write.table(getFG(i,c("MIC","CE"), diff,"up"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".CE.MIC.up.txt")))
        write.table(getFG(i,c("MIC","CE"), diff,"dn"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".CE.MIC.dn.txt")))
        write.table(getFG(i,"MIC", diff,"up"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".MIC.up.txt")))
        write.table(getFG(i,"MIC", diff,"dn"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".MIC.dn.txt")))
        write.table(getFG(i,"IR", diff,"up"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".IR.up.txt")))
        write.table(getFG(i,"IR", diff,"dn"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".IR.dn.txt")))
        write.table(getFG(i,"Alt5", diff,"up"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".Alt5.up.txt")))
        write.table(getFG(i,"Alt5", diff,"dn"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".Alt5.dn.txt")))
        write.table(getFG(i,"Alt3", diff,"up"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".Alt3.up.txt")))
        write.table(getFG(i,"Alt3", diff,"dn"), row.names=F, col.names=F, quote=F,
                    file=file.path(GOoutPath, paste0("AS_", contrTab$Contrast[i], ".Alt3.dn.txt")))      
    }
}
   
getBG <- function(x, kind, diff) {
    changing <- as.character(unique(v2g$EnsemblGeneID[which(diff[,x] & info$TYPE %in% kind)]))
    setdiff(unique(c(changing, as.character(v2g$EnsemblGeneID[!is.na(dpsi[,x]) & info$TYPE %in% kind]))), NA)
}
getFG <- function(x, kind, diff, way=c("up","dn", "change")) {
    if (way == "up") {
        return(setdiff(unique(v2g$EnsemblGeneID[which(dpsi[,x] >  0 & diff[,x] & info$TYPE %in% kind)]), NA))
    }
    if (way == "dn") {
        return(setdiff(unique(v2g$EnsemblGeneID[which(dpsi[,x] < 0 & diff[,x] & info$TYPE %in% kind)]), NA))
    }
    if (way == "change") {
        return(setdiff(unique(v2g$EnsemblGeneID[which(diff[,x] & info$TYPE %in% kind)]), NA))
    }
}


### Body ###

## Check input
pathSlots <- which(names(opt) %in% c("vastDir","outDir","sampleTab","contrTab","specEventTable","vast2EnsemblGene"))
suppressWarnings(opt[pathSlots] <- sapply(opt[pathSlots], normalizePath))

if (!file.exists(opt$vastDir) && !opt$continue)  {stop("vast-tools folder not found at ", opt$vastDir)}
if (!file.exists(opt$sampleTab))                 {stop(opt$sampleTab, " not found")}
if (!file.exists(opt$contrTab))                  {stop(opt$contrTab, " not found")}
if (!is.null(opt$specEventTable) && !file.exists(opt$specEventTable))   {
    stop(opt$specEventTable, " not found")
}
if (!opt$noGO && is.null(opt$vast2EnsemblGene))   {
    stop("If --noGO is not set, vast2EnsemblGene is required to generate GO input files")
}
if (!opt$noGO && !file.exists(opt$vast2EnsemblGene))   {
    stop(opt$vast2EnsemblGene, " not found")
}
if (!(opt$minRepFrac > 0 & opt$minRepFrac <= 1))      {stop("--minRepFrac must be in ]0, 1]")}
if (!(opt$filter %in% c("DEFAULT","STRICT","LEGACY"))) {stop("--filter must be DEFAULT|STRICT|LEGACY")}

if (!opt$continue) {
    vastMain <- sort(dir(opt$vastDir, pattern="INCLUSION_LEVELS_FULL-[[:alnum:]]{4,6}(-[[:alnum:]]{3,4})?.tab.*",
                         full.names=T),
                     decreasing=T)[1]
    if (is.na(vastMain)) {stop("vast-tools main table not found in ", opt$vastDir)}
} else {
    vastMain <- sort(dir(opt$outDir, pattern="INCLUSION_LEVELS_FULL-[[:alnum:]]{4,6}(-[[:alnum:]]{3,4})?_clean.*.tab.*",
                         full.names=T),
                     decreasing=T)[1]
    if (length(vastMain) != 1) {
        stop("Flag --continue is set but can't find filtered vast-tools inclusion table in ", opt$outDir)
    }
    if (!file.exists(file.path(opt$outDir, "AStables.Robj"))) {
        stop("Flag --continue is set but can't find ", file.path(opt$outDir, "AStables.Robj"))
    }

}

sampleTab <- read.csv(opt$sampleTab, as.is=T)
if (!all(c("Sample","Type") %in% names(sampleTab))) {
    stop("Sample table must contain columns Sample, Type")
}

contrTab <- read.csv(opt$contrTab, as.is=T)
if (!all(c("Experimental","Control","File") %in% names(contrTab))) {
    stop("Contrast table must contain columns Experimental, Control, File")
}

uqContr <- unique(unlist(lapply(1:nrow(contrTab), FUN=function(x) {contrTab[x,1:2]})))
foundContr <- uqContr %in% sampleTab$Type
if (!all(foundContr)) {
    stop("Sample type(s) ", paste(uqContr[!foundContr], collapse=", "), " has/have no samples associated")
}
if (is.null(contrTab$Contrast)) {
    contrTab$Contrast <- paste(contrTab$Experimental, contrTab$Control, sep=":")
}

diffFiles <- dir(opt$vastDir)
contrTab$found <- contrTab$File %in% diffFiles
if (!all(contrTab$found) && !opt$continue) {
    stop("DIFF table(s) not found: ", paste(contrTab$File[!contrTab$found], collapse=", "))
}

cols <- strsplit(opt$colors, split=",")[[1]]
if (!all(cols %in% colors())) {
    stop("Color(s) not defined: ", paste(cols[!(cols %in% colors())], collapse=", "))
}
if (length(uqContr) > length(cols)) {cols <- rep(cols, length(uqContr) %/% length(cols) + 1)}
sampleTab$col <- sapply(sampleTab$Type, FUN=function(x) {cols[which(uqContr == x)]})


## Read data
cat("Reading ", basename(vastMain), "...\n", sep="")
vast <- read.delim(vastMain)
if (any(duplicated(vast$EVENT))) {
    vast <- removeVastDupEvents(vast)
}

if (opt$continue) {
    vast <- reshapeVast(vast)
}

sampleTab$vastCol <- sapply(sampleTab$Sample, FUN=function(x) {
    out <- which(names(vast) == x)
    if (length(out) != 1) out <- NA
    out
})
if (any(is.na(sampleTab$vastCol))) {
    stop("Sample(s) present in sample table but not found in vast-tools output:\n",
         paste(sampleTab$Sample[is.na(sampleTab$vastCol)], collapse=", "))
}


if (!is.null(opt$vast2EnsemblGene)) {
    v2g <- read.delim(opt$vast2EnsemblGene)
    if (!all(c("EVENT","EnsemblGeneID") %in% names(v2g))) {
        stop("vast2EnsemblGene table needs columns 'EVENT', 'EnsemblGeneID'")
    }
}

if (!dir.exists(opt$outDir)) {dir.create(opt$outDir, recursive=T)}

if (opt$continue) {
    cat("Loading AS tables in AStables.Robj...\n")
    load(file.path(opt$outDir, "AStables.Robj"))
}


## Save log file
contTag <- ifelse(opt$continue, "_plots", "")
logName <- file.path(opt$outDir, paste0("vastResultsProcessing", contTag, ".log"))
saveLog(logName, opt, sampleTab, contrTab)


if (!opt$continue) {
    
## Qual score filtering
    cat("Quality filtering...\n")

    if (opt$filter == "DEFAULT") {
        covThresCE    <- c("SOK","OK","LOW","VLOW")
       	covThresOther <- c("SOK","OK","LOW","VLOW")
        covThresIR    <- 10
	balThresCE    <- c("OK","B1","B2","Bl","Bn")
    }
    if (opt$filter == "STRICT") {
        covThresCE    <- c("SOK","OK","LOW")
       	covThresOther <- c("SOK","OK","LOW")
        covThresIR    <- 15
	balThresCE    <- c("OK","B1","B2","Bl","Bn")
    }
    if (opt$filter == "LEGACY") {
        covThresCE    <- c("SOK","OK","LOW")
       	covThresOther <- c("SOK","OK","LOW")
        covThresIR    <- 15
	balThresCE    <- c("OK","B1")
    }

    for (i in seq(8, ncol(vast), 2)) { # for all quality columns...
        vast[,i - 1] <- ifelse(cleanAS(vast[,i], 
                                       covThresCE=covThresCE, covThresOther=covThresOther, covThresIR=covThresIR,
				       balThresCE=balThresCE,
                                       complexCol=vast$COMPLEX), vast[,i-1], NA)
    }
    gz1 <- gzfile(file.path(opt$outDir, paste0(sub("\\.tab(\\.gz)?", "_clean", basename(vastMain)),
                                               substr(format(Sys.time(), "%Y%m%d"),3,8),
                                               ".tab.gz")),
                  "w")
    write.table(vast[,-seq(8, ncol(vast), 2)], file=gz1, row.names=F, col.names=T, quote=F, sep='\t')
    try(close(gz1))

    
    ## Extract point estimates from DIFF files
    cat("Extracting average PSI...\n")
    diffres <- mclapply(file.path(opt$vastDir, contrTab$File), read.delim, mc.cores=opt$cores)
    names(diffres) <- contrTab$Contrast
    diffres <- mclapply(diffres, FUN=mergeToVast, events=vast$EVENT, mc.cores=opt$cores)
    
    means <- sapply(uqContr, diffPointEst, 
       minFrac=opt$minRepFrac, contrTab=contrTab, diffres=diffres, sampleTab=sampleTab, vast=vast)
    cat("\n")
    rownames(means) <- vast$EVENT
    
    cat("Extracting differential events...\n")
    diffsig <- sapply(names(diffres), FUN=function(x) {diffres[[x]]$MV.dPsi._at_0.95})
    rownames(diffsig) <- vast$EVENT

    
    ## Combine Alt5/Alt3 subevents
    if (opt$keepFullAlt) {
        means.long <- means
        diffsig.long <- diffsig
        
        type <- character(nrow(vast))
        type[vast$COMPLEX %in% c("C1","C2","C3","S","ANN")] <- "CE"
        type[grep("IR-", vast$COMPLEX)] <- "IR"
        type[vast$COMPLEX == "MIC"]  <- "MIC"
        type[type == "CE" & vast$LENGTH <= 27]  <- "MIC"
        type[vast$COMPLEX == "Alt5"] <- "Alt5"
        type[vast$COMPLEX == "Alt3"] <- "Alt3"
        
        info.long <- data.frame(vast[,1:6], TYPE=type)
    }

    means.shr <- altCombinePSI(vast=vast, dat=means)
    diffsig.shr <- altCombine.sig(vast=vast, diffsig=diffsig)
    
    reorder <- orderEvents(info=means.shr[,1:7])
    info    <- means.shr[reorder, 1:7]
    means   <- as.matrix(means.shr[reorder, 8:ncol(means.shr)])
    diffsig <- as.matrix(diffsig.shr[reorder,])

    
    ## Generate dpsi and diff tables
    dpsi <- sapply(1:nrow(contrTab), dpsiMake, contrTab=contrTab, uqContr=uqContr, means=means)
    colnames(dpsi) <- contrTab$Contrast
    
    if (opt$strictDiff) {
        diff05 <- 100 * diffsig >  5 & !is.na(dpsi)
        diff10 <- 100 * diffsig > 10 & !is.na(dpsi)
        diff15 <- 100 * diffsig > 15 & !is.na(dpsi)
    } else {
        diff05 <- 100 * diffsig > 0 & abs(dpsi) >  5
        diff10 <- 100 * diffsig > 0 & abs(dpsi) > 10
        diff15 <- 100 * diffsig > 0 & abs(dpsi) > 15
    }
    
    if (opt$keepFullAlt) {
        dpsi.long <- sapply(1:nrow(contrTab), dpsiMake, contrTab=contrTab, uqContr=uqContr, means=means.long)
        colnames(dpsi.long) <- contrTab$Contrast
        
        if (opt$strictDiff) {
            diff05.long <- 100 * diffsig.long >  5 & !is.na(dpsi.long)
            diff10.long <- 100 * diffsig.long > 10 & !is.na(dpsi.long)
            diff15.long <- 100 * diffsig.long > 15 & !is.na(dpsi.long)
        } else {
            diff05.long <- 100 * diffsig.long > 0 & abs(dpsi.long) >  5
            diff10.long <- 100 * diffsig.long > 0 & abs(dpsi.long) > 10
            diff15.long <- 100 * diffsig.long > 0 & abs(dpsi.long) > 15
        }
    }

    ## Save tables
    gz1 <- gzfile(file.path(opt$outDir, "ASmeans.tab.gz"), "w")
    write.table(data.frame(info, means),
                file=gz1,
                row.names=F, col.names=T, quote=F, sep='\t')
    close(gz1)
    
    gz1 <- gzfile(file.path(opt$outDir, "dPSI.tab.gz"), "w")
    write.table(data.frame(info, dpsi),
                file=gz1,
                row.names=F, col.names=T, quote=F, sep='\t')
    close(gz1)
    
    save(info, means, dpsi, diffsig, diff05, diff10, diff15,
         file=file.path(opt$outDir, "AStables.Robj")
         )
    
    if (opt$keepFullAlt) {
        save(info.long, means.long, dpsi.long, diffsig.long, diff05.long, diff10.long, diff15.long,
             file=file.path(opt$outDir, "AStables_fullAltSS.Robj")
             )
    }
    
    out <- data.frame(info, round(dpsi,2), diff15)[which(apply(diff15, MAR=1, FUN=function(x) {any(x)})),]
    names(out)[8:(ncol(out))] <- paste(rep(c("dPSI","signif"), each=ncol(dpsi)),
                                       rep(names(out)[8:(7 + ncol(dpsi))], 2), sep=".")
    gz1 <- gzfile(file.path(opt$outDir, "AS.DiffEvents_dPSI.15.tab.gz"), "w")
    write.table(out, file=gz1, row.names=F, col.names=T, quote=F, sep='\t')
    close(gz1)
    
    out <- data.frame(info, round(dpsi,2), diff10)[which(apply(diff10, MAR=1, FUN=function(x) {any(x)})),]
    names(out)[8:(ncol(out))] <- paste(rep(c("dPSI","signif"), each=ncol(dpsi)),
                                       rep(names(out)[8:(7 + ncol(dpsi))], 2), sep=".")
    gz1 <- gzfile(file.path(opt$outDir, "AS.DiffEvents_dPSI.10.tab.gz"), "w")
    write.table(out, file=gz1, row.names=F, col.names=T, quote=F, sep='\t')
    close(gz1)
    
    out <- data.frame(info, round(dpsi,2), diff05)[which(apply(diff05, MAR=1, FUN=function(x) {any(x)})),]
    names(out)[8:(ncol(out))] <- paste(rep(c("dPSI","signif"), each=ncol(dpsi)),
                                       rep(names(out)[8:(7 + ncol(dpsi))], 2), sep=".")
    gz1 <- gzfile(file.path(opt$outDir, "AS.DiffEvents_dPSI.05.tab.gz"), "w")
    write.table(out, file=gz1, row.names=F, col.names=T, quote=F, sep='\t')
    close(gz1)
    
}


### Analyses ###

## Numbers of changing events
rowCol <- setRowCol(ncol(dpsi))

pdf(file.path(opt$outDir, "ChangingEvents.dPSI.15_bars.pdf"),
    wid=1.8 + rowCol[2] * 2.9, hei=1.6 + rowCol[1] * 2.7)
changingEventTypes(dpsi=dpsi, sig=diff15, info=info, includeANN=T, main=contrTab$Contrast)
dev.off()

pdf(file.path(opt$outDir, "ChangingEvents.dPSI.10_bars.pdf"),
    wid=1.8 + rowCol[2] * 2.9, hei=1.6 + rowCol[1] * 2.7)
changingEventTypes(dpsi=dpsi, sig=diff10, info=info, includeANN=T, main=contrTab$Contrast)
dev.off()

pdf(file.path(opt$outDir, "ChangingEvents.dPSI.05_bars.pdf"),
    wid=1.8 + rowCol[2] * 2.9, hei=1.6 + rowCol[1] * 2.7)
changingEventTypes(dpsi=dpsi, sig=diff05, info=info, includeANN=T, main=contrTab$Contrast)
dev.off()


## Clustering of single samples by correlation
suppressWarnings(
    use <- 
        apply(vast[,seq(7,ncol(vast) - 1,2)], MAR=1,
              FUN=function(x) {length(which(is.na(x)))}) < 1/3 * nrow(sampleTab) &
        apply(vast[,seq(7,ncol(vast) - 1,2)], MAR=1, max, na.rm=T) -
        apply(vast[,seq(7,ncol(vast) - 1,2)], MAR=1, min, na.rm=T) > 10
)

pdf(file.path(opt$outDir, "PSI.cluster.corr.single.pdf"),
    wid=5 + 0.1 * nrow(sampleTab), hei=5.2 + 0.1 * nrow(sampleTab))
cors <- cor(vast[use,seq(7,ncol(vast) - 1,2)], use="p")
heatmap.2(cors, trace="n", col=heat.colors, na.col="grey90", margins=c(12,12),
          cex.main=0.5, main=paste("Correlation of PSI\n          (events with PSI range > 10)"),
          key.title="", key.xlab="Corr. coefficient",
          labRow=sampleTab$Sample, labCol=sampleTab$Sample)
dev.off()


## Plot MDS
distSing <- dist(t(vast[use,seq(7,ncol(vast) - 1,2)]))
mds <- cmdscale(distSing, k=2)
xlim <- range(mds[,1]) + c(-0.08,0.08) * (max(mds[,1]) - min(mds[,1]))
ylim <- range(mds[,2]) + c(-0.08,0.08) * (max(mds[,2]) - min(mds[,2]))
pdf(file.path(opt$outDir, "MDS_samplePSI.pdf"), wid=5.5, hei=6)
plot(mds[,c(1,2)], pch=20, col=sampleTab$col, xlim=xlim, ylim=ylim,
     main="Multi-dimensional scaling", xlab="DIM 1", ylab="DIM 2")
text(mds[,c(1,2)], labels=sub("\\.R1$", "", rownames(mds)), cex=0.7, pos=1, col=sampleTab$col)
dev.off()


## Clustering of differential events
if (ncol(dpsi) > 1) {
    change <- which(apply(diff10, MAR=1, FUN=function(x) {any(x)}) &
                    apply(dpsi, MAR=1, FUN=function(x) {length(which(is.na(x)))}) < 1/4 * nrow(sampleTab))

    pdf(file.path(opt$outDir, "dPSI.cluster.corr.pdf"),
        wid=5 + 0.1 * nrow(contrTab), hei=5 + 0.1 * nrow(contrTab))
    cors <- cor(dpsi[change,], use="p")
    heatmap.2(cors, trace="n", col=ifelse(min(cors) < 0, cm.colors, heat.colors), na.col="grey90", margins=c(16,16),
              cex.main=0.4, main="Correlation of dPSI\n(changing events with dPSI >= 10 only)",
              denscol="grey40", key.title="", key.xlab="Corr. coefficient",
              labRow=contrTab$Contrast, labCol=contrTab$Contrast)
    dev.off()
}


## Clustering by dPSI euclidian distance
if (ncol(dpsi) > 1) {
    try({
        distTreat <- dist(t(dpsi[change,]))
        distEvent <- dist(dpsi[change,])
        clustTreat <- hclust(distTreat, method="ward.D2")
        clustEvent <- hclust(distEvent, method="ward.D2")

        clustTreat <- as.dendrogram(clustTreat)
        clustEvent <- as.dendrogram(clustEvent)
        clustEvent <- reorder(clustEvent, wts=rowMeans(dpsi[change,], na.rm=T))

        typecol <- ifelse(info$TYPE == "CE", "deepskyblue", NA)
        typecol[info$TYPE == "MIC"]  <- "darkblue"
        typecol[info$TYPE == "IR"] <- "orange"
        typecol[info$TYPE == "Alt5"] <- "darkolivegreen"
        typecol[info$TYPE == "Alt3"] <- "darkolivegreen2"

        dpsi.plot <- as.matrix(dpsi[change,-3])
        dpsi.plot[dpsi.plot < -30] <- -30
        dpsi.plot[dpsi.plot >  30] <-  30

        mycols <- c("cyan4", cm.colors(9), "magenta4")


        pdf(file.path(opt$outDir, "AS.cluster.dPSI.pdf"),
            wid=3 + 0.1 * nrow(contrTab), hei=11)
        heatmap.2(as.matrix(dpsi[change,]), mar=c(22,3), main="dPSI of changing events",
                  trace="n", col=mycols, na.color="grey80", labRow="", labCol=contrTab$Contrast,
                  RowSideColors=typecol[change], lwid=c(1,2,1), lhei=c(1,6),
                  cexCol=1.4, density.info="none", sepcolor=NA, sepwidth=c(0,0), key.title="", key.xlab="dPSI",
                  Rowv=clustEvent, Colv=clustTreat,
                  breaks=c(-100, seq(-50,50,length.out=10), 100))
        dev.off()

        pdf(file.path(opt$outDir, "AS.cluster.legend.dPSI.pdf"), wid=5, hei=5)
        ypos <- barplot(rep(1,5), col=unique(typecol)[5:1], horiz=T, xaxt="n")
        text(0.5, ypos, c("Retained intron","Alternative 5'-ss","Alternative 3'-ss", "Microexon", "Cassette exon"),
             col="white", cex=2.3)
        dev.off()
    }, silent=T)
    if (!exists("clustEvent")) {warning("Too many NAs, event clustering failed")}
}



## Clustering of each event type separately by euclidian distance
if (ncol(dpsi) > 1 & exists("clustEvent")) {
    pdf(file.path(opt$outDir, "AS.cluster_types_dPSI10.pdf"),
         wid=3.5 + 0.1 * nrow(contrTab), hei=11)
    for (i in unique(info$TYPE)) 
        try({
               evType <- info$TYPE %in% i
               use <- which(1:nrow(dpsi) %in% change & evType)
               distTreat <- dist(t(dpsi[use,]))
               distEvent <- dist(dpsi[use,])
               clustTreat <- hclust(distTreat, method="ward.D2")
               clustEvent <- hclust(distEvent, method="ward.D2")
    
		clustTreat <- as.dendrogram(clustTreat)
        	clustEvent <- as.dendrogram(clustEvent)

        	labRow <- ""
        	if (i == "MIC") {labRow <- as.character(info$GENE[use])}
    
		heatmap.2(as.matrix(dpsi[use,]), mar=c(18, 6),
                  main=paste("dPSI of", length(use), "\nchanging", i, "events"), cex.main=1,
                  trace="n", col=mycols, na.color="grey80",
                  labRow=labRow, labCol=contrTab$Contrast,
                  lhei=c(1,6),
                  cexCol=1.2, density.info="none", sepcolor=NA, sepwidth=c(0,0), key.title="", key.xlab="dPSI",
                  Rowv=clustEvent, Colv=clustTreat, 
                  breaks=c(-100, seq(-50,50,length.out=10), 100))
    })
    dev.off()
}


## Scatter plots
if (nrow(contrTab) > 1 & (nrow(contrTab) <= maxPlotContr | opt$scatterForce)) {
    pdf(file.path(opt$outDir, "dPSI.scatter.pdf"), wid=16.5, hei=3.7)
    for (i in 1:(nrow(contrTab) - 1)) {
        for (j in 2:nrow(contrTab)) {
            if (i >= j) next
            scatterXY(contrTab$Contrast[i], contrTab$Contrast[j])
        }
    }
    dev.off()
}
    

## Percentage of changing events specific to certain situations
if (!is.null(opt$specEventTable)) {
    specEv <- read.delim(opt$specEventTable)
    if (!("EVENT" %in% names(specEv))) {stop("Specific event table is lacking column 'EVENT'")}
    if (length(which(as.character(specEv$EVENT) %in% sub("-[0-9]+/[0-9]+", "", info$EVENT))) / nrow(specEv) < 0.5) {
        stop("Less than half of the events in ", opt$specEventTable, " found in vast-tools results")
    }

    specRowCol <- setRowCol(ncol(specEv) - 1)
    pdf(file.path(opt$outDir, "SpecificEvents.dPSI.10.pdf"),
        wid=0.5 + specRowCol[2] * 3, hei=0.5 + specRowCol[1] * 2.9)
    for (i in 1:ncol(dpsi)) {
        changingSpecEvents(specEv=specEv, dpsi.i=dpsi[,i], sig.i=diff10[,i], info, main=contrTab$Contrast[i])
    }
    dev.off()
}


## Save files for GO analysis
if (!opt$noGO) {
    if (!dir.exists(file.path(opt$outDir, "GO", "input"))) {
        dir.create(file.path(opt$outDir, "GO", "input"), recursive=T)
    }

    v2g <- cleanV2G(v2g, info)
    saveGOfiles(v2g, info, opt, contrTab, diff10, dpsi)
}


## Finish
write.table(paste("Completed", strftime(Sys.time())),
            row.names=F, col.names=F, quote=F, sep='\t',
            file=logName, append=T)


