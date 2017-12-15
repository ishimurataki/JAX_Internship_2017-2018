plot.mediation <- function(med, 
                           col="firebrick4", 
                           chrlen=mouse.chrlen, ...){
  
  names(med) = toupper(names(med))
  
  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(med))
  
  if (!("GMB" %in% names(med)))
    med$GMB <- gmb.coordinates(med$CHR, med$POS, chrlen=chrlen)
  
  # reorganize chr-lengths as in gmb.coordinates
  max.chr <- max(as.numeric(med$CHR[grep("[0-9]+", med$CHR)]))
  unique.chr <- levels(factor(factor(med$CHR, levels=c(1:max.chr, "X", "Y", "M"))))
  chrlen = chrlen[unique.chr]
  
  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6
  
  ### Create X-axis gmb values for mediators (e.g. proteins/RNA) genome positions
  gmb.mediation.extended = c(min(med$GMB)-30,med$GMB,max(med$GMB)+30) #Used to extend the x-axis of the plot a bit
  chrlen = c(0, cumsum(chrlen))
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)
  
  #### Create the plot
  par(font = 2, font.lab = 2, font.axis = 2, xaxs="i",las = 1, mar=c(3, 4, 3, 1) + 0.1)
  plot(gmb.mediation.extended, c(0,med$LOD,0), col = 0, ylim = c(0, max(med$LOD)*1.05),ylab = "Conditioned LOD", xaxt = "n", xlab = "", ...)
  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5)-1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4],
       border = "grey60", lty=3, col = rgb(0.96,0.96,0.96))
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  points(med$GMB, med$LOD, type = "p", pch=21,lty = 1,lwd=1, cex=0.9, col=col) #dark red points
  text(chrmid, 0.97 * usr[4], names(chrlen)[-1], cex=1)
}

gmb.coordinates <- function(chr, pos, chrlen = mouse.chrlen) {
  
  # chr must be numeric or X or Y or M
  chr <- as.character(chr)
  stopifnot(grepl("[0-9]+", chr) | chr %in% c("X", "Y", "M") )
  
  # length of all chromosomes must be given
  stopifnot(chr %in% names(chrlen))
  
  # pos must be finite, non-negative, numberic
  stopifnot(is.finite(pos) & is.numeric(pos) & pos>=0)
  
  # length of 'chr' vector equals length of 'pos' vector
  stopifnot(length(chr) == length(pos))
  
  # number of autosomes
  max.chr <- max(as.numeric(chr[grep("[0-9]+", chr)]))
  
  # all chromosomes, ordered 1..max.chr,X,Y,M
  unique.chr <- levels(factor(factor(chr, levels=c(1:max.chr, "X", "Y", "M"))))
  
  # chrlen ordered as unique.chr
  chrlen = chrlen[unique.chr]
  
  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6
  if (max(pos) > MAXPOS) pos <- pos / 10^6
  
  # cumulative chr. lengths (=shift)
  chrcumsum <- c(0, cumsum(chrlen))
  names(chrcumsum) <- c(names(chrlen), "End") # shift by 1
  
  gmb <- pos + chrcumsum[chr]
  return(gmb)
}

mediation.scan <- function(target, 
                           mediator, 
                           annotation, 
                           qtl.geno, 
                           covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff", 
                                    "lod-ratio"), 
                           verbose=TRUE) {
  
  # calculates log10-Likelihood of linear model y ~ 1 + X
  LL <- function(y, X) {
    -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
  }
  
  # check input
  stopifnot(NROW(target) == NROW(mediator))
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(qtl.geno) == NROW(target))
  stopifnot(is.null(covar) | NROW(target) == NROW(covar))
  stopifnot(!any(is.na(covar)))
  stopifnot(!any(is.na(qtl.geno)))
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(qtl.geno)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(c("CHR", "POS") %in% toupper(names(annotation)))
  method = match.arg(method)
  
  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- ncol(mediator) # number of points to scan
  if (is.null(covar)) covar <- cbind(rep(1, N)) # if no covariates, use just intercept
  LOD <- rep(NA, N) # prepare output
  
  if (method == "double-lod-diff") {
    no.na <- !is.na(target)
    LOD0 <- LL(target[no.na], cbind(covar, qtl.geno)[no.na,]) - LL(target[no.na], covar[no.na,])
  }
  
  # for-loop comparing M0: target~covar+mediator[,i] vs M1: target~covar+mediator[,i]+qtl.geno
  for (i in 1:N) {
    if (verbose & i %% 1000 == 0) print(i)
    no.na <- !is.na(target) & !is.na(mediator[,i])
    loglik0 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i]))
    loglik1 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i], qtl.geno[no.na,]))
    
    if (method == "ignore" | (method == "double-lod-diff" & all(no.na))) {
      # "double-lod-diff" for no missing observation is identical to "ignore"
      LOD[i] <- loglik1 - loglik0
    } else {
      loglik2 <- LL(target[no.na], covar[no.na,])
      loglik3 <- LL(target[no.na], cbind(covar[no.na,], qtl.geno[no.na,]))
      
      if (method == "lod-diff") {
        LOD[i] <- loglik3 - loglik2 - (loglik1-loglik0)
      } else if (method == "double-lod-diff") {
        LOD[i] <- LOD0 - (loglik3 - loglik2 - (loglik1-loglik0))
      } else if (method == "lod-ratio") {
        LOD[i] <- (10^loglik1-10^loglik0) / (10^loglik3 - 10^loglik2)
      }
    }
  }
  
  output <- annotation
  output$LOD <- LOD
  class(output) <- c("mediation", "data.frame")
  return(output)
}
