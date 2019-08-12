#' IMF significance test
#'
#' Tests the hypothesis that IMFs do not contain more information than white
#'    noise.
#'
#' @param object \code{mimf} object or array containing (potentially
#'    multivariate) IMFs.
#' @param tt Numeric vector of custom time indices.
#' @param type The type of significance test to perform. Either the test of 
#'    Wu and Huang (2004) (\code{type = "wu"}) or the test of Flandrin et al. 
#'    (2004) (\code{type = "flandrin"}). 
#' @param alternative Character giving if the alternative hypothesis is
#'    one-sided (\code{alternative = "one.sided"}) or two-sided
#'    (\code{alternative = "two.sided"}).
#' @param H Numeric value giving the Hurst exponent necessary for the test of
#'    Flandrin et al. (2004). Must be between 0 and 1.
#'
#' @details
#'    In both tests, the first (highest frequency) IMF is used to estimate the
#'    noise level of the input signal and is thus never significant.
#'
#' @return
#'    An object of class \code{imftest} containing the following elements:
#'      \item{object}{The input \code{mimf} object name.}
#'      \item{alternative}{The alternative hypothesis.}
#'      \item{type}{The type of test performed.}
#'      \item{logp}{Matrix containing the log mean-period of each IMF.}
#'      \item{loge}{Matrix containing the log energy (variance) of each IMF.}
#'      \item{p.value}{Matrix containing the p-value of the test for each IMF.}
#'
#' @references
#'    Wu, Z., Huang, N.E., 2004. A study of the characteristics of white noise 
#'    using the empirical mode decomposition method. \emph{Proceedings of the 
#'    Royal Society of London. Series A: Mathematical, Physical and Engineering 
#'    Sciences} 460, 1597-1611.
#'
#'    Flandrin, P., Goncalves, P., Rilling, G.G., 2004. Detrending and 
#'    denoising with empirical mode decompositions. Presented at the \emph{12th 
#'    European Signal Processing Conference}, IEEE, Vienna, Austria.
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Decompose both temperature and relative humidity with NA-MEMD
#'    # Adding two noise variables 
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(123)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(10:11, 12:13), 
#'      new.names = c("C10", "C11"))
#'    
#'    # Apply the test of Wu and Huang (2004)
#'    testres <- imf.test(cmimfs)
#'    testres
#'    plot(testres)
#'
#' @export
imf.test <- function(object, tt = NULL, type = c("wu","flandrin"), alternative = c("one.sided","two.sided"), H = 0.5)
{
    stopifnot(inherits(object,c("array","matrix","numeric","integer")))
    type <- match.arg(type)
    if (inherits(object,c("numeric","integer"))){
       mimfs <- array(object,c(length(object),1,1), dimnames = c(names(object), list(NULL,NULL)))
    } else {
       mimfs <- object
    }
    if(length(dim(mimfs))==2) mimfs <- array(mimfs,c(dim(mimfs),1), dimnames = c(dimnames(object), list(NULL)))
    N <- dim(mimfs)[1]
    K <- dim(mimfs)[2]
    P <- dim(mimfs)[3]
    if (is.null(tt)){
       tt <- if (inherits(object,"mimf")) attr(object,"tt") else 1:N
    }
    imfener <- apply(mimfs,c(2,3),function(x)sum(x^2)/length(x))
    if (type == "wu"){
      alternative <- match.arg(alternative)
      imfener <- imfener / N
      for (j in 1:P) mimfs[,,j] <- mimfs[,,j]/sum(imfener[,j])# rescale
      imfener <- apply(mimfs,c(2,3),function(x)sum(x^2)/length(x))
      imfener <- 0.5636*imfener / matrix(imfener[1,],ncol=P,nrow=K,byrow=T)
      loge <- log(imfener)
      imfperiods <- period(object=mimfs,tt=tt)
      R <- is.na(imfperiods)
      logp <- log(imfperiods)
      p.value <- matrix(NA,nrow=K,ncol=P,dimnames=dimnames(mimfs)[2:3])
      p.value[1,] <- rep(1, P)
      for (j in 1:P){
          p.value[R[,j],j] <- 0
          for (k in which(!R[,j])[-1]){
              lneVec <- seq(0,-3-logp[k,j]^2,length.out=5000)
              dlne <- lneVec[2] - lneVec[1]
              lnePDF <- energy.spread(lneVec,0.12-0.934*logp[k,j],N)
              lneCDF <- cumsum(lnePDF)/sum(lnePDF)
              cdfind <- which.min(abs(lneVec-loge[k,j]))
              p.value[k,j] <- switch(alternative,
                  one.sided = lneCDF[cdfind],
                  two.sided = 2 * min(lneCDF[cdfind],1-lneCDF[cdfind])
              )           
          }
      }
      logp <- logp/log(2)
      loge <- loge/log(2)
    } else {
      alternative <- "one.sided"
      logp <- matrix(1:K,nrow = K, ncol = P, byrow = FALSE)
      loge <- log2(imfener)
      coefs <- switch(sign(H-0.5) + 2, c(0.487,0.458,-2.435,0.452,-1.951), c(0.719, 0.474, -2.449, 0.460, -1.919), c(1.025, 0.497, -2.331, 0.495, -1.833))
      H <- switch(sign(H-0.5) + 2, 0.2, 0.5, 0.8)
      Wk <- sapply(imfener[1,], function(x) (x/coefs[1])*(2.01+0.2*(H-0.5)+0.12*(H-0.5)^2)^(-2*(1-H)*(1:K)))
      Wk[1,] <- imfener[1,]
      u95 <- log2(Wk) + matrix(2^(coefs[2]*(1:K) + coefs[3]), nrow = K, ncol = P, byrow=FALSE)
      u99 <- log2(Wk) + matrix(2^(coefs[4]*(1:K) + coefs[5]), nrow = K, ncol = P, byrow=FALSE)
      p.value <- matrix(1, nrow = K, ncol = P)
      p.value[loge > u95] <- 0.025
      p.value[loge > u99] <- 0.005
    }
    dimnames(logp) <- dimnames(mimfs)[2:3]
    dimnames(loge) <- dimnames(mimfs)[2:3]
    cc <- stats::complete.cases(logp)
    out <- list(object = as.list(match.call())[[2]], alternative = alternative, type = type, logp = logp[cc,], loge = loge[cc,], p.value = p.value[cc,])
    if (type == "flandrin") out <- c(out, list(H = H))
    class(out) <- "imftest"
    return(out)    
}

#--- Subfunction of imf.test to compute the energy spread of the test of Wu and Huang (2004) ---
energy.spread <- function(y, ybar, N)
# y: vector, the signal for which to compute the energy
# ybar: vector, signal mean
# N: the number of observations

# Value: A vector giving the energy spread
{
    E <- exp(y)
    Ebar <- exp(ybar)
    rho <- 0.5*N*Ebar*log(N)-(E/Ebar-y)*N*Ebar/2
    rho <- rho - max(rho)
    return(exp(rho)) 
}

#' @rdname imf.test
#'
#' @param x An \code{imftest} object resulting from a call to \code{imf.test}.
#' @param digits The number of digits to print p-values. 
#' @param ... Formatting arguments to be passed to the function 
#'    \code{\link[base]{format}}
#'
#' @export
print.imftest <- function(x, digits = 3, ...)
{
    signif.star <- with(x,{
        mat <- matrix("",nrow=nrow(p.value),ncol=ncol(p.value))
        mat[p.value < .1] <- "."
        mat[p.value < .05] <- "*"
        mat[p.value < .01] <- "**"
        mat[p.value < .001] <- "***"
        mat
    })
    pval.string <- formatC(x$p.value, digits = digits, format="f", ...)
    pval.string[x$p.value < 0.001] <- "<0.001"
    printmat <- matrix(paste(pval.string,signif.star,sep=" "), nrow = nrow(x$p.value), dimnames = dimnames(x$p.value))
    if (is.null(dimnames(printmat)[[1]])) dimnames(printmat)[[1]] <- sprintf("C%s",1:nrow(printmat))
    if (is.null(dimnames(printmat)[[2]])) dimnames(printmat)[[2]] <- sprintf("Var%s",1:ncol(printmat))
    cat(sprintf("IMF significance test: '%s'\n\n",as.character(x$x)))
    cat(sprintf("Type: %s\n",x$type))
    cat(sprintf("Alternative: %s\n",x$alternative))
    print(printmat,quote=F)
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#' Graphical IMF significance test
#'
#' Graphically displays the result of an IMF significance test by plotting
#'    the log energy vs. log period of IMFs.
#'
#' @param x An \code{imftest} x resulting from a call to 
#'    \code{link{imf.test}}.
#' @param level The significance level at which to draw the rejectance region
#'    lines.
#' @param select.var Character of numeric vector giving a subset of variables 
#'    to plot.
#' @param select.imf Character of numeric vector giving a subset of IMFs 
#'    to plot.
#' @param alternative The alternative hypothesis, either "one.sided" or 
#'    "two.sided".
#' @param nline Numeric value giving the resolution of the rejectance region 
#'    lines.
#' @param col The colors of the points.
#' @param col.signif The color to use for significant IMFs exclusively. If NULL,
#'    \code{col} is used.
#' @param pch The plotting character of IMFs.
#' @param pch.signif Specific \code{pch} for significant IMFs only.
#' @param imf.names Character vector of IMF labels to plot next to the points.
#' @param col.names Color of the IMF labels.
#' @param col.lines Colors of the significance lines.
#' @param lty.lines Significance lines type.
#' @param lwd.lines Significance lines width.
#' @param add Logical value. If TRUE, the plot is added to an existing one.
#' @param ... Additional graphical parameters. See \code{\link[graphics]{par}}.
#'
#' @references
#'    Wu, Z., Huang, N.E., 2004. A study of the characteristics of white noise 
#'    using the empirical mode decomposition method. \emph{Proceedings of the 
#'    Royal Society of London. Series A: Mathematical, Physical and Engineering 
#'    Sciences} 460, 1597–1611.
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Decompose both temperature and relative humidity with NA-MEMD
#'    # Adding two noise variables 
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(123)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(10:11, 12:13), 
#'      new.names = c("C10", "C11"))
#'    
#'    # Apply the test of Wu and Huang (2004)
#'    testres <- imf.test(cmimfs)
#'    plot(testres)
#'
#' @export
plot.imftest <- function(x, level = 0.05, select.var = NULL, select.imf = NULL, alternative = x$alternative, nline = 500, col = "black", col.signif = NULL, pch = 4, pch.signif = 8, imf.names = NULL, col.names = NULL, col.lines = "black", lty.lines = 2, lwd.lines = 1, add = FALSE, ...)
{
    if (!alternative %in% c("one.sided","two.sided")){
       alternative <- "one.sided"
       warning("Unknown input for 'alternative'. 'one.sided' is used as a default")
    }
    col.lines <- rep_len(col.lines,2+(alternative=="two.sided"))
    lty.lines <- rep_len(lty.lines,2+(alternative=="two.sided"))
    lwd.lines <- rep_len(lwd.lines,2+(alternative=="two.sided"))
    logp <- x$logp
    loge <- x$loge
    p.value <- x$p.value
    dots <- list(...)
    argsPlot <- dots[names(dots) %in% unique(names(c(formals(graphics::plot.default),graphics::par())))]
    argsText <- dots[names(dots) %in% unique(names(formals(graphics::text.default)))]
    dims <- dim(x$p.value)
    K <- dims[1]
    p <- dims[2]
    N <- length(attr(get(as.character(x$object),pos=1),"tt"))
    if (x$type == "wu"){
      sigline <- matrix(NA,nrow=nline,ncol=2+(alternative=="two.sided" || x$type=="flandrin"))
      sigline[,1] <- seq(1,ceiling(max(logp, na.rm=T)),length.out=nline)
      for (i in 1:nline){
          lneVec <- seq(0,-3-sigline[i,1]^2,length.out=5000)
          dlne <- lneVec[2] - lneVec[1]
          lnePDF <- energy.spread(lneVec,-sigline[i,1],N)
          lneCDF <- cumsum(lnePDF)/sum(lnePDF)
          if (alternative == "one.sided"){
             idSig <- which(lneCDF>level)[1]-1
             sigline[i,2] <- lneVec[idSig] + dlne*(lneCDF[idSig+1]-level)/(lneCDF[idSig+1] - lneCDF[idSig])
          } else {
             idSig <- which(lneCDF>(level/2))[1]-1
             sigline[i,2] <- lneVec[idSig] + dlne*(lneCDF[idSig+1]-level)/(lneCDF[idSig+1] - lneCDF[idSig])
             idSig <- which(lneCDF>(1-(level/2)))[1]-1
             sigline[i,3] <- lneVec[idSig] + dlne*(lneCDF[idSig+1]-level)/(lneCDF[idSig+1] - lneCDF[idSig])
          }        
      }
      sigline[,-1] <- sigline[,-1] + matrix(0.066 * sigline[,1] + 0.12,nrow=nline,ncol=1+(alternative=="two.sided"))
      sigline <- sigline/log(2)
    }
    if (p > 1){
      if (is.numeric(select.var)) select.var <- select.var[select.var <= p]
      if (is.character(select.var)) select.var <- grep(pattern=select.var,colnames(logp))
      if (is.null(select.var)) select.var <- 1:p 
    } else {
      select.var <- 1
    }
    pv <- length(select.var)
    if (pv==0) stop("Incorrect variable specification")
    if (K > 1){
      if (is.numeric(select.imf)){
         if (any(select.imf > K)) warning(sprintf("select.imf contains values greater than the number of IMFs (%s). Values removed.", K))
         select.imf <- select.imf[select.imf <= K]
      } 
      if (is.character(select.imf) && K>1) select.imf <- grep(pattern=select.imf,rownames(logp))
      if (is.null(select.imf)) select.imf <- 1:K 
    } else {
      select.imf <- 1
    }
    Ks <- length(select.imf)
    argsPlot <- within(argsPlot,{
       if (!exists("ylab")) ylab <- "logE"
       if (!exists("xlab")) xlab <- "logT"
    })
    argsText <- within(argsText,{
       if (!exists("pos")) pos <- 3
    })
    if (is.null(imf.names)){
       if (!is.null(rownames(p.value))) imf.names <- rownames(p.value) else imf.names <- select.imf 
    }
    if (pv > 1) graphics::par(ask=TRUE)
    for (j in select.var){
        if (x$type == "flandrin"){
           sigline <- matrix(NA,nrow=K,ncol=3)
           sigline[,1] <- 1:K
           H <- x$H
           coefs <- switch(sign(H-0.5) + 2, c(0.487,0.458,-2.435,0.452,-1.951), c(0.719, 0.474, -2.449, 0.460, -1.919), c(1.025, 0.497, -2.331, 0.495, -1.833))
           Wk <- ((2^loge[1,j])/coefs[1])*(2.01+0.2*(H-0.5)+0.12*(H-0.5)^2)^(-2*(1-H)*(sigline[,1]))
           Wk[1] <- 2^loge[1,j]
           sigline[,2] <- log2(Wk) + 2^(coefs[2]*(sigline[,1]) + coefs[3])
           sigline[,3] <- log2(Wk) + 2^(coefs[4]*(sigline[,1]) + coefs[5])
        }
        colj <- col
        if (!is.null(col.signif)){
           colj <- rep(colj,Ks)
           colj[p.value[,j] < level] <- col.signif
        }
        pchj <- pch
        if (!is.null(pch.signif)){
           pchj <- rep(pchj,Ks)
           pchj[p.value[,j] < level] <- pch.signif
        }
        argsPlotj <- within(argsPlot,{
            x <- logp[select.imf,j]
            y <- loge[select.imf,j]
            col <- colj
            pch <- pchj
            if(!exists("main") && !is.null(colnames(p.value))) main <- colnames(p.value)[j]
            if(!exists("ylim")) ylim <- range(c(loge,sigline[,-1]))
            if(!exists("xlim")) xlim <- range(logp)
        })
        argsTextj <- within(argsText,{
            x <- logp[select.imf,j]
            y <- loge[select.imf,j]
            labels <- imf.names
            col <- if(!is.null(col.names)) col.names else colj 
        })
        if (add && (grDevices::dev.cur() > 1)){
           do.call(graphics::points,argsPlotj)
        } else {
           do.call(graphics::plot,argsPlotj)
        }
        do.call(graphics::text,argsTextj)
        graphics::matlines(sigline[,1],sigline[,-1], col = col.lines[-1], lwd = lwd.lines[-1], lty = lty.lines[-1])
        if (x$type == "wu"){
           graphics::abline(a = 0.12, b = -0.934, col = col.lines[1], lwd = lwd.lines[1], lty = lty.lines[1])
        } else {
           graphics::lines(1:K,log2(Wk), col = col.lines[1], lwd = lwd.lines[1], lty = lty.lines[1])
        }
    }
}
