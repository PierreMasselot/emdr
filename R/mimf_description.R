#' Basic summary of an mimf object
#'
#' \code{period} returns the mean period and \code{mean_amplitude} returns the 
#'    mean amplitude of each imf on \code{object}. \code{summary}
#'    wraps them up into a single object.
#'
#' @param object A \code{mimf} object resulting from a call to \code{memd}. 
#' @param tt A vector of custom time indices.
#' @param x An object of class \code{summary.mimf} resulting from a call to
#'    \code{summary.mimf}.
#' @param ... Additional parameters to be passed to the function
#'    \code{\link[base]{format}} for the print method.
#'
#' @details \code{period} computes the mean period of an IMF by counting its 
#'    number of zero crossing. The period then corresponds to twice the mean
#'    time difference between two zero crossings.
#'
#'    \code{mean_amplitude} computes the root mean square amplitude of the IMF, 
#'    i.e. the square root of the sum of its squared values.
#'
#'    \code{summary} computes the mean periods, amplitudes and a reconstruction
#'    error, i.e. the mean relative difference between the original series
#'    and the sum of IMFs. The \code{print} method allows pretty printing
#'    in the console.
#'
#' @return \code{period} and \code{mean_amplitude} both return a \emph{nimfs} x 
#'    \emph{nvariables} matrix containing the mean period and mean amplitude
#'    respectively.
#'
#'    \code{summary} return a list containing the mean periods, mean
#'    amplitudes, and the relative difference between the original
#'    signal and its reconstruction when EEMD is used.
#'
#' @seealso \code{\link{as.mimf}} to convert an array into a \code{mimf} object.
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Decompose both temperature and relative humidity with NA-MEMD
#'    # Adding two noise variables 
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(3)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(11:12, 13:19, 20:21), 
#'      new.names = c("C11", "C12", "r"))
#'
#'    summary(cmimfs)
#'
#' @export
period <- function(object, tt = attr(object, "tt"))
{
    computePeriod <- function(x,tt){
        ind <- zero.crossings(x)
        zc <- tt[ind]
        return(2*mean(diff(zc)))
    }
    if (inherits(object,c("numeric","factor"))) object <- as.matrix(object)
    if (is.null(tt)) tt <- 1:nrow(object)
    nas <- apply(object,1,function(x)any(is.na(x)))
    if (length(dim(object))>2) object <- object[!nas,,,drop = F] else object <- object[!nas,,drop = F]
    tt <- tt[!nas]    
    periodArray <- as.matrix(apply(object,2:length(dim(object)),computePeriod,tt))
    return(periodArray)      
}

#' @rdname period
#'
#' @export
mean_amplitude <- function(object, tt = attr(object, "tt"))
{
    if (inherits(object,c("numeric","factor"))) object <- as.matrix(object)
    if (inherits(object,"data.frame")) object <- data.matrix(object)
    if (is.null(tt)) tt <- 1:dim(object)[1]
    nas <- apply(object,1,function(x)any(is.na(x)))
    if (length(dim(object))>2) object <- object[!nas,,] else object <- object[!nas,]
    tt <- tt[!nas]
    amps <- apply(object,2:length(dim(object)),function(x){mean(compute.mean.enveloppe(x,tt)$amp)})
    return(amps*2)
}

#' @rdname period
#'
#' @export
summary.mimf <- function(object, ...)
{
     X <- attr(object,"x")
     tt <- attr(object, "tt")
     imfperiods <- period(object, tt)
     imfperiods[nrow(imfperiods),] <- Inf
     imfamp <- mean_amplitude(object, tt)
     Xreconstruct <- apply(object,(1:(length(dim(object))))[-2],sum)
     recon.resid <- X - Xreconstruct
     relError <- apply(recon.resid^2,2,mean)/diag(stats::var(X))
     relError[relError<.Machine$double.eps] <- 0
     out <- list(periods = imfperiods, amplitude = imfamp, reconstruction.error = relError, call = attr(object,"call"))
     class(out) <- "summary.mimf"
     return(out)
}

#' @rdname period
#'
#' @export
print.summary.mimf <- function(x, ...)
{
     if (is.null(dimnames(x$periods)[[1]])) dimnames(x$periods)[[1]] <- sprintf("C%s",1:nrow(x$periods))
     if (is.null(dimnames(x$periods)[[2]])) dimnames(x$periods)[[2]] <- sprintf("Var%s",1:ncol(x$periods))
     cat("Call:\n")
     print(x$call)
     cat("\nIMF mean periods:\n")
     print(format(t(x$periods), ...), quote=F)
     cat("\nIMF mean amplitude:\n")
     print(format(t(x$amplitude), ...), quote=F)
     cat("\nRelative reconstruction error:\n")
     cat(paste(sprintf("%s %%",format(x$reconstruction.error, ...)),collapse="\t"))
     cat("\n")
} 

#' Plot IMFs
#'
#' Method to display the (M)IMFs obtained by the function \code{\link{memd}}.
#'
#' @param x \code{mimf} x or array to plot.
#' @param tt Vector containing custom time indices for the IMFs. If NULL, looks
#'    for the \code{tt} attribute of \code{x}.
#' @param select.var Character or numeric vector giving a subset of variables 
#'    for which to plot the IMFs. 
#' @param select.imf Character or numeric vector giving a subset of IMFs to 
#'    plot.
#' @param input Logical. If TRUE, the top panel shows the original signal.
#'    Only considered if \code{x} is a \code{mimf} object.
#' @param input.lab The label of the panel showing the input signal.
#' @param imf.lab Character vector giving labels for the IMFs. NULL displays the
#'    dimnames of \code{x} and NA removes the labels.
#' @param grid Character giving the type of grid to plot. "zeroline" 
#'    (the default) only draws the zeroline, "complete" draws a complete grid 
#'    and "none" draws no grid at all.
#' @param grid.col The grid color.
#' @param space Numeric value giving the margin between two panels in the plot.
#' @param ... Other graphical parameters. See 
#'    \code{\link[graphics]{par}}.
#' 
#' @details One panel is drawn for each IMF. In the multivariate case, 
#'    by default all IMF's variables are scaled and displayed on the same panel.
#'    To obtain the true amplitude of each variable, they must be plotted
#'    separately.
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Decompose both temperature and relative humidity with NA-MEMD
#'    # Adding two noise variables 
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'
#'    # Plot the two variables on the same graph
#'    plot(mimfs)
#'
#'    # Plot the two variables separately
#'    plot(mimfs, select.var = "temp", col = "red")
#'    plot(mimfs, select.var = "rhum", col = "blue")
#'
#' @export                                                 
plot.mimf <- function(x, tt = NULL, select.var = NULL, select.imf = NULL, input = TRUE, input.lab = "X", imf.lab = NULL, grid = c("zeroline","complete","none"), grid.col = "lightgray", space = 1, ...)
{
    grid <- match.arg(grid)
    stopifnot(inherits(x,c("matrix","array")))
    if (length(dim(x))==2){
       dm <- dimnames(x)
       dim(x) <- c(dim(x),1)
       dimnames(x) <- c(dm,list(NULL))
    } 
    if (is.null(select.var)) select.var <- 1:dim(x)[3]
    if (is.null(select.imf)) select.imf <- 1:dim(x)[2]
    n <- dim(x)[1]
    K <- length(select.imf)
    p <- length(select.var)
    if (is.null(tt)){
       tt <- if (inherits(x,"mimf")) attr(x,"tt") else 1:n
    }
    if (is.null(imf.lab)) imf.lab <- dimnames(x)[[2]][select.imf]
    imf.lab <- rep_len(imf.lab, K)
    if (inherits(x,"mimf") && input){      
       K <- K + 1
       imf.lab <- c(input.lab,imf.lab)
       mimfs <- array(NA,c(n,K,p),list(NULL,imf.lab,dimnames(x)[[3]][select.var]))
       mimfs[,1,] <- attr(x,"x")[,select.var]
       mimfs[,-1,] <- x[,select.imf,select.var]
    } else {
       mimfs <- x[,select.imf,select.var,drop=F]
       dimnames(mimfs)[[2]] <- imf.lab
    }
    dots <- list(...)
    checkLPar <- function(param){
         if (inherits(param,"matrix")){
            if (all(dim(param) == c(p,K))){
               outparam <- t(param)
            } else {
               if (all(dim(param) == c(K,p))){
                  outparam <- param
               } else {
                  stop(sprintf("Dimension of the matrix provided to parameter %s do not fit number of IMFs to plot",deparse(substitute(param)))) 
               }
            }            
         } else {
            if (length(param) == K) outparam <- matrix(param,K,p)
            if (length(param) == p) outparam <- matrix(param,K,p,byrow=T)
            if (!length(param) %in% c(p,K)) outparam <- matrix(rep_len(param,p),K,p,byrow=T)
         }
         return(outparam)     
    }
    argsMatplot <- dots[names(dots) %in% c("type","lty","lend","pch","col","cex","lwd")]
    argsMatplot <- lapply(argsMatplot,checkLPar)    
    argsPlot <- dots[names(dots) %in% setdiff(unique(names(c(formals(graphics::plot.default),formals(graphics::plot.xy)))),names(argsMatplot))]
    argsPar <- dots[names(dots) %in% setdiff(names(graphics::par()),c(names(argsPlot),names(argsMatplot)))]
    defPar <- list(mfrow = c(K,1), oma = c(max(5-space,0),4,4,2) + .1, mar = c(space,0,0,0), tcl = 0.5)
    argsPar <- c(argsPar, defPar[!names(defPar) %in% names(argsPar)])
    do.call(graphics::par,argsPar)
    for (k in 1:K){
        argsPlotK <- within(argsPlot,{
            y <- mimfs[,k,]
            if (p > 1){
               y <- scale(y)
            } 
            ylab <- dimnames(mimfs)[[2]][k]
            if (!exists("type")) type <- "l"
            x <- tt
            xpd <- NA
            xaxt <- "n"
            yaxt <- "n"
            if (k < K) xlab <- ""
            else if (!exists("xlab")) xlab <- "tt"
            main <- ""            
        })
        argsK <- lapply(argsMatplot, "[", i =k, j = )
        do.call(graphics::matplot,c(argsPlotK[!names(argsPlotK) %in% names(argsMatplot)],argsK))
        if ((!argsPlot$xaxt == "n") || is.null(argsPlot$xaxt)){
           if (inherits(tt, "Date")){
             graphics::axis.Date(1, x = tt, labels = FALSE)
             graphics::axis.Date(3, x = tt, labels = FALSE)
             if (k == K) graphics::axis.Date(1, x = tt, tcl = -0.5)
           } else {
             if (inherits(tt, "POSIXt")){
               graphics::axis.POSIXct(1, x = tt, labels = FALSE)
               graphics::axis.POSIXct(3, x = tt, labels = FALSE)
               if (k == K) graphics::axis.POSIXct(1, x = tt, tcl = -0.5)
             } else {
               graphics::axis(1, labels = FALSE)
               graphics::axis(3, labels = FALSE)
               if (k == K) graphics::axis(1, tcl = -0.5)
             }
           } 
           yaxp <- graphics::par("yaxp")
           if (yaxp[3] %% 2 == 1) yaxp[3] <- yaxp[3] - 1
           ytcks <- graphics::axTicks(2, axp = yaxp)
           if(graphics::strheight(0) > diff(ytcks[1:2])) ytcks <- graphics::axTicks(2,axp = c(graphics::par("yaxp")[1:2],2))
           if (p == 1) graphics::axis(2,at=ytcks) else graphics::axis(2,at=ytcks, labels=FALSE)   
        }
        if (grid == "complete") grid(col=grid.col)
        if (grid == "zeroline") graphics::abline(h=0,lty="dotted",col=grid.col)
    }
    if (is.null(dots$main) && !is.null(dimnames(mimfs)[[3]]) && p == 1) dots$main <- dimnames(mimfs)[[3]] 
    if (!is.null(dots$main)) graphics::title(main=dots$main,outer=TRUE)
}

#' Adds IMFs
#' 
#' Sums one or several IMFs and returns a \code{mimf} object. Useful for
#'    when several IMFs have similar frequencies.
#'
#' @param object \code{mimf} object or array.
#' @param select List of vector giving the IMFs to be summed. Each vector gives
#'    the IMFs to be summed into a single IMF. Elements of \code{select} should
#'    not overlap.
#' @param new.names Character vector giving names ti the new IMFS. If NULL
#'    default names are generated.
#'
#' @details When noise variables are added to the signal to perform the NA-MEMD,
#'    it often results in a too large number of IMFs with very similar 
#'    oscillation periods. \code{combine.mimf} allows gathering htese IMFs into 
#'    a single one. The IMFs to comnine should be carefully selected by 
#'    checking their mean period and examining the plot.
#'
#' @return
#'    Returns a \code{mimf} object.
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Decompose both temperature and relative humidity with NA-MEMD
#'    # Adding two noise variables 
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(3)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'
#'    # Plot resulting IMFs
#'    plot(mimfs)
#'
#'    # Sum MIMFs with similar frequencies
#'    cmimfs <- combine.mimf(mimfs, list(11:12, 13:19, 20:21), 
#'      new.names = c("C11", "C12", "r"))
#'    plot(cmimfs)
#' 
#' @export 
combine.mimf <- function(object, select, new.names = NULL)
{
    mimfs <- object   
    if(length(dim(mimfs))==2) mimfs <- array(mimfs,c(dim(mimfs),1))
    dims <- dim(mimfs)  
    K <- dims[2]
    if(any(unlist(select) > K)) stop("'select' contains invalid IMF number")
    if(length(unique(unlist(select)))!=length(unlist(select))) stop("Attempt to combine twice the same IMF")
    if(any(unlist(lapply(select,diff)) > 1)) stop("Attempting to combine non consecutive IMFs")
    keptimfs <- setdiff(1:K,unlist(select))
    Kk <- length(keptimfs)
    Kc <- length(select)
    newdims <- dims
    newdims[2] <- Kk + Kc
    newind <- sapply(select,min)
    if (Kc > 1) newind <- newind - c(0, cumsum(sapply(select[-Kc],length) - 1))
    keptnewind <- setdiff(1:newdims[2],newind)    
    new.mimfs <- array(NA,newdims)
    dimnames(new.mimfs)[c(1,3)] <- dimnames(mimfs)[c(1,3)]
    new.mimfs[,keptnewind,] <- mimfs[,keptimfs,]
    for (i in 1:Kc){
        nwimf <- apply(mimfs[,select[[i]],, drop=F],c(1,3),sum)
        new.mimfs[,newind[i],] <- nwimf
    }
    if (!is.null(dimnames(object)[[2]])){
     dimnames(new.mimfs)[[2]] <- rep("",newdims[2])
     dimnames(new.mimfs)[[2]][keptnewind] <- dimnames(object)[[2]][keptimfs]
    }
    if (is.null(new.names)){
      dimnames(new.mimfs)[[2]][newind] <- sapply(select,function(x) paste(dimnames(object)[[2]][x],collapse = "+"))
    } else {
      dimnames(new.mimfs)[[2]][newind] <- new.names
    }
    attributes(new.mimfs) <- c(attributes(new.mimfs),list(x = attr(object,"x"), tt = attr(object,"tt"), call = attr(object,"call")))
    class(new.mimfs) <- class(object)
    return(new.mimfs)
}