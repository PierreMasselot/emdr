#' EMD-R1 preparation
#'
#' Prepare a data.frame to be used in a regression function for EMD-R1, 
#'    including the predictors MIMF and the response.
#'
#' @param x \code{mimf} object to be prepared for regression.
#' @param y Vector containing the response variable of the regression.
#' @param covariates List containing eventual covariates to be added to the 
#'    output data.frame.
#' @param tt Vector of custom time indices. Useful to lag irregular time series.
#' @param lag Vector of lagging values for the IMFs in \code{x}. If NA (the
#'    default) lags are automatically chosen by the function 
#'    \code{\link{choose_lag}}.
#' @param lag.max The maximum lag to consider in the automatic lag choice. See
#'    \code{\link{choose_lag}}.
#' @param fill.na Logical indicating if irregularly sampled IMFs should be
#'    filled with NAs.
#'
#' This function take a \code{mimf} object, a response vector and prepare them
#'    to be easily used in any regression function. It essentially transform
#'    the MIMF array into a matrix containing all IMFs and lags these
#'    IMFs if necessary.
#'
#' Also allows the inclusion of non-IMF predictors in \code{covariates} that
#'    are included in the output data.frame.
#'
#' @return A data.frame containing all necessary variables for a regression
#'    function. The first column contains the response vector \code{y}. The
#'    following columns contain the IMFs to be used and the last columns 
#'    contain the eventual covariates.
#'
#' @references
#'  Masselot, P., Chebana, F., Belanger, D., St-Hilaire, A., Abdous, B., 
#'    Gosselin, P., Ouarda, T.B.M.J., 2018. EMD-regression for modelling 
#'    multi-scale relationships, and application to weather-related 
#'    cardiovascular mortality. \emph{Science of The Total Environment} 
#'    612, 1018-1029. 
#'
#' @examples
#'    library(dlnm)
#'    
#'    # Predictor decomposition
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(123)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(10:11, 12:13), 
#'      new.names = c("C10", "C11"))
#'
#'    # Response variable
#'    Y <- chicagoNMMAPS$resp[attr(cmimfs, "tt")]
#'
#'    # Data preparation: includes the day-of-week variable as potential
#'    # confounder
#'    dataR1 <- pimf(cmimfs, Y, covariates = list(dow = 
#'      chicagoNMMAPS$dow[attr(cmimfs, "tt")]))
#'    
#'    # Apply the Lasso
#'    library(glmnet)
#'    lasso.res <- cv.glmnet(data.matrix(dataR1[,-1]), dataR1[,1], 
#'      family = "poisson")
#'
#'    # Compute sensitivity and plot results
#'    amps <- mean_amplitude(dataR1[,2:25])
#'    betas <- coef(lasso.res)
#'    s <- sensitivity(amps, coefs = betas[2:25]) 
#'    plot_emdr(matrix(s, ncol = 2, byrow = FALSE), periods = period(cmimfs), 
#'    show.coef = "nonzero", col = c("red", "blue"), pch = 16:17)
#'    abline(h = 0, lty = 2)
#'
#' @export
pimf <- function(x, y, covariates = NULL, tt = attr(x,"tt"), lag = NA, lag.max = 0.5, fill.na = FALSE)
{
    if (length(dim(y)) == 2) y <- y[,1]
    dims <- dim(x)
    if (length(dim(x))>2){
       dm <- dimnames(x)
       if (is.null(dm[[3]])) dm[[3]] <- 1:dims[3]
       dim(x) <- c(dims[1],prod(dims[2:3]))
       dimnames(x) <- list(dm[[1]],apply(expand.grid(dm[[2]],dm[[3]]), 1, paste, collapse="."))           
    }
    stopifnot(nrow(x) == length(y))
    x <- fill.mimf(x,tt)
    yf <- fill.mimf(y,tt)
    lagged <- array(NA,dim(x))
    lag <- rep_len(lag,ncol(x))
    for (i in which(!is.na(lag))){
        lagged[,i] <- c(rep(NA,lag[i]),x[1:(dims[1]-lag[i]),i])
    }
    if (sum(is.na(lag)) > 0){
       if (is.null(y)) stop("y missing")
       tolag <- x[,is.na(lag)]
       attr(tolag,"tt") <- attr(x,"tt")
       auto <- choose_lag(tolag, yf, lag.max = lag.max)
       lagged[,is.na(lag)] <- auto
       lag[is.na(lag)] <- attr(auto,"lag")
    }
    dimnames(lagged) <- dimnames(x)
    result <- data.frame(y = yf, lagged)
    if (!is.null(covariates)){
       covsf <- sapply(covariates, fill.mimf, tt)
       result <- cbind(result,covsf)
    }
    if (!fill.na){
       nas <- apply(result,1,function(x) any(is.na(x)))
       result <- result[!nas,]
       tt <- attr(x,"tt")[!nas]
    }
    attributes(result) <- c(attributes(result),list(lag = lag, tt = tt))
    return(result)
}

fill.mimf <- function(x, tt = attr(x,"tt"))
{
    steps <- diff(tt)
    tstep <- min(steps)
    where <- which((steps - tstep) > tstep)
    final <- array(NA, dim = c((max(tt)-(min(tt) - tstep))/tstep,max(c(dim(x)[2],1),na.rm=T), max(c(dim(x)[3],1),na.rm=T)), dimnames = dimnames(x))
    final[tt/tstep,,] <- x
    final <- drop(final)
    class(final) <- class(x)
    attr(final,"tt") <- seq(min(tt),max(tt),tstep)
    attr(final,"added.ind") <- setdiff(1:length(attr(final,"tt")),as.integer(tt/tstep))
    attr(final,"call") <- c(attr(x,"call"),match.call())
    attributes(final) <- c(attributes(final), attributes(x)[setdiff(names(attributes(x)),names(attributes(final)))])
    return(final)
}

#' Lag IMFs
#'
#' This function lags the all series in \code{x} according to their cross-
#'    correlation with \code{y}.
#'
#' @param x Array containing the IMFs to lag.
#' @param y vector containing the response with which to maximize the
#'    correlation.
#' @param lag.max The maximum lag to consider relatively to the mean period of
#'    each IMF.
#'
#' @details In this function, the lag is determined automatically for each IMF 
#'    in \code{x}. The lag chosen is the one maximizing the cross-correlation
#'    (in absolute value) between the IMF and \code{y}.
#'
#'    For each IMF in \code{x}, the maximum lag considered in the search 
#'    is lag.max * IMF period (rounded). For example, for an IMF of mean period
#'    10 and \code{lag.max = 0.5}, the maximum lag considered is 5.
#'
#' @return An array with the same dimensions as \code{x} containing the lagged
#'    IMFs. In addition contains an attribute \code{lag} storing the lags chosen
#'    for each IMF.
#'
#' @examples
#'    choose_lag(cos(3:22), cos(1:20))
#'    choose_lag(cos(2:21), cos(1:20))
#'
#' @export  
choose_lag <- function(x, y, lag.max = 0.5)
{  
  stopifnot(length(dim(y))==0)
  n <- length(y)
  xper <- period(x)
  xarr <- as.matrix(x)
  dims <- dim(xarr)
  if (length(dims) == 2){
     dims <- dim(xarr) <- c(dims,1)
  }  
  xlag <- array(NA, dims)
  mlag <- array(NA, dims[-1])
  inds <- expand.grid(lapply(dims[-1],function(x) 1:x))
  for (i in 1:nrow(inds)){
      xi <- xarr[,inds[i,1],inds[i,2]]
      peri <- xper[inds[i,1],inds[i,2]]
      if(is.na(peri) || is.infinite(peri)) peri <- dims[1]
      ccfi <- stats::ccf(y,xi,lag.max=lag.max*peri,plot=F, 
        na.action = stats::na.pass)
      lagi <- which.max(abs(ccfi$acf[ccfi$lag>=0])) - 1
      mlag[inds[i,1],inds[i,2]] <- lagi
      xlag[,inds[i,1],inds[i,2]] <- c(rep(NA,lagi),xarr[1:(dims[1]-lagi),inds[i,1],inds[i,2]])
  }
  attributes(xlag) <- attributes(x)
  attr(xlag,"lag") <- drop(mlag)
  return(xlag)    
}