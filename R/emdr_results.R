#' Response sensitivity to IMFs
#'
#' This function computes the sensitivity to each IMF, which is basically a
#'    regression coefficient scaled by the IMF mean amplitude.
#'
#' @param amplitudes A matrix of mean amplitudes as computed by 
#'    \code{\link{mean_amplitude}}. One column corresponds to a variable and
#'    one line to an IMF.
#' @param model The result from a regression function. Must have a 
#'    \code{\link[stats]{coef}} method. If not, use the argument \code{coefs}
#'    instead.
#' @param coefs User provided coefficient matrix. Must have the same dimensions
#'    as \code{amplitudes}.
#' @param ... Additional arguments for the \code{\link[stats]{coef}} method
#'    used by \code{model}.
#'
#' @details The sensitivy term hereby designates regression coefficients
#'    scaled according to the corresponding IMF's mean amplitude. It estimates
#'    the amplitude of response's variations explained by the IMF.
#'
#'    The function uses the results from a regression model to compute 
#'    the sensitivity.
#'    If the resulting object contains a \code{\link[stats]{coef}} method
#'    it is used to extract the necessary coefficients. If this is not the 
#'    case, the \code{coefs} argument must be used instead.
#'
#' @return A matrix of sensitivities, with the same dimensions as 
#'    \code{amplitudes}.
#'
#' @seealso \code{link{coef.emdr2}} to extract coefficients from an \code{emdr2}
#'    object.
#'
#' @references
#'  Masselot, P., Chebana, F., Bélanger, D., St-Hilaire, A., Abdous, B., 
#'    Gosselin, P., Ouarda, T.B.M.J., 2018. EMD-regression for modelling 
#'    multi-scale relationships, and application to weather-related 
#'    cardiovascular mortality. \emph{Science of The Total Environment} 
#'    612, 1018–1029. 
#'
#' @examples
#'    ## EMD-R1
#'    library(dlnm)
#'    library(glmnet)
#'    
#'    # Predictor decomposition
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(3)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(11:12, 13:19, 20:21), 
#'      new.names = c("C11", "C12", "r"))
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
#'    amps <- mean_amplitude(dataR1[,-1])
#'    betas <- coef(lasso.res)
#'    s <- sensitivity(amps, coefs = betas[-1]) 
#'
#'    ## EMD-R2
#'    dat <- chicagoNMMAPS[,c("death", "temp", "rhum")]
#'
#'    set.seed(3475)
#'    mimfs <- memd(dat, l = 2, wn.power = .1)
#'    cmimfs <- combine.mimf(mimfs, list(12:14, 15:19, 20:22), 
#'      new.names = c("C12", "C13", "r"))
#'
#'    # EMD-R2 with glm
#'    lm.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs)
#'    betas.R2 <- coef(lm.R2)
#'    amps <- mean_amplitude(cmimfs)
#'    sensitivity.R2 <- sensitivity(amps[,-1], coefs = betas.R2[,-1])
#'
#' @export
sensitivity <- function(amplitudes, model = NULL, coefs = NULL, ...)
{    
    if (is.null(coefs)){
       if (is.null(model)) stop("At least 'model' or 'coef' must be provided")
       coefs <- tryCatch(stats::coef(model, ...), error = function(w) "error")
       if (identical(coefs, "error")) stop("Method coef() must exist for the class of object 'model'. If not extract them manually and use the argument 'coef'")
       if (length(dim(coefs)) == 2){
          inds <- pmatch(colnames(amplitudes), colnames(coefs))
          coefs <- coefs[,inds]
       } else {
          inds <- pmatch(names(amplitudes), names(coefs))
          coefs <- coefs[inds]       
       }
       if (any(is.na(coefs))) warning("Some names of 'amplitudes' have not been found in the coefficients of 'model'")
    } else {
       if (length(amplitudes) != length(coefs)) stop("Inconsistent lengths between 'coefs' and 'amplitudes'")
    }
    return(amplitudes * coefs)
}


#' Plot coefficients of EMD-regression
#'
#' Plot the coefficients resulting from an EMD-regression acording to the
#'    mean period of the corresponding IMFs.
#'
#' @param x The coefficient matrix to plot. Can also contain sensitivities
#'    (see \code{\link{sensitivity}}).
#' @param periods Matrix containting the mean period of IMFs correspondind to
#'    the coefficients in \code{x}. See \code{\link{period}}. If NULL,
#'    the period are taken as the two to the power of the IMF's order.
#' @param lower,upper Matrices containing lower and upper confidence limits.
#' @param ci.args A list of arguments to be passed to the function
#'    \code{\link[graphics]{arrows}} for drawing confidence intervals.
#' @param period.log2 Logical. If TRUE, a log2 transformation is applied to the
#'    x axis.
#' @param trend.label The label to be displayed for the trend
#'    component's coefficient.
#' @param show.coef Character giving restrictions for the coefficients to draw.
#'    \code{show.coef = "all"} (the default) draws all coefficients, 
#'    \code{show.coef = "nonzero"} plots only nonzero coefficients, 
#'    \code{show.coef = "significant"} draws only coefficients for which the
#'    confidence interval given in \code{lower} and \code{lower} excludes 
#'    the zero value.
#' @param col,pch color and point type for drawn coefficients. If matrices one
#'    value corresponds to one coefficient and if vectors a value per variable
#'     is used.
#' @param line.pars List of parameters for the line separating the trend
#'    coefficients from the other. See \code{\link[graphics]{abline}}.
#' @param ... Other arguments to be passed to the plot. See 
#'    \code{\link[graphics]{par}}.
#' 
#' @examples
#'    ## EMD-R1
#'    library(dlnm)
#'    library(glmnet)
#'    
#'    # Predictor decomposition
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(3)
#'    mimfs <- memd(X, l = 2) # Takes a couple of minutes
#'    cmimfs <- combine.mimf(mimfs, list(11:12, 13:19, 20:21), 
#'      new.names = c("C11", "C12", "r"))
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
#'    amps <- mean_amplitude(dataR1[,2:27])
#'    betas <- coef(lasso.res)
#'    s <- sensitivity(amps, coefs = betas[2:27]) 
#'
#'    plot_emdr(matrix(s, ncol = 2, byrow = FALSE), periods = period(cmimfs), 
#'    show.coef = "nonzero", col = c("red", "blue"), pch = 16:17)
#'    abline(h = 0, lty = 2) 
#'
#'    ## EMD-R2
#'    dat <- chicagoNMMAPS[,c("death", "temp", "rhum")]
#'
#'    set.seed(3475)
#'    mimfs <- memd(dat, l = 2, wn.power = .1)
#'    cmimfs <- combine.mimf(mimfs, list(12:14, 15:19, 20:22), 
#'      new.names = c("C12", "C13", "r"))
#'
#'    # EMD-R2 with glm
#'    lm.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs)
#'    betas.R2 <- coef(lm.R2)
#'    amps <- mean_amplitude(cmimfs)
#'    sensitivity.R2 <- sensitivity(amps[,-1], coefs = betas.R2[,-1])
#'
#'    plot_emdr(sensitivity.R2, periods = period(cmimfs)[,-1],  
#'      col = c("red", "blue"), pch = 16:17)
#'    abline(h = 0, lty = 2)
#'
#' @export
plot_emdr <- function(x, periods = NULL, lower = NULL, upper = NULL, ci.args = list(), period.log2 = TRUE, trend.label = "Trend", show.coef = c("all","nonzero","significant"), col = NULL, pch = NULL, line.pars = list(), ...)
{
    dots <- list(...)
    if (is.null(dim(x))){
       K <- length(x)
       P <- 1
    }  else {
       K <- dim(x)[1]
       P <- dim(x)[2]
    }
    betas <- as.vector(x)
    if (is.null(periods)) periods <- rep(2^(1:K),P) else periods <- as.vector(periods)
    lower <- as.vector(lower)
    upper <- as.vector(upper)
    show.coef <- match.arg(show.coef)
    if (show.coef == "significant" && (is.null(lower) || is.null(upper))) show.coef <- "nonzero"
    nb <- length(betas)
    if (length(as.vector(periods)) != nb){
       periods <- rep_len(periods,nb)
       warning("periods length is different than number of coefficients.")    
    }
    trends <- is.na(periods) | is.infinite(periods)
    periods[trends] <- 2*max(periods,na.rm=T)
    if (period.log2) periods <- log2(periods)    
    plot.pars <- dots[names(dots) %in% c(names(formals(graphics::plot.default)), names(graphics::par(no.readonly = TRUE)))]
    plot.pars$x <- 0
    plot.pars$y <- 0
    plot.pars$col <- "white"
    plot.pars$xaxt <- "n"    
    plot.def <- list(ylim = range(c(betas,lower,upper)), xlim = range(periods), xlab = "Period", ylab = "Coefficient")
    plot.pars <- c(plot.pars, plot.def[!names(plot.def) %in% names(plot.pars)])
    do.call(graphics::plot,plot.pars)
    treat.args <- function(arg){
      if (is.null(arg)){
         arg <- rep(1:P,each=K)
      } else {
         if (is.matrix(arg)){
            arg <- as.vector(arg)
         } else {
            arg <- rep(rep_len(arg,P), each = K)
         }       
      }
    }
    tokeep <- switch(show.coef,
        all = rep(TRUE,nb),
        nonzero = betas!=0,
        significant = sign(lower)==sign(upper)
    )        
    if (!is.null(lower) && !is.null(upper)){
       arrows.def <- list(angle = 90, length = 0.05, code = 3)
       ci.args <- c(ci.args, arrows.def[!names(arrows.def) %in% names(ci.args)])
       ci.args[names(ci.args) %in% c("col","lty","lwd")] <- lapply(ci.args[names(ci.args) %in% c("col","lty","lwd")], treat.args)
       ci.args$x0 <- ci.args$x1 <- periods
       ci.args$y0 <- lower
       ci.args$y1 <- upper       
       ci.args[names(ci.args) %in% c("col","lty","lwd","x0","x1","y0","y1")] <- lapply(ci.args[names(ci.args) %in% c("col","lty","lwd","x0","x1","y0","y1")],"[",tokeep)
       do.call(graphics::arrows,ci.args)
    }
    pars.list <- list(x = periods ,y = betas, col = treat.args(col), pch = treat.args(pch))
    pars.list <- lapply(pars.list,"[",tokeep)    
    do.call(graphics::points,pars.list)    
    axis.pars <- plot.pars[grep("axis",names(plot.pars))]
    axis.pars$side <- 1
    axis.pars$at <- axis.pars$labels <- graphics::axTicks(1)
    if (period.log2) axis.pars$labels <- 2^axis.pars$labels
    if (any(trends[tokeep])){
       sortper <- sort(unique(periods))
       sep <- mean(sortper[length(sortper) + (-1:0)])
       line.pars$v <- sep
       do.call(graphics::abline,line.pars)
       axis.pars$at <- c(axis.pars$at[axis.pars$at < sep], 
         mean(c(sep,graphics::par("usr")[2]))) 
       axis.pars$labels <- c(axis.pars$labels[axis.pars$at < sep], 
         trend.label)
       pars.list$x[trends[tokeep]] <- seq(sep,graphics::par("usr")[2],length.out=sum(trends[tokeep])+2)[-c(1,sum(trends[tokeep])+2)]
    }
    do.call(graphics::axis,axis.pars)
}
