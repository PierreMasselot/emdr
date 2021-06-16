#' EMD-R2
#'
#' Performs EMD-R2 by using the given regression function in each MIMF.
#'
#' @param formula A formula describing the models to be fit for each MIMF.
#' @param mimf \code{mimf} object containing both the response and predictor
#'    IMFs.
#' @param covariates A list or data.frame containing eventual covariates
#'    which are not IMFs.
#' @param tt A vector of time indices for \code{mimf}. The default is to use
#'    the \code{tt} attribute.
#' @param reg.fun The regression function to perform either as a function or
#'    as a character. Should include \code{formula} and \code{data} arguments.
#' @param reg.args A list containing additional arguments to be passed to 
#'    \code{reg.fun}.
#' @param pimf.args A list of arguments to be passed to the function 
#'    \code{\link{pimf}} to prepare the IMFs for the regression.
#' @param xy.args A character vector of length two indicating the argument
#'    names corresponding to the predictors and response in the case 
#'    \code{reg.fun} does not include arguments \code{formula} and \code{data}.
#'    The default of \code{xy.args = c("x", "y")} should cover the majority
#'    of cases.
#'
#' @details \code{emdr2} wraps the chosen regression function (\code{reg.fun}),
#'    and applies it once per MIMF in \code{mimf}. Variables given in
#'    \code{covariates} are used with each IMF.
#'
#'    \code{reg.fun} can include either a couple \code{formula} and \code{data}
#'    arguments, either \code{x} and \code{y} arguments. Particular
#'    argument names for the response or predictor variables can be given
#'    in \code{xy.args}. 
#'
#' @return A \code{emdr2} object which is a list of length nimf containing all
#'    results of \code{reg.fun} on individual MIMFs.
#'
#' @seealso \code{\link{extract.emdr2}} to extract a particular element of the
#'    result from \code{reg.fun} for each submodel. \code{\link{coef.emdr2}}
#'    to specifically extract the coefficients. \code{\link{predict.emdr2}}
#'    to obtain predictions from the whole EMD-R2 model.
#'
#' @references
#'  Masselot, P., Chebana, F., Belanger, D., St-Hilaire, A., Abdous, B., 
#'    Gosselin, P., Ouarda, T.B.M.J., 2018. EMD-regression for modelling 
#'    multi-scale relationships, and application to weather-related 
#'    cardiovascular mortality. \emph{Science of The Total Environment} 
#'    612, 1018-1029. 
#'
#' @examples
#'  
#'    library(dlnm)
#'    library(glmnet)
#'
#'    dat <- chicagoNMMAPS[,c("death", "temp", "rhum")]
#'
#'    mimfs <- memd(dat)
#'    cmimfs <- combine.mimf(mimfs, list(12:13, 14:17, 18:19), 
#'      new.names = c("C12", "C13", "r"))
#'    amps <- mean_amplitude(cmimfs)
#'
#'    # EMD-R2 with glm
#'    lm.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs)
#'    betas.R2 <- coef(lm.R2)
#'    sensitivity.R2 <- sensitivity(amps[,-1], coefs = betas.R2[,-1])
#'    plot_emdr(sensitivity.R2, periods = period(cmimfs)[,-1],  
#'      col = c("red", "blue"), pch = 16:17)
#'    abline(h = 0, lty = 2)
#'    
#'    # EMD-R2 with lasso
#'    lasso.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs, 
#'      reg.fun = "cv.glmnet")
#'    betas.R2 <- coef(lasso.R2, s = "lambda.1se")
#'    sensitivity.R2 <- sensitivity(amps[,-1], coefs = betas.R2[,-1])
#'    plot_emdr(sensitivity.R2, periods = period(cmimfs)[,-1], 
#'      show.coef = "nonzero", pch = 16:17)
#'    abline(h = 0, lty = 2)
#'  
#'
#' @export
emdr2 <- function(formula, mimf, covariates = NULL, tt = attr(mimf,"tt"), reg.fun = "glm", reg.args = list(), pimf.args = list(), xy.args = c("x", "y"))
{                                                                        
    nimf <- dim(mimf)[2]
    vars <- all.vars(formula)
    if ("." %in% vars){
       vars <- vars[-which(vars==".")]
       toadd <- dimnames(mimf)[[3]][!dimnames(mimf)[[3]] %in% vars]
       vars <- c(vars, toadd)
    }
    reg.formals <- formals(reg.fun)
    result <- vector("list",nimf)
    for (i in 1:nimf){
        def.pimf <- list(x = mimf[,i,vars[-1]], y = mimf[,i,vars[1]], tt = tt, covariates = covariates)
        pimf.args <- c(def.pimf, pimf.args[!names(pimf.args) %in% names(def.pimf)])
        iprepared <- do.call(pimf, pimf.args)
        if ("formula" %in% names(reg.formals)){
           iform <- stats::as.formula(sprintf("y ~ %s", paste(names(iprepared)[-1], collapse = " + ")))
           def.args <- list(formula = iform, data = iprepared)
        } else {
           if (length(xy.args) != 2) stop(sprintf("'xy.args' must be provided since there is no 'formula' argument in '%s'", reg.fun))
           def.args <- list()
           def.args[[xy.args[1]]] <- as.matrix(iprepared[,-1])
           def.args[[xy.args[2]]] <- iprepared[,1]
        }
        reg.args <- c(def.args, reg.args[!names(reg.args) %in% names(def.args)])
        result[[i]] <- do.call(reg.fun,reg.args)
    }
    names(result) <- dimnames(mimf)[[2]]
    class(result) <- "emdr2"
    return(result)    
}

#' Extract results from an emdr2 object
#'
#' \code{extract.emdr2} extracts particular elements from an \code{emdr2} object 
#'    such as performance criteria or residuals.
#' \code{coef.emdr2} specifically extract regression coefficients.
#'
#' @param object An \code{emdr2} object.
#' @param what A character giving the name of the element to extract.
#' @param select A character of numeric vector giving a subset of MIMFs 
#'    for which to extract the element.
#'
#' @return \code{extract.emdr2}: a vector, matrix or list of the extracted 
#'    elements. The class of the output depend on the type of element.
#'
#' @seealso \code{\link{emdr2}} to produce an \code{emdr2} object.
#'
#' @examples
#'    library(dlnm)
#'
#'    dat <- chicagoNMMAPS[,c("death", "temp", "rhum")]
#'
#'    mimfs <- memd(dat)
#'    cmimfs <- combine.mimf(mimfs, list(12:13, 14:17, 18:19), 
#'      new.names = c("C12", "C13", "r"))
#'
#'    lm.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs)
#'    betas.R2 <- coef(lm.R2)
#'    aic.R2 <- extract.emdr2(lm.R2, what = "aic")
#'
#' @export  
extract.emdr2 <- function(object, what, select = NULL)
{
    if (missing(what)) stop("Argument 'what' missing")
    if (!(what %in% ls(object[[1]]))) stop("Invalid 'what' argument")
    if (is.null(select)) select <- 1:length(object)
    vals.list <- sapply(object[select],"[[",what)
    return(vals.list)
}

#' @rdname extract.emdr2
#'
#' @param method The coefficient method linked to the regression function
#'    used. The default should cover most cases.
#' @param ... Additional arguments to be passed to the \code{method}.
#'
#' @return \code{coef.emdr2}: a nimfs x nvariable matrix of coefficients.
#'
#' @export
coef.emdr2 <- function(object, method = "coef", ...)
{
    dots <- list(...)
    K <- length(object)
    coef.list <- vector("list",K)
    for (i in 1:K){
        met.args <- c(list(object[[i]]), dots)
        coef.list[[i]] <- drop(as.matrix(do.call(method, met.args)))
    }
    coef.mat <- do.call(rbind,coef.list)
    rownames(coef.mat) <- names(object)
    colnames(coef.mat) <- names(coef.list[[1]])
    return(coef.mat)
}


#' Predictions from emdr2
#'
#' Method predict for \code{emdr2} objects.
#'
#' @param object An \code{emdr2} object.
#' @param newmimfs Array of new MIMFs values for prediction. Number of IMFs and 
#'    variables must match those of the \code{mimf} argument given in the
#'    \code{\link{emdr2}} function. If NULL, the data used to fit the model
#'    are used.
#' @param newcovariates Data.frame of new covariates value for prediction.
#' @param method The name of a prediction method linked to the regression
#'    function used for \code{emdr2}.
#' @param ... Additional arguments for the predict method.
#'
#' @details If new data are given for prediction, both \code{newmimfs} and 
#'    \code{newcovariates} must be given (provided covariates have been used in 
#'    the model).
#'
#' @return A vector of predicted values.
#'
#' @examples 
#'    library(dlnm)
#'
#'    dat <- chicagoNMMAPS[,c("death", "temp", "rhum")]
#'
#'    mimfs <- memd(dat)
#'    cmimfs <- combine.mimf(mimfs, list(12:13, 14:17, 18:19), 
#'      new.names = c("C12", "C13", "r"))
#'
#'    lm.R2 <- emdr2(death ~ temp + rhum, mimf = cmimfs)
#'    deathHat <- predict(lm.R2)
#'
#' @export
predict.emdr2 <- function(object, newmimfs = NULL, newcovariates = NULL, method = "predict", ...)
{
    dots <- list(...)
    K <- length(object)
    pred.list <- vector("list",K)
    if (length(dim(newmimfs)) == 2) dim(newmimfs) <- c(dim(newmimfs),1)
    for (i in 1:K){
        met.args <- list(object[[i]])
        if (!is.null(newmimfs)) met.args <- c(met.args, list(cbind(newmimfs[,i,], newcovariates)))
        met.args <- c(met.args, dots)
        pred.list[[i]] <- drop(do.call(method, met.args))
    }
    yhat <- Reduce("+",pred.list)
    return(yhat)
}


