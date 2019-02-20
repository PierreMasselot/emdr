#' MIMFs
#'
#' \code{is.imf} tests if its argument is an IMF.
#'
#' \code{as.mimf} turns an array into a \code{mimf} object.
#'
#' @param x A vector for \code{is.imf} or an array for \code{as.mimf}.
#' @param complete Logical. If TRUE computes the sum of all IMFs in \code{x}
#'    as the recomposed original serie.
#'
#' @details
#'    \code{is.imf} checks if the vector \code{x} respects the mathematical
#'    definition of an IMF, i.e. if its difference between the number of extrema
#'    zero-crossing is at most 1. In addition, computes the relative amplitude
#'    of the mean enveloppe.
#'
#'    \code{as.mimf} takes an array and transofrm it in a R object \code{mimf}
#'    regardless of its mathematical properties. It is made by adding the 
#'    necessary attributes sch as \code{tt}. It is mostly a convenient function.
#'
#' @return 
#'    \code{is.imf} returns a logical value indicating if the \code{x} is an
#'    IMF. 
#'
#'    \code{as.mimf} returns a \code{mimf} object.
#'
#' @export
as.mimf <- function(x, complete = T){    
    if (inherits(x,c("numeric","factor"))) x <- as.matrix(x)
    if (inherits(x,"data.frame")) x <- data.matrix(x)
    X <- NULL
    if (complete) X <- apply(x,(1:length(dim(x)))[-2],sum)
    attributes(x) <- c(attributes(x),list(X = X, tt = 1:dim(x)[1], nb.iter = rep(NA,dim(x)[2]), call = NULL))
    class(x) <- c(class(x),"mimf")
    return(x)
}

#' @rdname as.mimf
#'
#' @export
is.imf <- function(x)
{
  menv <- compute.mean.enveloppe(x)
  cond1 <- with(menv, abs(nzc-nex))
  cond2 <- with(menv, median(sqrt(rowSums(meanenv^2)) / amp))
  out <- cond1 <= 1 && menv$nzc > 2
  attributes(out) <- list(C1 = cond1, C2 = cond2)
  class(out) <- c("is.imf",class(out)) 
  return(out)
}

#--- Print method for the IMF evaluation ---
print.is.imf <- function(obj)
# obj : object resulting from is.imf function
{    
    cat(sprintf("%s\n",obj))
    cat(sprintf("Difference between number of extrema and zero crossings is %i \n",attr(obj,"C1")))
    cat(sprintf("Relative mean enveloppe is %.0f%% of the amplitude \n",attr(obj,"C2")*100))
}

