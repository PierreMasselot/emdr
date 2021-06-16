#' Multivariate empirical mode decomposition
#'
#' Decomposes a multivariate series through empirical mode decomposition (EMD).
#'    Includes the possibility to perform ensemble EMD (EEMD) as well as 
#'    noise-assisted multivariate EMD (NA-MEMD).
#'
#' @param x A numeric matrix containing the multivariate signal. Each column 
#'    corresponds to a variable.
#' @param tt A numeric vector the same size of \code{row(x)} containing time 
#'    indices for the data. Allows for irregularly sampled signals.
#' @param ndirections Integer. The number of projections necessary to compute 
#'    the multivariate enveloppe. Should be at least twice the number of
#'    variables (\code{ncol(x)}).
#' @param stopping Character indicating the stopping criterion of the sifting
#'    process. When \code{stopping = "absmean"} (the default) the criterion of 
#'    Rilling et al. (2003) based on the mean enveloppe is used. When 
#'    \code{stopping = "S"}, the stopping criterion of Huang et al. (2003)
#'    based on the number of iteration is used. See details.
#' @param tol A numeric vector givving the tolerance for the stopping criterion.
#'    A vector of length 3 when \code{stopping = "absmean"} and a single value
#'    when \code{stopping = "S"}. See details.
#' @param max.iter Integer giving the maximum number of iterations for the 
#'    sifting process.
#' @param max.mimfs Integer giving the maximum number of IMFs to extract. If 
#'    NULL (the default), IMFs are extracted until one extremum is left in
#'    the signal.
#' @param l Integer giving the number of gaussian white noise variables to add 
#'    for the noise-assisted MEMD.
#' @param Ne Integer giving the ensemble number for EEMD.
#' @param wn.power Numeric value > 0 giving the relative standard deviation of 
#'    noise variables (see details).
#' @param keep.noise Logical. If set to TRUE, the white noise channels created
#'    by the algorithm are returned as well. Note that it is useful only
#'    when \code{l > 0}.
#' @param memd.stop Numeric value. The algorithm stops and consider that the
#'     IMF has been reached when the number of extrema of the remaining
#'      signal is below \code{memd.stop}. Note that in the multivariate case,
#'       the number of extrema of all the projections must below this value.         
#'
#' @details The EMD algorithm iteratively estimates IMF beginning with the 
#'    highest frequency to the lowest frequency. Each time an IMF is estimated,
#'    it is retrieved from the signal and the next IMF is estimated. The
#'    algorithm continues until the remaining signal contains only one (or none)
#'    extremum. This remaining signal is then considered as the trend.
#'
#'    Each IMF is estimated through a sifting process, which consists in 
#'    fitting the envelopes of the signal on its local extrema. The local mean
#'    is then computed as the mean of the envelopes and is retrieved from the
#'    signal to yield an IMF. If this IMF is not satisfying enough, the same 
#'    process is repeated, until the stopping criterion is met. See references
#'    below for the full details of the algorithm.
#'
#'    The stopping criterion, given in \code{stopping}, is arguably the most 
#'    important parameter of the algorithm. Two criteria are currently 
#'    implemented. The one of Rilling et al. (2003) (\code{stopping = "absmean"})
#'    stops the sifting process when the relative mean envelope is below a 
#'    predetermined threshold. In this case, the parameter \code{tol} must be a
#'    vector of length 3. \code{tol[1]} gives the threshold  and \code{tol[2]} 
#'    gives the maximum proportion of the signal allowed to be above this 
#'    threshold. \code{tol[3]} indicates another threshold that none of the
#'    mean enveloppe is allowed to exceed. 
#'    The second criteria implemented is the S of Huang et al. (2003)
#'    (\code{stopping = "S"}). It stops the sifting process when the difference
#'    between the number of local extrema and zero-crossings is at most 1, for
#'    \code{tol} steps. In this case, \code{tol} is a single integer value.
#'
#'    The function also includes the possibility to perform noise-assisted 
#'    extensions to manage the mode-mixing issue. Ensemble EMD (EEMD) is
#'    performed by setting the number of ensembles \code{Ne} higher than 1.
#'    NA-MEMD is performed by setting the number of added white noise variables
#'    \code{l} higher than 0. When both are given, the NA-MEMD is performed.
#'    In all cases, the argument \code{wn.power} indicates the amplitude of
#'    added white noise.
#'
#'    At the end, IMFs with similar frequencies can be obtained. They can be
#'    summed using the function \code{\link{combine.mimf}}.
#'
#' @return An object of class \code{mimf}, i.e an array with dimension 
#'    \emph{nindividuals} x \emph{nimfs} x \emph{nvariables} which also 
#'    contains the attributes:
#'      \item{x}{The original multivariate signal.}
#'      \item{tt}{The vector of time indices.}
#'      \item{nb.iter}{A vector containing the number of sifting
#'        iterations necessary to estimate each IMF. If \code{Ne > 1},
#'        a \code{Ne} x nimfs matrix.}
#'      \item{call}{The function call.}
#'
#' @references 
#'    Huang, N.E., Shen, Z., Long, S.R., Wu, M.C., Shih, H.H., Zheng, Q., 
#'      Yen, N.-C., Tung, C.C., Liu, H.H., 1998. The empirical mode 
#'      decomposition and the Hilbert spectrum for nonlinear and non-stationary 
#'      time series analysis. \emph{Proceedings of the Royal Society of London. 
#'      Series A: Mathematical, Physical and Engineering Sciences} 454, 903-995.
#'
#'    Rehman, N., Mandic, D.P., 2010. Multivariate empirical mode decomposition.
#'      \emph{Proceedings of the Royal Society A: Mathematical, Physical and 
#'      Engineering Science} 466, 1291-1302.
#' 
#'    Rehman, N.U., Park, C., Huang, N.E., Mandic, D.P., 2013. EMD Via MEMD: 
#'      Multivariate Noise-Aided Computation of Standard EMD. 
#'      \emph{Advances in Adaptive Data Analysis} 05, 1350007.
#'
#' @examples
#'    library(dlnm)
#'    
#'    X <- chicagoNMMAPS[,c("temp", "rhum")]
#'    set.seed(3)
#'    
#'    # EMD
#'    imfs <- memd(X[,1])
#'    
#'    # EEMD
#'    imfs <- memd(X[,1], Ne = 100, wn.power = .05)
#'    
#'    # MEMD
#'    imfs <- memd(X) # Takes a couple of minutes
#'    
#'    # NA-MEMD
#'    imfs <- memd(X, l = 2, wn.power = .02) # Takes a couple of minutes
#'    
#'    # Plot resulting (M)IMFs
#'    plot(imfs)
#'
#' @export
memd <- function(x, tt = 1:NROW(x), ndirections = 64, 
  stopping = c("absmean", "S"), tol = c(0.075, 0.075, 0.75), max.iter = 50, 
  max.mimfs = NULL, l = 0, Ne = 1, wn.power = 0.02, keep.noise = FALSE, 
  memd.stop = ifelse(p + l > 1, 3, 2))
{
    if (is.data.frame(x)) x <- data.matrix(x)
    if (is.vector(x)) x <- matrix(x, ncol = 1)
    stopifnot(is.matrix(x))
    stopping <- match.arg(stopping)
    cc <- stats::complete.cases(x)
    x <- x[cc,,drop=F]
    tt <- tt[cc]
    n <- nrow(x)
    p <- ncol(x)
    # Standardization of variables
    Xsc <- scale(x)
    nimfs <- ifelse(is.null(max.mimfs), ceiling(2 * log2(n)), max.mimfs)
    mimfs <- array(NA, dim = c(n, nimfs, p+l, Ne))
    imfcount <- vector("numeric", Ne)
    siftcounts <- matrix(NA, Ne, nimfs)
    for (e in 1:Ne){
        # Addition of white noise variables 
        if (l > 0){
          wn <- matrix(stats::rnorm(n * l, 0, wn.power), ncol = l)
          Xwn <- cbind(Xsc, wn)
        } else {
          if (Ne > 1){
            wn <- matrix(stats::rnorm(n * p, 0, wn.power), ncol = p)
            Xwn <- Xsc + wn
          } else {
            Xwn <- Xsc
          }
        }
        #--- MEMD decomposition -------  
        #Hammersley sequence for projections
        if (p + l > 1){
           projmat <- get.projections(ndirections, p + l)
        } else {
           projmat <- matrix(1)
        }
        R <- Xwn
        imfcount[e] <- 0
        repeat {
              tostop <- compute.mean.enveloppe(R, tt, projmat) 
              if(all(tostop$nex < memd.stop)) break     
              if(!is.null(max.mimfs)) if(imfcount == max.mimfs) break 
              siftp <- sifting(R, tt, projmat, stopping, tol, max.iter)
              imfcount[e] <- imfcount[e] + 1
              if (imfcount[e] == nimfs){
                mimfs <- abind::abind(mimfs[,-nimfs,,, drop = F], 
                  array(NA ,dim=c(n, nimfs, p+l, Ne)), 
                  mimfs[,nimfs,,, drop = F],
                  along = 2)
                siftcounts <- cbind(siftcounts, matrix(NA, Ne, nimfs))
                nimfs <- 2 * nimfs
              }
              mimfs[,imfcount[e],,e] <- siftp$mimf
              siftcounts[e, imfcount[e]] <- siftp$count
              R <- R - siftp$mimf
        }
        mimfs[,nimfs,,e] <- R 
        if (Ne > 1){
           print(sprintf("Ensemble %i / %i",e,Ne))
           utils::flush.console()
        }
    }
    # Aggregating of all mimfs found
    nmimfs <- max(imfcount)
    mimfs <- mimfs[,c(1:nmimfs, nimfs),,, drop = F]
    if (!keep.noise){
      mimfs <- mimfs[,,1:p,, drop = F]
      varnames <- colnames(x)
    } else {
      if (is.null(colnames(x))){
        varnames <- sprintf("V%i", 1:p)
      } else {
        varnames <- colnames(x)
      }
      varnames <- c(varnames, sprintf("noise%i", 1:l))
    }
    finalmimfs <- apply(mimfs,c(1,2,3), mean, na.rm = T)
    for (j in 1:p){
      finalmimfs[,,j] <- finalmimfs[,,j] * attr(Xsc,"scaled:scale")[j]
    }
    finalmimfs[,nmimfs + 1,1:p] <- finalmimfs[,nmimfs + 1,1:p] + 
      matrix(attr(Xsc,"scaled:center"), n, p, byrow = T)
    dimnames(finalmimfs) <- list(NULL, 
      c(sprintf("C%s", 1:(dim(finalmimfs)[2]-1)), "r"),
      varnames)
    finalmimfs <- drop(finalmimfs)
    siftcounts <- siftcounts[,1:nmimfs]
    if (Ne == 1) siftcounts <- unlist(siftcounts)
    attributes(finalmimfs) <- c(attributes(finalmimfs),
      list(x = x, tt = tt, nb.iter = siftcounts, call = match.call()))
    class(finalmimfs) <- c("mimf","array")
    return(finalmimfs)
}

#' Internal functions for MEMD
#'
#' \code{get.projections} creates projection directions through Hammersley 
#'    sequence.
#' \code{sifting} performs the sifting process to obtain one MIMF.
#' \code{compute.mean.enveloppe} computes the mean enveloppe for sifting.
#'
#' @param ndirections Integer giving the number of projections.
#' @param p Integer giving the number of variables to decompose.
#'
#' @return \code{get.projections}: a \code{ndirection} x \code{p} matrix giving 
#'    the projection vectors.
get.projections <- function(ndirections, p){
    hammseq <- spatstat.geom::Hammersley(ndirections,rev(randtoolbox::get.primes(p-1)),raw=T)
    b <- 2*hammseq-1
    tht <- atan2(sqrt(as.matrix(apply(b[,1:(p-1),drop=F]^2,1,cumsum))[,(p-1):1,drop=F]),b[,p:2,drop=F])
    dir.vec <- cbind(rep(1,ndirections),matrix(apply(tht,1,function(x)cumprod(sin(x))),ncol=p-1))
    dir.vec[,1:(p-1)] <- cos(tht) * dir.vec[,1:(p-1)]
    return(dir.vec)
}

#' @rdname get.projections
#'
#' @param x A numeric matrix containing the multivariate signal. Each column 
#'    corresponds to a variable.
#' @param tt A numeric vector the same size of \code{row(x)} containing time 
#'    indices for the data. Allows for irregularly sampled signals.
#' @param projmat A \code{ndirection} x \code{p} giving the
#'    projection vectors to compute the enveloppe.
#' @param stopping Character indicating the stopping criterion of the sifting
#'    process. When \code{stopping = "absmean"} (the default) the criterion of 
#'    Rilling et al. (2003) based on the mean enveloppe is used. When 
#'    \code{stopping = "S"}, the stopping criterion of Huang et al. (2003)
#'    based on the number of iteration is used.
#' @param tol A numeric vector givving the tolerance for the stopping criterion.
#'    A vector of length 3 when \code{stopping = "absmean"} and a single value
#'    when \code{stopping = "S"}.
#' @param max.iter Integer giving the maximum number of iterations for the 
#'    sifting process.
#'
#' @return \code{sifting}: a list with the following elements:
#'  \item{mimf}{A matrix containing the resulting mimf.} 
#'  \item{count}{The number of iterations for completing the sifting process.} 
sifting <- function(x, tt, projmat, stopping, tol, max.iter){
    n <- nrow(x)
    hk <- x
    count <- 0
    if (stopping == "S") stpcount <- 0 #for the S stopping criterion
    repeat {
       menv <- compute.mean.enveloppe(hk, tt, projmat)       
       if (stopping == "absmean"){
          sx <- sqrt(rowSums(menv$meanenv^2)) / menv$amp
          if (!(mean(sx > tol[1]) > tol[2] || any(sx > tol[3]))) break
       }
       if (stopping == "S"){
          if (any(abs(menv$nex - menv$nzc) <= 1)){
             stpcount <- stpcount + 1
          } else {
             stpcount <- 0
          }
          if (stpcount >= tol) break
       } 
       if (count == max.iter) break
       hk <- hk - menv$meanenv
       count <- count + 1   
    }
    return(list(mimf = hk, count = count)) 
}

#' @rdname get.projections
#'
#' @param hk Matrix containing the current MIMF prototype. 
#'
#' @return \code{compute.mean.enveloppe}: a list with elements:
#'    \item{meanenv}{The mean enveloppe.}
#'    \item{amp}{The mean amplitude of the signal.}
#'    \item{nzc}{The number of zero crossings.}
#'    \item{nex}{The number of local extrema.}
compute.mean.enveloppe <- function(hk, tt = 1:n, projmat = diag(p)){
       hk <- as.matrix(hk)
       n <- nrow(hk)
       p <- ncol(hk)
       # projection along all direction vectors
       hkproj <- hk %*% t(projmat)
       # symmetry to manage boundary issues (Rilling et al. 2003)
       hkext <- rbind(hk[n:2,,drop=F],hk,hk[n:2-1,,drop=F])
       hkprojext <- rbind(hkproj[n:2,,drop=F],hkproj,hkproj[n:2-1,,drop=F]) 
       ttext <- c(tt[1] - rev(cumsum(diff(tt))), tt, tt[n] + cumsum(rev(diff(tt))))
       meanenv <- matrix(0,nrow=n,ncol=ncol(hk))
       amp <- vector("numeric",n)
       nzc <- vector("numeric",ncol(hkproj))
       nex <- vector("numeric",ncol(hkproj))
       for (i in 1:ncol(hkproj)){           
           # computes enveloppes
           extseq <- find.extrema(hkprojext[,i])
           nex[i] <- length(extseq$indmin[extseq$indmin >= n & extseq$indmin < 2*n])
           nzc[i] <- length(zero.crossings(hkproj[,i]))
           minindex <- extseq$indmin
           maxindex <- extseq$indmax
           if (length(minindex) > 0 && length(maxindex) > 0){
              minenv <- apply(hkext,2,function(x){stats::spline(ttext[minindex],x[minindex],xout = ttext)$y})
              maxenv <- apply(hkext,2,function(x){stats::spline(ttext[maxindex],x[maxindex],xout = ttext)$y})
              # computes local mean
              meanenv <- meanenv + ((maxenv[n:(2*n-1),] + minenv[n:(2*n-1),]) /2)
              amp <- amp + (sqrt(rowSums((maxenv-minenv)^2))[n:(2*n-1)] / 2)
           } else {
              meanenv <- meanenv + hk
              amp <- amp + 0
           }
       }
       amp <- amp / ncol(hkproj)
       meanenv <- meanenv / ncol(hkproj)
       return(list(meanenv = meanenv, amp = amp, nzc = nzc, nex = nex))
}

#' Find local extrema and zero-crossings of a data series
#' 
#' \code{find.extrema} gives the local extrema of \code{x} and 
#'    \code{zero.crossings} give the zero-crossings.
#'
#' @param x The signal.
#'
#' @return For \code{find.extrema}, a list giving the indices of local minima, 
#'    local maxima and the total number of extrema.
find.extrema <- function(x){
    n <- length(x)
    # when consecutive points have the same value, keep only one (the one on the middle)
    adoub <- which(diff(x) != 0)
    inddoub <- which(diff(adoub) != 1) + 1
    d <- adoub[inddoub] - adoub[inddoub - 1]
    adoub[inddoub] <- adoub[inddoub] - floor(d/2)
    adoub <- c(adoub,n)
    x1 <- x[adoub]
    # Search the extrema among the remaining points
    d2x <- diff(sign(diff(x1)))
    indmin <- adoub[which(d2x > 0) + 1] 
    indmax <- adoub[which(d2x < 0) + 1]
    nextr <- c(length(indmin),length(indmax))
    names(nextr) <- c("nmin","nmax")
    return(list(indmin = indmin, indmax = indmax, nextrema = nextr))
}

#' @rdname find.extrema
#'
#' @return For \code{zero.crossings}, a vector giving the time indices following 
#'    each zero crossing.
zero.crossings <- function(x){
    sx <- sign(x)
    indzc <- which(abs(diff(sx)) > 1) + 1
    #Check if there are (possibly consecutive) zero values
    if (any(sx==0)){
      indz <- which(sx==0)
      if (any(diff(indz)==1)){
        dz <- diff(c(0,sx==0,0))
        stz <- which(dz==1)
        endz <- which(dz==-1) - 1
        indz <- round((stz+endz)/2)
      }
    } else {
      indz <- NULL
    }
    return(sort(c(indzc,indz)))
}