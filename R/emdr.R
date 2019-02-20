#' emdr: Empirical mode decomposition based regression
#'
#' The \code{emdr} package provides functions to decompose time series
#' and use the resulting components (called intrinsic mode functions, IMFs)
#' in a regression analysis.
#'
#' Functions in the package are roughly divided in three sections:
#'
#' \itemize{
#'  \item MEMD decomposition and description of the resulting IMFs;
#'  \item preparing
#'    MIMFs being used in a regression function to predict a non-decomposed
#'    response variables (denoted EMD-R1);
#'  \item  predicting one variable's
#'    IMF from other variables' IMFs (denoted EMD-R2).
#' }
#' 
#'
#' @section MEMD functions:
#' An analysis usually begins with a call to \code{\link{memd}} to decompose
#' a multivariate time-series into IMFs. This function supports different
#' modifications from the original algorithm, i.e. ensemble EMD and 
#' noise-assisted MEMD to address the issue of mode-mixing.
#' The result from \code{\link{memd}} is an object of class \code{mimf}
#' which can then be analyzed by looking at 
#' summarized characteristics through a call to 
#' \code{\link{summary.mimf}} or visually using the \code{\link{plot.mimf}}
#' method. Finally, \code{\link{imf.test}} performs an IMF significance test
#' and the method \code{\link{plot.imftest}} shows the result.
#'
#' @section EMD-R1:
#' The EMD-R1 design regresses a non-decomposed response against predictors'
#' IMFs. Thus, any regression function can be used to perform EMD-R1. The
#' function \code{\link{pimf}} prepares a \code{mimf} object as a 
#' \code{data.frame} to be used in a regression function. It includes lagging
#' the IMFs and adding non-IMF covariates. Since several IMFs can be correlated,
#' it is advised to consider the Lasso regression which can be performed by
#' the function \code{\link[glmnet]{glmnet}} in the package \code{glmnet}.
#' Resulting coefficients can then be standardized by the function 
#' \code{\link{sensitivity}} and displayed by the function 
#' \code{\link{plot_emdr}}.
#'
#' @section EMD-R2:
#' In the EMD-R2 design, the response variable is also decomposed and each of 
#' its IMFs is regressed against predictors' IMFs of similar frequencies.
#' After a call to \code{\link{memd}} to jointly decompose the response and
#' predictors, the resulting object can be used in the function 
#' \code{\link{emdr2}}. the result is a list of submodels for each IMF.
#' The function \code{\link{extract.emdr2}} extracts any element from each
#' submodel with the \code{\link{coef.emdr2}} method for coefficients
#' specifically. As for EMD-R1, these coefficients can then be standardized by 
#' the function \code{\link{sensitivity}} and displayed by the function 
#' \code{\link{plot_emdr}}.
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
#'  Masselot, P., Chebana, F., Belanger, D., St-Hilaire, A., Abdous, B., 
#'    Gosselin, P., Ouarda, T.B.M.J., 2018. EMD-regression for modelling 
#'    multi-scale relationships, and application to weather-related 
#'    cardiovascular mortality. \emph{Science of The Total Environment} 
#'    612, 1018-1029. 
#'
#' @docType package
#' @name emdr
NULL