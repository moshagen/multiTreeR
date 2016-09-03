
#' doMT: estimate mpt model using multiTree
#'
#' Estimate an MPT model using multiTree
#'
#' @param fEqn eqn file location.
#' @param data either a string pointing to a mdt file or a N*M data matrix or a N*M data frame.
#' @param catlabels A string vector of category labels (matching the heading of the data file), only required for data matrix.
#' @param restrictions = either a string pointing to a restrictions file or a list containing restrictions e.g. list("a=b", "a=.5").
#' @param fOut filename where outputfile should be written.
#' @param lambda lambda value for PD^lamda used in minimization and fit assessment.
#' @param numRep how many times parameter estimation should be repeated (to avoid local minima).
#' @param maxIter maximum number of iterations.
#' @param convergence convergence criterion.
#' @param CI size of confidence interval in \%, e.g. 95 for 95\%
#' @param bootstrap samples number of bootstrap samples. 0 = no bootstrap.
#' @param jacobian whether to compute the jacobian.
#' @param fia whether to compute the fisher information approximation of the minimum description length.
#' @return mt results object (list)
#' @examples
#' \dontrun{
#' mtres <- doMT('2HT.eqn', '2HT.mdt', restrictions=list("g = 0.5"))
#' mtres$paramEst
#' summary(mtres)
#' }
#' @export
#'
doMT  <- function(fEqn, data, catlabels = NULL,
                  restrictions = NULL, fOut = "",
                  lambda = 0.0, numRep=3.0, maxIter = 10000.0, convergence = 1.0E-10, CI=95.0,
                  bootstrapSamples = 0.0, jacobian = FALSE, fia=FALSE){


  # unit tests
  #   fEqn <- 'X:/00 projekte/multinomial modeling/multitree/multiTreeR/tmp.eqn'
  #   data <- 'X:/00 projekte/multinomial modeling/multitree/multiTreeR/tmp.mdt'
  #   catlabels <- NULL
  #   restrictions <- NULL
  #   fOut <- ""
  #   lambda <- 0
  #   numRep <- 3
  #   maxIter <- 10000.0
  #   convergence <- 1.0E-8
  #   CI <- 95
  #   bootstrapSamples <- 0
  #   jacobian <- F
  #   fia <- T


  # do some input validation
  mdtfile <- T
  if(!is.character(data)){
    mdtfile <- F
    if(is.matrix(data) || is.data.frame(data)){
      if(is.null(colnames(data)) && is.null(catlabels)){
        warning("no column names in data file and none provided as catlabels")
        return();
      }
      if(!is.null(colnames(data)) && is.null(catlabels)){
        catlabels <- colnames(data)
      }
      # must be R numeric (not R integer as returned by rbinom, rmultinom!)
      data <- matrix(as.numeric(unlist(data)),
                     nrow(data), ncol(data),
                     dimnames = dimnames(data))
    }else{
      warning("data file must be a matrix or a dataframe")
      return();
    }
  }else{
    if(!file.exists(data)){
      warning("mdt file does not exist")
      return();
    }
  }
  restrfile <- F
  if(is.null(restrictions)){
    restrfile <- T
    restrictions = "";
  }else if(is.character(restrictions)){
    restrfile <- T
    if(!file.exists(restrictions)){
      warning("restrictions file does not exist")
      return();
    }
  }else{
    restrictions <- unlist(restrictions)
  }
  bOut = FALSE
  if(fOut != "" & nchar(fOut) > 0){
    bOut = TRUE
  }

  # set up options
  dOptions <- c(lambda, numRep, maxIter, convergence, CI, bootstrapSamples)
  bOptions <- c(jacobian, fia)

  # initialize: choose constructor depending on provided arguments
  if(mdtfile){
    if(restrfile){
      mtR <- .jnew("multitree/start/StartmnR", fEqn, data, restrictions, bOut, fOut, dOptions, bOptions)
    }else{
      mtR <- .jnew("multitree/start/StartmnR", fEqn, data, .jarray(restrictions), bOut, fOut, dOptions, bOptions)
    }
  }else{
    if(restrfile){
      mtR <- .jnew("multitree/start/StartmnR", fEqn, .jarray(data,dispatch=T), catlabels, restrictions, bOut, fOut, dOptions, bOptions)
    }else{
      mtR <- .jnew("multitree/start/StartmnR", fEqn, .jarray(data,dispatch=T), catlabels, .jarray(restrictions), bOut, fOut, dOptions, bOptions)
    }
  }

  if(is.jnull(mtR) || is.null(mtR)){
    warning("ERROR. failed to init mt")
    return()
  }

  # do ana
  .jcall(mtR, "V", "runAna")

  # get results (call .jmethods(mtR) to see all public methods)
  paramEst <- .jcall(mtR, "[[D", "getParamEst", simplify = T)
  paramLabels <- .jcall(mtR, "[S", "getParamLabels", simplify = T)
  ll <- .jcall(mtR, "[D", "getLogLik", simplify = T)
  fit <- .jcall(mtR, "[D", "getFitStatistic", simplify = T)
  df <- .jcall(mtR, "[I", "getDf", simplify = T)
  p <- .jcall(mtR, "[D", "getPasymptotic", simplify = T)
  p[p == -1] <- NA  # zero df
  aic <- .jcall(mtR, "[D", "getAIC", simplify = T)
  bic <- .jcall(mtR, "[D", "getBIC", simplify = T)
  aicDelta <- .jcall(mtR, "[D", "getAICdelta", simplify = T)
  bicDelta <- .jcall(mtR, "[D", "getBICdelta", simplify = T)

  if(fia){
    fia <- .jcall(mtR, "[D", "getFIA", simplify = T)
    cfia <- .jcall(mtR, "[D", "getCfia", simplify = T)
  }else{
    fia <- cfia <- NULL
  }

  res <- list(paramEst=paramEst,
              paramSE=matrix(NA, nrow(paramEst), ncol(paramEst)),
              paramCI=matrix(NA, nrow(paramEst), ncol(paramEst)),
              logLik = ll, fit=fit, df=df, p.value=p, AIC = aic,
              AIC.delta = aicDelta, BIC = bic, BIC.delta = bicDelta,
              fisherInformation = NULL,
              FIA = fia, cFIA = cfia,
              jacobian = NULL)
  colnames(res$paramEst) <- colnames(res$paramSE) <- colnames(res$paramCI) <- paramLabels


  # SE/CI etc. might fail:
  try({

    res$fisherinf <- .jcall(mtR, "[[[D", "getFisherInfo", simplify = T)

    ############# SE error #####################
    # Error in validObject(.Object) :
    #   invalid class “jobjRef” object: invalid object for slot
    # "jobj" in class "jobjRef": got class "NULL", should be or extend class "externalptr"
    res$paramSE <- .jcall(mtR, "[[D", "getParamSE", simplify = T)
    res$paramSE[res$paramSE == -1] <- NA # constant parameters

    res$paramCI <- .jcall(mtR, "[[D", "getParamCI", simplify = T)
    res$paramCI[is.na(res$paramSE)] <- NA # constant parameters

    if(jacobian){
      res$jacobian <- .jcall(mtR, "[[[D", "getJacobian", simplify = T)
      colnames(res$jacobian) <- rownames(res$jacobian) <- paramLabels
    }
  }, silent = TRUE)
  if(all(is.na(res$paramSE))){
    warning("Standard errors could not be estimated.")
  }


  class(res) <- "multiTreeR"
  return(res)
}





