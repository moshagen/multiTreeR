
#' doMT: estimate mpt model using multitree
#' @param feqn eqn file location
#' @param data either a string pointing to a mdt file or a N*M data matrix or a N*M data frame
#' @param catlabels A string vector of category labels (matching the heading of the data file), only required for data matrix
#' @param restrictions = either a string pointing to a restrictions file or a list containing restrictions e.g. list("a=b", "a=.5")
#' @param fOut filename where outputfile should be written
#' @param lambda lambda value for PD^lamda used in minimization and fit assessment
#' @param numRep how many times parameter estimation should be repeated (to avoid local minima)
#' @param maxIter maximum number of iterations
#' @param convergence convergence criterion
#' @param CI size of confidence interval in %, e.g. 95 for 95%
#' @param bootstrap samples number of bootstrap samples; 0 = no bootstrap
#' @param jacobian whether to compute the jacobian
#' @param fia whether to compute the fisher information approximation of the minimum description length
#' @return mt results object (list)
#' @examples
#' \dontrun{
#' mtres <- doMT('2HT.eqn', '2HT.mdt', restrictions=list("g = 0.5"))
#' mtres$paramEstm
#' summary(mtres)
#' }
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
      if(is.data.frame(data)){
        data <- as.matrix(data)
      }
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
  res.paramEst <- .jcall(mtR, "[[D", "getParamEst", simplify = T)
  res.paramEstSE <- .jcall(mtR, "[[D", "getParamSE", simplify = T)
  res.paramEstSE[res.paramEstSE == -1] <- NA # constant parameters
  res.paramEstCI <- .jcall(mtR, "[[D", "getParamCI", simplify = T)
  res.paramEstCI[is.na(res.paramEstSE)] <- NA # constant parameters
  res.ll <- .jcall(mtR, "[D", "getLogLik", simplify = T)
  res.fit <- .jcall(mtR, "[D", "getFitStatistic", simplify = T)
  res.df <- .jcall(mtR, "[I", "getDf", simplify = T)
  res.p <- .jcall(mtR, "[D", "getPasymptotic", simplify = T)
  res.p[res.p == -1] <- NA  # zero df
  res.aic <- .jcall(mtR, "[D", "getAIC", simplify = T)
  res.bic <- .jcall(mtR, "[D", "getBIC", simplify = T)
  res.aicDelta <- .jcall(mtR, "[D", "getAICdelta", simplify = T)
  res.bicDelta <- .jcall(mtR, "[D", "getBICdelta", simplify = T)
  res.fisherinf <- .jcall(mtR, "[[[D", "getFisherInfo", simplify = T)
  res.paramLabels <- .jcall(mtR, "[S", "getParamLabels", simplify = T)
  if(fia){
    res.fia <- .jcall(mtR, "[D", "getFIA", simplify = T)
    res.cfia <- .jcall(mtR, "[D", "getCfia", simplify = T)
  }
  if(jacobian){
    res.jacobian <- .jcall(mtR, "[[[D", "getJacobian", simplify = T)
  }


  colnames(res.paramEst) <- colnames(res.paramEstSE) <- colnames(res.paramEstCI) <- res.paramLabels

  res <- list(paramEstm=res.paramEst, paramSE=res.paramEstSE, paramCI=res.paramEstCI,
              logLik = res.ll, fit=res.fit, df=res.df, p.value=res.p, AIC = res.aic,
              AIC.delta = res.aicDelta, BIC = res.bic, BIC.delta = res.bicDelta,
              fisherInformation = res.fisherinf, FIA = NULL, cFIA = NULL, jacobian = NULL)

  if(fia){
    res$FIA <- res.fia
    res$cFIA <- res.cfia
  }

  if(jacobian){
    colnames(res.jacobian) <- rownames(res.jacobian) <- res.paramLabels
    res$jacobian <- jacobian
  }

  class(res) <- "multiTreeR"
  return(res)
}


#' summary writes summary given mt result object
#' @param mtres the mt result object (\code{\link{doMT}})
summary.multiTreeR <- function(mtres){

  plab <- colnames(mtres$paramEstm)

  # calc CIs
  ci.u <- mtres$paramEstm + mtres$paramCI
  ci.l <- mtres$paramEstm - mtres$paramCI
  colnames(ci.u) <- paste('CI upper (',plab,')', sep='')
  colnames(ci.l) <- paste('CI lower (',plab,')', sep='')
  ci <- cbind(ci.l, ci.u)
  ci <- ci[, c(matrix(1:ncol(ci), nrow = 2, byrow = T))]

  # prepare fitstring
  fit <- cbind(mtres$fit, mtres$df, mtres$p.value, mtres$logLik, mtres$AIC, mtres$BIC, mtres$AIC.delta, mtres$BIC.delta)
  fit <- fit[, c(matrix(1:ncol(fit), nrow = 8, byrow = T))]
  colnames(fit) <- c('PD^lambda', 'df', 'p', 'logLik', 'AIC', 'BIC', 'delta AIC', 'delta BIC')

  # prepare FIA string
  if(!is.null(mtres$FIA)){
    fia <- cbind(mtres$FIA, mtres$cFIA)
    colnames(fia) <- c('FIA', 'cFIA')
  }

  # prep aggr param
  pagg <-t(apply(res$paramEstm, 2, FUN=function(x){c(mean(x),sd(x),quantile(x, c(.025, .50, .975)))}))
  colnames(pagg) <- c('Mean', 'SD', '2.5%', 'Median', '97.5%')

  cat("\nModel Fit:\n\n")
  print(round(fit, 5))

  if(!is.null(mtres$FIA)){
    cat("\n\nMinimum Description Length (Fisher Information Approximation):\n\n")
    print(round(fia, 5))
  }

  cat("\n\nMean parameter estimates:\n\n")
  print(round(pagg, 5))


  cat("\n\nParameter estimates:\n\n")
  print(round(mtres$paramEstm, 5))

  cat("\n\nStandard errors:\n\n")
  print(round(mtres$paramSE, 5))

  cat("\n\nConfidence intervals:\n\n")
  print(round(ci, 5))

}



