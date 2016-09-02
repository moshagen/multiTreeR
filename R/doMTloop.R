
#' doMT: estimate mpt model using multiTree
#'
#' Estimate an MPT model using multiTree, but using a loop for data (ensures that FIA is computed correctly for different N per tree)
#'
#' @inheritParams doMT
#' @export
#'
doMTloop  <- function(fEqn, data, catlabels = NULL,
                      restrictions = NULL, fOut = "",
                      lambda = 0.0, numRep=3.0, maxIter = 10000.0, convergence = 1.0E-10, CI=95.0,
                      bootstrapSamples = 0.0, fia=FALSE){

  if(! (is.data.frame(data) | is.matrix(data))){
    stop("data must be a data frame/matrix.")
  }

  ll <- list()
  for(i in 1:nrow(data)){
    ll[[i]] <- doMT(fEqn=fEqn, data[i,,drop=FALSE], catlabels,
                    restrictions, fOut,
                    lambda, numRep, maxIter , convergence , CI,
                    bootstrapSamples , jacobian=FALSE , fia)
  }

  res <- list(paramEst=combn(ll, "paramEst"),
              paramSE=combn(ll, "paramSE"),
              paramCI=combn(ll, "paramCI"),
              logLik = sapply(ll, function(xx) xx$logLik),
              fit=sapply(ll, function(xx) xx$fit),
              df=sapply(ll, function(xx) xx$df),
              p.value=sapply(ll, function(xx) xx$p.value),
              AIC = sapply(ll, function(xx) xx$AIC),
              AIC.delta = sapply(ll, function(xx) xx$AIC.delta),
              BIC = sapply(ll, function(xx) xx$BIC),
              BIC.delta =sapply(ll, function(xx) xx$BIC.delta),
              fisherInformation = NULL,
              FIA = sapply(ll, function(xx) xx$FIA),
              cFIA = sapply(ll, function(xx) xx$cFIA),
              jacobian = NULL)

  ## Note: fisherInformation missing...


  return(res.list)
}


combn <- function(list, what){
  tmp <- lapply(list, function(x) x[[what]])
  mat <- do.call("rbind",tmp)
  mat
}

