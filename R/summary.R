
#' summary writes summary given mt result object
#' @param mtres the mt result object (\code{\link{doMT}})
#' @param digits number digits to round
#' @param ci whether to show confidence intervals
#' @export
#'
summary.multiTreeR <- function(mtres, digits=4, ci=FALSE){

  plab <- colnames(mtres$paramEst)

  if(ci){
    # calc CIs
    ci.u <- mtres$paramEst + mtres$paramCI
    ci.l <- mtres$paramEst - mtres$paramCI
    colnames(ci.u) <- paste('CI upper (',plab,')', sep='')
    colnames(ci.l) <- paste('CI lower (',plab,')', sep='')
    ci <- cbind(ci.l, ci.u)
    ci <- ci[, c(matrix(1:ncol(ci), nrow = 2, byrow = T))]
  }

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
  pagg <-t(apply(res$paramEst, 2, FUN=function(x){c(mean(x),sd(x),quantile(x, c(.025, .50, .975)))}))
  colnames(pagg) <- c('Mean', 'SD', '2.5%', 'Median', '97.5%')

  cat("\nModel Fit:\n\n")
  print(round(fit, digits))

  if(!is.null(mtres$FIA)){
    cat("\n\nMinimum Description Length (Fisher Information Approximation):\n\n")
    print(round(fia, digits))
  }

  cat("\n\nMean parameter estimates:\n\n")
  print(round(pagg, digits))


  cat("\n\nParameter estimates:\n\n")
  print(round(mtres$paramEst, digits))

  cat("\n\nStandard errors:\n\n")
  print(round(mtres$paramSE, digits))

  if(ci){
    cat("\n\nConfidence intervals:\n\n")
    print(round(ci, digits))
  }
}
