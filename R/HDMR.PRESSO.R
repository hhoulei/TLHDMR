#'
#' @title HDMR.PRESSO
#'
#' @description Removing invalid IVs for TL.HDMR.
#'
#' @param fdata1 Conditional GWAS summary statistics for the target dataset. Output from function \emph{TransCondi}.
#' @param Methodname One of the penalty "Lasso","ElasticNet","MCP" and "SCAD".
#' @param cv.lambd A range of parameter \eqn{\lambda}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#'
#' @return A list object, \emph{removeIV} is the index of invalid IVs. \emph{fdata2} is the new datasets removing invalid IVs.
#'
#' @export
#'
HDMR.PRESSO <- function(fdata1,Methodname,cv.lambd){

  methth <- Methodname

  bx_corMatrxi <- cor(fdata1$betaGX_new)

  giv <- nrow(fdata1$betaGX_new)
  Piv <- ncol(fdata1$betaGX_new)
  Pvalue <- NULL

  Pvalue <- NULL
  for(gi in 1:giv){

    #cat('gi=',gi,'\n')

    fdata_once <- fdata1

    fdata_once$betaGX_new <- fdata_once$betaGX_new[-gi,]
    fdata_once$betaGY_new <- fdata_once$betaGY_new[-gi]
    fdata_once$sebetaGX_new <- fdata_once$sebetaGX_new[-gi,]
    fdata_once$sebetaGY_new <- fdata_once$sebetaGY_new[-gi]

    if(methth=='MCP'){
      res_cisMVMR.MCP <- HDMR(fdata1,cv.lambd,penalname='MCP',alpha=0)
      thtaj <- res_cisMVMR.MCP
    }else if(methth=='SCAD'){
      res_cisMVMR.SCAD <- HDMR(fdata1,cv.lambd,penalname='SCAD',alpha=0)
      thtaj <- res_cisMVMR.SCAD
    }else if(methth=='Lasso'){
      res.LDAlasso <- HDMR(fdata1,cv.lambd,penalname='Lasso',alpha=1)
      thtaj <- res.LDAlasso
    }else if(methth=='ElasticNet'){
      res.elas <- HDMR(fdata1,cv.lambd,penalname='ElasticNet',alpha=0.5)
      thtaj <- res.elas
    }

    predgammaj <- sum(fdata1$betaGX_new[gi]*thtaj)
    RSSobsj <- fdata1$betaGY_new[gi]-predgammaj

    Aj <- fdata1$betaGX_new[gi,]
    seAj <- fdata1$sebetaGX_new[gi,]
    bx_cor_Matrxionce <- bx_corMatrxi
    for(ro in 1:Piv){
      for(co in 1:Piv){
        bx_cor_Matrxionce[ro,co] <- seAj[ro]*seAj[co]*bx_corMatrxi[ro,co]
      }
    }

    nboot <- 100
    Pvalj <- 0
    for(kp in 1:nboot){
      Ajrandom <- mvrnorm(1,Aj,bx_cor_Matrxionce)
      gammarandom <- rnorm(1,predgammaj,fdata1$sebetaGY_new[gi])
      RSSexpjk <- gammarandom-sum(Ajrandom*thtaj)
      Pvalj <- Pvalj + as.numeric(RSSexpjk>RSSobsj)
    }
    Pvalj <- Pvalj/nboot

    Pvalue <- c(Pvalue,Pvalj)

  }

  Pcut <- (0.05/giv)

  fdata2 <- fdata1
  fdata2$betaGX <- fdata2$betaGX[Pvalue>=Pcut,]
  fdata2$betaGY <- fdata2$betaGY[Pvalue>=Pcut]
  fdata2$betaGX_new <- fdata2$betaGX_new[Pvalue>=Pcut,]
  fdata2$betaGY_new <- fdata2$betaGY_new[Pvalue>=Pcut]
  fdata2$sebetaGX_new <- fdata2$sebetaGX_new[Pvalue>=Pcut,]
  fdata2$sebetaGY_new <- fdata2$sebetaGY_new[Pvalue>=Pcut]

  return(list(fdata2=fdata2,
              Pvalue=Pvalue,
              removeIV=which(Pvalue<Pcut)))

}


