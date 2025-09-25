#'
#' @title HDMR.TSD
#'
#' @description Transferable Source Detection for TL.HDMR.
#'
#' @param fdata1 Conditional GWAS summary statistics for the target dataset. Output from function \emph{TransCondi}.
#' @param fdata_NK1 A list of conditional GWAS summary statistics for source datasets. Output from function \emph{TransCondi}.
#' @param cv.lambd A range of parameter \eqn{\lambda}.
#' @param Methodname One of the penalty "Lasso","ElasticNet","MCP" and "SCAD".
#' @param r.A0 The number of folds for target dataset in the TSD algorithm.
#' @param C0 The threshold in the TSD algorithm.
#'
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#'
#' @return A list object, \emph{Ind.loc} is the index of the informative source datasets.
#'
#' @export
#'
HDMR.TSD <- function(fdata1,fdata_NK1,cv.lambd,Methodname,r.A0,C0){

  #r.A0 <- 3
  g0 <- nrow(fdata1$betaGX_new)
  ind.A0 <- 1:g0
  n.AK <- length(fdata_NK1)

  loss.A0 <- NULL
  loss.Ak <- matrix(NA,nrow=r.A0,ncol=n.AK)
  for(r in 1:r.A0){

    loc <- sample(ind.A0,floor(g0/r.A0))
    ind.A0 <- ind.A0[!(ind.A0%in%loc)]
    betaGX_A0 <- fdata1$betaGX_new[-loc,]
    betaGY_A0 <- fdata1$betaGY_new[-loc]

    betaGX_Ar <- fdata1$betaGX_new[loc,]
    betaGY_Ar <- fdata1$betaGY_new[loc]

    #### calculate beta.A0

    beta.A0 <- Method.A0(betaGX_A0,betaGY_A0,cv.lambd,Methodname)
    Q.A0 <- sum((betaGY_Ar-betaGX_Ar%*%beta.A0)^2)/(2*g0)
    loss.A0 <- c(loss.A0,Q.A0)

    #### calculate beta.Ak

    for(nk in 1:n.AK){

      betaGX_Ak <- fdata_NK1[[nk]]$betaGX_new
      betaGY_Ak <- fdata_NK1[[nk]]$betaGY_new

      betaGX_AK <- rbind(betaGX_A0,betaGX_Ak)
      betaGY_AK <- c(betaGY_A0,betaGY_Ak)
      beta.Ak.nk <- TL.HDMR(betaGX_A0,
                            betaGY_A0,
                            betaGX_AK,
                            betaGY_AK,
                            cv.lambd,
                            Methodname)
      beta.Ak.nk <- beta.Ak.nk$betall
      Q.Ak <- sum((betaGY_Ar-betaGX_Ar%*%beta.Ak.nk)^2)/(2*g0)

      loss.Ak[r,nk] <- Q.Ak
    }


  }

  loss.A0.mean <- mean(loss.A0)
  loss.A0.sd <- sqrt(sum(0.5*(loss.A0-loss.A0.mean)^2))
  loss.Ak.mean <- apply(loss.Ak,2,mean)

  diff.loss <- abs(loss.Ak.mean-loss.A0.mean)
  Ind.loc <- which(diff.loss<=C0*max(loss.A0.sd,0.01))

  return(list(loss.A0=loss.A0,
              loss.Ak=loss.Ak,
              loss.A0.sd=loss.A0.sd,
              diff.loss=diff.loss,
              Ind.loc=Ind.loc))
}


Method.A0 <- function(betaGX_A0,betaGY_A0,cv.lambd,Methodname){

  if(Methodname=='Lasso' | Methodname=='ElasticNet'){

    if(Methodname=='Lasso') alpha=1
    if(Methodname=='ElasticNet') alpha=0.5

    cv.init <- cv.glmnet(betaGX_A0,
                         betaGY_A0,
                         alpha = alpha,
                         lambda=cv.lambd)
    w.kA <- glmnet(betaGX_A0,
                   betaGY_A0,
                   alpha = alpha,
                   lambda=cv.init$lambda.min)
    beta <- as.numeric(w.kA$beta)

  }else
    if(Methodname=='MCP' | Methodname=='SCAD'){

      penalname <- Methodname

      cv.init <- ncvreg(betaGX_A0,
                        betaGY_A0,
                        family = "gaussian",
                        penalty = penalname,
                        lambda=cv.lambd,
                        max.iter = 10000000)
      cv.w.kA <- as.matrix(cv.init$beta)
      cv.w.kA <- cv.w.kA[-1,]
      beta <- cv.w.kA[,cv.init$convex.min]

    }
  return(beta)
}




