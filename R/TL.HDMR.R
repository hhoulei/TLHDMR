#'
#' @title TL.HDMR
#'
#' @description Transfer Learning for High-Dimensional Mendelian Randomization (Lei Hou)
#'
#' @param betaGX_A0 A covariance-corrected matrix of G-X conditional association coefficients (nSNP*nexposure) from the target population, corresponding to the \eqn{\Omega_{Y}^{-1/2}\hat{A}^{(0)}} in the Hou et al.
#'
#' @param betaGY_A0 A covariance-corrected vector of G-Y conditional association coefficients (nSNP*1) from the target population, corresponding to the \eqn{\Omega_{Y}^{-1/2}\hat{\Gamma}^{(0)}} in the Hou et al.
#'
#' @param betaGX_AK A covariance-corrected matrix of G-X conditional association coefficients (nSNP*nexposure) from the source population, corresponding to the \eqn{\Omega_{Y}^{-1/2}\hat{A}^{(k)}} in the Hou et al.
#'
#' @param betaGY_AK  A covariance-corrected vector of G-Y conditional association coefficients (nSNP*1) from the source population, corresponding to the \eqn{\Omega_{Y}^{-1/2}\hat{\Gamma}^{(k)}} in the Hou et al.
#'
#' @param cv.lambd A range of parameter \eqn{\lambda}.
#'
#' @param Methodname One of the penalty "Lasso","ElasticNet","MCP" and "SCAD".
#'
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#'
#' @return A list object, the results of causal effect estimation of each exposure on the outcome as well as the value lambda.
#'
#'
#' @export
#'
TL.HDMR <- function(betaGX_A0,betaGY_A0,
                    betaGX_AK,betaGY_AK,
                    cv.lambd,Methodname){

  p <- ncol(betaGX_AK)
  g <- nrow(betaGX_AK)
  g0 <- nrow(betaGX_A0)

  if(Methodname=='Lasso' | Methodname=='ElasticNet'){

    if(Methodname=='Lasso') alpha=1
    if(Methodname=='ElasticNet') alpha=0.5

    ### all data containing the source data
    cv.init <- cv.glmnet(betaGX_AK,
                         betaGY_AK,
                         alpha = alpha,
                         #lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/g)
                         lambda=cv.lambd*sqrt(2*log(p)/g))
    # cv.lam.const <- cv.init$lambda/sqrt(2*log(p)/g)
    # cv.w.kA <- as.matrix(cv.init$glmnet.fit$beta)
    # cv.w.kA<- t(apply(cv.w.kA,1,function(x) x*(abs(x)>=cv.init$lambda)))


    lam.const <- cv.init$lambda.min/sqrt(2*log(p)/g)
    loc.min <- which(cv.init$lambda==cv.init$lambda.min)
    w.kA <- as.numeric(glmnet(betaGX_AK,
                              betaGY_AK,
                              alpha = alpha,
                              lambda=lam.const*sqrt(2*log(p)/g))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/g))

    ### only source data
    # cv.delta.kA <- glmnet(x=betaGX_A0,
    #                       y=betaGY_A0-betaGX_A0%*%w.kA,
    #                       alpha = alpha,
    #                       lambda=cv.lam.const*sqrt(2*log(p)/g0))
    # cv.delta.kA <- as.matrix(cv.delta.kA$beta)
    # cv.delta.kA<- t(apply(cv.delta.kA,1,function(x) x*(abs(x)>=cv.lam.const*sqrt(2*log(p)/g0))))

    delta.kA <- as.numeric(glmnet(x=betaGX_A0,
                                  y=betaGY_A0-betaGX_A0%*%w.kA,
                                  alpha = alpha,
                                  lambda=lam.const*sqrt(2*log(p)/g0))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/g0))

    # combine
    beta.kA <- w.kA + delta.kA
    # cv.beta.kA <- cv.w.kA + cv.delta.kA

  }else
    if(Methodname=='MCP' | Methodname=='SCAD'){

      penalname <- Methodname

      ### all data containing the source data
      cv.init <- ncvreg(betaGX_AK,
                        betaGY_AK,
                        family = "gaussian",
                        penalty = penalname,
                        lambda=cv.lambd*sqrt(2*log(p)/g),
                        max.iter = 10000000)

      cv.lam.const <- cv.init$lambda/sqrt(2*log(p)/g)
      cv.w.kA <- as.matrix(cv.init$beta)
      cv.w.kA <- cv.w.kA[-1,]
      cv.w.kA<- t(apply(cv.w.kA,1,function(x) x*(abs(x)>=cv.init$lambda)))

      loc.min <- cv.init$convex.min
      lam.const <- cv.init$lambda[loc.min]/sqrt(2*log(p)/g)
      w.kA <- cv.w.kA[,loc.min]

      ### only source data
      cv.delta.kA <- ncvreg(betaGX_A0,
                            betaGY_A0-betaGX_A0%*%w.kA,
                            family = "gaussian",
                            penalty = penalname,
                            lambda=cv.lam.const*sqrt(2*log(p)/g0),
                            max.iter = 10000000)
      cv.delta.kA <- as.matrix(cv.delta.kA$beta)
      cv.delta.kA <- cv.delta.kA[-1,]
      cv.delta.kA<- t(apply(cv.delta.kA,1,function(x) x*(abs(x)>=cv.lam.const*sqrt(2*log(p)/g0))))

      delta.kA <- cv.delta.kA[,loc.min]

      # combine
      beta.kA <- w.kA + delta.kA
      #cv.beta.kA <- cv.w.kA + cv.delta.kA

    }

  #cv.beta.kA1  <- apply(cv.beta.kA,2,function(x) Metrics(as.numeric(x!=0),trueX))

  return(list(betall=beta.kA,
              lam.const=lam.const))

}
