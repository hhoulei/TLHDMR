#'
#' @title HDMR
#'
#' @description High-dimensional MR method.
#'
#' @param fdata1 Output from function \emph{TransCondi}.
#' @param cv.lambd A range of parameter \eqn{\lambda}.
#' @param Methodname One of the penalty "Lasso","ElasticNet","MCP" and "SCAD".
#' @param alpha The tuning parameter of "Lasso" and "ElasticNet". If the penalty term is "MCP" or "SCAD", ignore this parameter.
#'
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#'
#' @return A vector, the results of causal effect estimation of each exposure on the outcome.
#'
#' @export
#'
HDMR <- function(fdata1,cv.lambd,penalname,alpha){

  if(penalname %in% c('MCP','SCAD')){

    betaGX_new <- fdata1$betaGX_new
    betaGY_new <- fdata1$betaGY_new

    ################# ncvreg #####################

    ######
    cvfit1 <- ncvreg(betaGX_new,
                     betaGY_new,
                     family = "gaussian",
                     penalty = penalname,
                     lambda=cv.lambd,
                     max.iter = 10000000)
    loc1 <- cvfit1$convex.min
    Estbeta1 <- cvfit1$beta[,loc1]

    # cv.est <- cvfit1$beta
    # cv.est <- as.matrix(cv.est)
    # cv.est <- cv.est[-1,]
    # cv1  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))

    ######
    resall <- Estbeta1[-1]

    P <- ncol(betaGX_new)
    names(resall) <- paste0('X',1:P)


  }else if(penalname %in% c("Lasso","ElasticNet")){

    P <- nrow(fdata1$betaGX_new)

    cvfit <- cv.glmnet(fdata1$betaGX_new, fdata1$betaGY_new,
                       alpha = alpha,
                       lambda=cv.lambd)
    # cv.est <-  cvfit$glmnet.fit$beta
    # cv.est <- as.matrix(cv.est)
    # cv.elas  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))

    fit <- glmnet(fdata1$betaGX_new, fdata1$betaGY_new,
                  alpha = alpha,
                  lambda = cvfit$lambda.min)
    resall <- as.vector(fit$beta)
    #names(Estbeta.elas) <- paste0('X',1:P)

  }

  return(resall)
}


