#'
#' @title TransCondi
#'
#' @description Transform marginal genetic associations into conditional associations
#'
#' @param fdata A list of G-X and G-Y observed (marginal) associations.
#' Bx_obs is the G-X association coefficient (\eqn{\hat{\beta}_{p}});
#' By_obs is the G-Y association coefficient (\eqn{\hat{\beta}_{Y}});
#' Sigy_condi is the covariance matrix of \eqn{\hat{\beta}_{Y}},  corresponding to the \eqn{\Sigma_{Y}=R\sigma_{Y}\sigma_{Y}^{T}} in the Hou et al.;
#' LD_mat is the LD matrix of IVs (\eqn{R});
#' LD_inv is the inverse of LD matrix (\eqn{R^{-1}}).
#'
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#'
#' @return A list object, including G-X and G-Y conditional genetic associations (\eqn{\hat{\Gamma}} and \eqn{\hat{A}}) and their covariance matrix.
#'
#' @export
#'
TransCondi <- function(fdata){

  ### transform marginal genetic associations into conditional associations ###

  betaGX <- fdata$LD_inv %*% fdata$Bx_obs
  betaGY <- fdata$LD_inv %*% fdata$By_obs
  #sebetaGX <-  LD_inv %*% (LD_mat * crossprod(t(Sigx_condi))) %*% LD_inv
  sebetaGY <- fdata$LD_inv %*% (fdata$LD_mat * crossprod(t(fdata$Sigy_condi))) %*% fdata$LD_inv

  sebetaGY_new <- Matrix12(sebetaGY)
  betaGY_new <- sebetaGY_new %*% betaGY
  betaGX_new <- sebetaGY_new %*% betaGX

  sebetaGX_new <- matrix(1,nrow=nrow(betaGX_new),ncol=ncol(betaGX_new))
  sebetaGY_new <- matrix(1,nrow=nrow(betaGY_new),ncol=ncol(betaGY_new))

  P <- ncol(fdata$Bx_obs)
  colnames(sebetaGX_new) <- paste0('X',1:P)
  colnames(betaGX_new) <- paste0('X',1:P)

  return(list(betaGX_new=betaGX_new,
              betaGY_new=betaGY_new,
              sebetaGX_new=sebetaGX_new,
              sebetaGY_new=sebetaGY_new
  ))

}

Matrix12 <- function(Sigma){

  tol = 1e-10
  eig <- eigen(Sigma, symmetric = TRUE)
  V <- eig$vectors  # 特征向量矩阵
  D <- eig$values   # 特征值向量
  D[D < tol] <- 0
  D_inv_sqrt <- ifelse(D > tol, 1 / sqrt(D), 0)
  result <- V %*% diag(D_inv_sqrt) %*% t(V)

  if (is.complex(result)) result <- Re(result)

  return(result)
}

