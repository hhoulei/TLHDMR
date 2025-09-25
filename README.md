# TLHDMR

Transfer Learning for Mendelian Randomization with Pleiotropic and Correlated High-Dimensional Exposures and Its Application to Trans-ethnic Populations  

***Installation***  
`devtools::install_github("hhoulei/TLHDMR")`  

***Toy Example***  

`library(MASS)`  
`library(glmnet)`    
`library(ncvreg)`  
`library(TL.HDMR)`  

`data("fdata")`  
`data("fdata_NK")`  

`fdata1 <- TransCondi(fdata)`  

`fdata_NK1 <- list()`  
`for(oj in 1:length(fdata_NK)){`  
`  once <- fdata_NK[[oj]]`  
`  fdata_NK1[[oj]] <- TransCondi(once)`  
`}`  

`cv.lambd <- exp(c(seq(-10,-2,1)))`  

`##################### remove invalid IV #############`  
`fdata1j <- HDMR.PRESSO(fdata1,`  
`                       Methodname='MCP',`  
`                       cv.lambd)`  

`fdata1j <- fdata1j$fdata2`  
`fdata_NK1o <- list()`  
`for(kpi in 1:length(fdata_NK1)){`  
`  cat('kpi=',kpi,'\n')`  
`  onon <- HDMR.PRESSO(fdata_NK1[[kpi]],`  
`                      Methodname='MCP',`  
`                      cv.lambd)`  
`  fdata_NK1o[[kpi]] <- onon$fdata2`  
`}`  

`##################### select informative source datasets #############`  

`res.select <- HDMR.TSD(fdata1,`  
`                       fdata_NK1o,`  
`                       cv.lambd,`  
`                       Methodname='MCP',`  
`                       r.A0=3,C0=1.6)`  

`sel.NK<- res.select$Ind.loc`  
`fdata_NK2o <- fdata_NK1o[sel.NK]`  

`##################### transfer learning #############`  
`cat('Transfer learning...','\n')`  
`#> Transfer learning...`  
`DATA_condi <- function(fdata1,fdata_NK0){`  
  
`  betaGX_A0 <- fdata1$betaGX_new`  
`  betaGY_A0 <- fdata1$betaGY_new`  
  
`  betaGX_AK <- betaGX_A0`  
`  betaGY_AK <- betaGY_A0`  
`  for(oak in 1:length(fdata_NK0)){`  
`    betaGX_AK <- rbind(betaGX_AK,fdata_NK0[[oak]]$betaGX_new)`  
`    betaGY_AK <- c(betaGY_AK,fdata_NK0[[oak]]$betaGY_new)`  
`  }`  
  
`  return(list(betaGX_A0=betaGX_A0,`  
`              betaGY_A0=betaGY_A0,`  
`              betaGX_AK=betaGX_AK,`  
`              betaGY_AK=betaGY_AK))`  
`}`  
`fdata_condi <- DATA_condi(fdata1j,fdata_NK2o)`  

`cisMVMR.Trans.MCP.presso <- TL.HDMR(fdata_condi$betaGX_A0,fdata_condi$betaGY_A0,`  
`                                               fdata_condi$betaGX_AK,fdata_condi$betaGY_AK,`  
`                                               cv.lambd,Methodname='MCP')`  
`cisMVMR.Trans.MCP.presso`  

***Citation***:  
Transfer Learning for Mendelian Randomization with Pleiotropic and Correlated High-Dimensional Exposures and Its Application to Trans-ethnic Populations  
Lei Hou, Xiaohua Zhou, Fuzhong Xue, Hao Chen

Please contact houlei@bicmr.pku.edu.cn for any questions. We will continue to update this R package and reduce the problems that may be encountered during its installation.

