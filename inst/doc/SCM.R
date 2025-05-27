## ----fit_mcd, echo=TRUE, eval=FALSE-------------------------------------------
# gam_scm(formula, family = mvn_scm(d = 2, param = NULL, nb = 1), optimizer = NULL, data = list(), aGam = list())

## ----loadingData, echo=TRUE---------------------------------------------------
library(SCM)
data(GEF14_d4)
d <- 4

## ----meanmodelformula, echo=TRUE----------------------------------------------
my_k = 15
my_bs = "cr"
mformula <- list(load_h17 | load_h18 | load_h19 | load_h20  ~ dow + s(doy, k = my_k, bs = my_bs),
                 load_h17 ~ load24_h17 + s(temp95_h17),
                 load_h18 ~ load24_h18 + s(temp95_h18),
                 load_h19 ~ load24_h19 + s(temp95_h19),
                 load_h20 ~ load24_h20 + s(temp95_h20))

## ----covmodelformula, echo=TRUE, eval=TRUE, include=TRUE----------------------
my_k2 = 10
my_bs = "tp"
mformula <- c(mformula, list( 
                 Th_11 | Th_22 | Th_33 | Th_44 | Th_12 | Th_23 | Th_34  ~ dow + s(doy, k = my_k, bs = my_bs),
                 Th_11 ~ s(temp95_h17),
                 Th_22 ~ s(temp95_h18),
                 Th_33 ~ s(temp95_h19),
                 Th_44 ~ s(temp95_h20)))

## ----fit1, echo=TRUE, eval=TRUE, include=TRUE---------------------------------
fit1 <- gam_scm(mformula, family = mvn_scm(d=4), data = GEF14_d4)

## ----gam plot, echo=TRUE, eval=FALSE, include = TRUE--------------------------
# plot(fit1, scale = FALSE, pages = 1)

## ----summary, echo=TRUE, eval=FALSE, include = TRUE---------------------------
# summary(fit1, intercept = FALSE)

## ----residuals, echo=TRUE, eval=FALSE-----------------------------------------
# head(residuals(fit1, type = "deviance"))

## ----pred1, echo=TRUE, eval=FALSE---------------------------------------------
# head(predict(fit1))

## ----pred2, echo=TRUE, eval=FALSE---------------------------------------------
# head(predict(fit1, type = "response"))

## ----pred3, echo=TRUE, eval=TRUE----------------------------------------------
fit1$family$put_cflag(FALSE) 
Sigma_pred <- predict(fit1, type = "response")
head(Sigma_pred)
##        [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
## 1 115.57999 109.29417 107.12426 107.55051 345.00041 380.41531 404.21131
## 2 103.62242  95.79812  92.15967  90.95859 194.88233 211.45315 218.47222
## 3  90.72186  82.53147  78.76657  77.47152  75.38906  82.02693  87.68686
## 4  91.27840  82.81718  78.84653  77.25713  82.63282  88.90992  93.91295
## 5  90.34403  81.23424  76.69174  74.78123  74.13510  80.50152  87.18804
## 6 101.00277  93.29391  88.97594  87.38839 146.50549 154.48771 162.74744
##        [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
## 1 430.83223 361.22653 370.78073 391.54811 380.85143 402.90381 416.52786
## 2 224.53788 202.14192 204.03489 214.30886 205.08170 216.00432 220.79684
## 3  93.52467  78.10238  79.82161  84.38060  81.34185  86.35829  90.14268
## 4  98.78637  85.00924  86.13507  90.80895  86.98921  92.20045  95.88335
## 5  94.53048  76.63749  78.63748  83.23227  80.58371  85.71751  90.29103
## 6 172.51746 149.85170 152.78919 158.12093 155.89525 161.74611 166.90203

## ----pred4, echo=TRUE, eval=TRUE----------------------------------------------
Sigma_mat(Sigma_pred[,-c(1:d)])[[1]]
##          [,1]     [,2]     [,3]     [,4]
## [1,] 345.0004 361.2265 370.7807 380.8514
## [2,] 361.2265 380.4153 391.5481 402.9038
## [3,] 370.7807 391.5481 404.2113 416.5279
## [4,] 380.8514 402.9038 416.5279 430.8322

