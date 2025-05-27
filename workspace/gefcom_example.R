

jacobian2 <- function(eta, jj){

  res <- matrix(0, nrow(eta), no_eta) # Jacobian matrix
  Cm <- internal()$mat2vec(d) #Matrix C (see the UK paper)

  cor_flag <-  as.integer(getcflag())

  if ( jj <= d ) {
    res[, jj] <- 1 # Mean vector case
  } else { # Covariance/Correlation matrix case
    #jj <- jj - d
    idx_jj <- which(Cm == (jj - d), arr.ind = TRUE) # identify the row and the column of the C matrix associated the linear predictor index
    S_r <- as.numeric(idx_jj[1, 1]) - 1 # extract the row
    S_c <- as.numeric(idx_jj[1, 2]) - 1 # extract the column

    if ( param == 1 ) { #mcd
      rc_idx_s <- rc_idx_t <- rep(NA, d * (d - 1)/2)
      count <- 1
      for ( j in (d + 1) : (d * (d + 1)/2) ) {  # identify the rows and the columns of the C matrix associated the linear predictor index
        rc_idx_s[count] <- as.numeric(which(Cm == j, arr.ind = TRUE)[1, 1]) - 1
        rc_idx_t[count] <- as.numeric(which(Cm == j, arr.ind = TRUE)[1, 2]) - 1
        count <- count + 1
      }
      internal()$jacobian_mcd(eta, res, d, S_r, S_c, rc_idx_s, rc_idx_t, cor_flag)
    }

    if ( param == 2 ) internal()$jacobian_logm(eta, res, d, S_r, S_c, cor_flag) #logm
  }

  if(param == 1){
    idx_jac <- function(j, d){
      if (j <= (d + 1) ) idx <- j
      if (j > (d + 1) &  j <= 2*d) {
        lidx <- j - d + choose(j - d, 2)
        idx <- rep(0, lidx)
        idx[1 : (j-d)] <- (d + 1): j
        idx[(j - d + 1) : lidx] <- t(Gm)[1 : (j - d -1), 1 : (j - d -1)][upper.tri(t(Gm)[1:(j - d -1),1:(j - d -1)], diag = TRUE)] + 1

      }
      if (j > 2*d) {
        lidx <- w[j - 2*d] + 1 + choose(w[j - 2*d] + 1, 2)
        idx <- rep(0, lidx)
        idx[1 : (w[j - 2*d] + 1)] <- d + (1 : (w[j - 2*d]+1))
        idx[(w[j - 2*d] + 2): lidx] <- t(Gm)[1 : w[j - 2*d], 1 : w[j - 2*d]][upper.tri(t(Gm)[1:w[j - 2*d],1:w[j - 2*d]], diag = TRUE)] + 1
      }
      return( idx )
    }

    eta_idx <- idx_jac(jj, d)
    eta_deriv <- as.matrix(res[,eta_idx])
  }

  if(param == 2){
    if(jj <= d){
      eta_deriv <- as.matrix(res[,jj])
    } else {
      eta_idx <- (d+1) : no_eta
      eta_deriv <- res[,eta_idx ]
    }

  }

  return(list(DmuDeta = eta_deriv, eta_idx = eta_idx))
} ## end jacobian2


load("GEF14_mvn_d6_17_22.RData")

d <- 4
my_k = 15
my_bs = "tp"
form_full <- list(load_h17 | load_h18  ~ dow + s(doy, k = my_k, bs = my_bs),
                  load_h19 | load_h20  ~ dow + s(doy, k = my_k, bs = my_bs),
                  Th_11 | Th_22 | Th_33 | Th_44  ~ dow + s(doy, k = my_k, bs = my_bs),
                  Th_12 | Th_23 | Th_34 ~ s(doy, k = my_k, bs = my_bs),
                  load_h17 ~ load24_h17 + s(temp95_h17),
                  load_h18 ~ load24_h18 + s(temp95_h18),
                  load_h19 ~ load24_h19 + s(temp95_h19),
                  load_h20 ~ load24_h20 + s(temp95_h20),
                  Th_11 ~ s(temp95_h17),
                  Th_22 ~ s(temp95_h18),
                  Th_33 ~ s(temp95_h19),
                  Th_44 ~ s(temp95_h20),
                  Th_24 ~ dow,
                  Th_13 ~ s(doy))

fit <- gam_scm(form_full, family = mvn_scm(d, nb = 1, param = "mcd"), data = GEF14_data_d6)

ETA <- matrix(rnorm(14*2),2,14, byrow=TRUE)

library(mgcViz)
fit <- getViz(fit, nsim = 100)

plot(ALE(fit, "temp95_h17", oind = 10, type = "response")) + l_fitLine() + l_ciLine()

# hist(residuals(fit)[ , 2])
# hist(residuals(fit, type = "response")[ , 2])
#
# cov(residuals(fit, type = "response"))
# cov(fit$family$rd(predict(fit)) - fit$fitted.values[ , 1:4])

sim <- simulate(fit, nsim = 10)

res <- lapply(sim, function(x) x - fit$fitted.values[ , 1:4])

hist(sapply(res, function(x) x[ , 2]))
hist((fit$y - fit$fitted.values[ , 1:4])[ , 2])

str(simulate(fit, nsim = 2))

# Histogram of observed and simulated residuals for 1st response variable
check0D(fit, type = "response", trans = function(x) x[ , 1]) + l_hist() + l_dens1D() + l_rug()

check0D(fit, type = "response", trans = function(x) cor(x[ , 1], x[ , 2])) + l_hist()  + l_vline() + l_dens1D() + l_rug()

# Binned variance of observed and simulated residuals of 2nd response variable
check1D(fit, "doy", type = "response", trans = function(x) x[ , 2]) + l_gridCheck1D(gridFun = sd)

# Binned covariance between observed and simulated residuals of 2nd and 3rd response variable
check1D(fit, "doy", type = "response", trans = function(x) x[ , 2]*x[ , 3]) + l_gridCheck1D()

# Binned correlation between observed and simulated residuals of 2nd and 3rd response variable
check1D(fit, "doy", type = "pearson", trans = function(x) x[ , 2]*x[ , 3]) + l_gridCheck1D()

check1D(fit, "doy", type = "response", trans = function(x) x[ , 2]) + l_densCheck()

# MEASURE OF TAIL DEPENDENCE???
check0D(fit, type = "response", trans = function(x) cor(x[ , 1], x[ , 2]))


pl <- check1D(fit, "doy", type = "response")
pl + l_gridCheck1D(gridFun = function(.x) skewness(.x[ , 1], .x[ , 4]))

pl + l_gridQCheck1D(qu = 0.5)

pl + l_densCheck()

library(microbenchmark)
microbenchmark(ALE(fit, "temp95_h17", oind = 10, type = "response"),
               ALE(fit2, "temp95_h17", oind = 10, type = "response"),
               ALE(fit, "temp95_h17", oind = 10, type = "link"),
               ALE(fit2, "temp95_h17", oind = 10, type = "link"), times = 20)

fit$family$jacobian <- jacobian2
environment(fit$family$jacobian2) <- environment(fit$family$jacobian)
fit$family$jacobian2(ETA, jjjj)

environment(fit$family$rd) <- environment(fit$family$ll)
fit$family$rd <- function (mu, wt, scale)
{
  no_eta <- ncol(mu)
  d <- -3/2 + sqrt(9/4 + 2 * no_eta)
  out <- matrix(NA, nrow(mu), d)
  if (param == 1) {
    for (i in 1:nrow(mu)) {
      Dm2 <- diag(exp(-mu[i,(d+1):(2*d)]))
      T <- matrix(0, d, d)
      count <- 2*d + 1
      for(j in 2:d){
        for(k in 1:(j-1)){
          T[j,k] <- mu[i, count]
          count <- count + 1
        }
      }
      diag(T) <- rep(1,d)
      Sigma <- solve(t(T) %*% Dm2 %*% T)
      C <- chol(Sigma)
      u <- rmvn(1, rep(0, d), diag(rep(1, d)))
      out[i,] <- t(mu[i, 1 : d] + t(C) %*% t(u))

    }
  }
  if (param == 2) {
    for (i in 1:nrow(mu)) {
      Sigma <- internal()$logM_Sigma(mu[i, ], d)
      C <- t(chol(Sigma))
      u <- rmvn(1, rep(0, d), diag(rep(1, d)))
      out[i, ] <- t(mu[i, 1:d] + C %*% t(u))
    }
  }
  return(out)
}

fit$family$residuals <- function(object, type = c("response", "deviance")) { #by defualt deviance residuals
  type <- match.arg(type)

  if(type == "pearson"){
    n <- dim(object$fitted.values)[1]
    res <- object$y - object$fitted.values[, 1 : d]
  }

  if ( type == "deviance" ) {
    n <- dim(object$fitted.values)[1]
    res <- matrix(0, n, d)
    if ( param == 1 ) SCM:::internal()$res_dev_mcd(object$fitted.values, object$y, res) #Deviance residuals for mcd
    if ( param == 2 ) SCM:::internal()$res_dev_logm(object$fitted.values, object$y, res) #Deviance residuals for logm
  } else {
    res <- object$y - object$fitted.values[, 1 : d]
  }
  res
} ## residuals


