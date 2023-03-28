#Required Packages
#install_github("cran/EcoHydRology")
#library("EcoHydRology")
#install.packages(c("tcltk", "lubridate", "dplyr", "tidyr", "pastecs", "extRemes", "svDialogs", "fExtremes", "evd", "copula", "RcppRoll"))
#install.packages(c("zoo", "ismev", "VGAM", "SciViews", "lmomco", "copBasic", "hexbin", "scatterplot3d", "rgl", "graphics", "gumbel", "SuppDists", "openxlsx", "gnFit"))
#install.packages(c("ggforce", "ggpubr", "ggExtra", "ggridges", "ggplot2", "tidyverse", "hrbrthemes", "viridis", "kit", "eva", "VineCopula", "POT", "stringr", "plotrix", "mev", "logspline"))
lapply(c("tcltk", "lubridate", "dplyr", "tidyr", "pastecs", "extRemes", "svDialogs", "fExtremes", "evd", "copula", "RcppRoll"), require, character.only = TRUE)
lapply(c("ggforce", "ismev", "VGAM", "SciViews", "lmomco", "copBasic", "hexbin", "scatterplot3d", "rgl", "graphics", "gumbel", "SuppDists", "openxlsx"), require, character.only = TRUE)
lapply(c("ggpubr", "ggExtra", "ggridges", "ggplot2", "tidyverse", "hrbrthemes", "viridis", "kit", "eva", "zoo", "VineCopula", "POT", "stringr", "plotrix", "mev", "logspline", "gnFit"), require, character.only = TRUE)

#To clean
rm(list = ls())
graphics.off()

#Loading data (Summary_Final.csv)
#data = read.csv(file.choose())
data = read.csv(tk_choose.files(default = "", caption = "Select csv file"))

#Define functions
firstup = function(x){
  
  substr(x, 1, 1) = toupper(substr(x, 1, 1))
  x
  
}
gevrProfLoc = function(z, conf = 0.95, plot = TRUE, limi_1 = 0.93, limi_2 = 1.08, opt = c("Nelder-Mead"), nub = 100){
  
  if (z$gumbel | !z$stationary) 
    stop("Object cannot be from a Gumbel and/or a nonstationary fit!")
  data <- as.matrix(z$data)
  theta <- as.numeric(z$par.ests)
  opt <- match.arg(opt)
  sol <- c(theta[2], theta[3])
  gevrLikShape <- function(a, sh) {
    if (a[2] <= 0) {
      #out <- .Machine$double.xmax
      
      out <- dgevr(data, loc = sh, scale = a[1], shape = a[2], 
                   log.d = TRUE)
      out <- -sum(out)
      if (out == Inf) 
        out <- .Machine$double.xmax
    }
    else {
      out <- dgevr(data, loc = sh, scale = a[1], shape = a[2], 
                   log.d = TRUE)
      out <- -sum(out)
      if (out == Inf) 
        out <- .Machine$double.xmax
    }
    out
  }
  cutoff <- qchisq(conf, 1)
  prof <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    2 * (lmax - lci) - cutoff
  }
  proff <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    lci
  }
  prof <- Vectorize(prof)
  proff <- Vectorize(proff)
  suppressWarnings(out1 <- uniroot(prof, c(theta[1] - 1e-06, theta[1]), extendInt = "downX"))
  suppressWarnings(out2 <- uniroot(prof, c(theta[1], theta[1] + 1e-06), extendInt = "upX"))
  CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  if(plot){
    
    prof1 <- function(sh){proff(sh)}
    curve(prof1, from = CI[1]*limi_1, to = CI[2]*limi_2, n = nub, xlab = "Location", ylab = "")
    title(main = "Profile Likelihood Function")
    abline(v = theta[1], col = "blue")
    abline(v = CI[1], col = "blue")
    abline(v = CI[2], col = "blue")
    abline(h = sum(dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE))-cutoff/2, lty = 2)
    
  }
  out <- list(theta[1], CI, conf)
  names(out) <- c("Estimate", "CI", "ConfLevel")
  out
}
gevrProfScale = function(z, conf = 0.95, plot = TRUE, limi_1 = 0.93, limi_2 = 1.08, opt = c("Nelder-Mead"), nub = 100){
  
  if (z$gumbel | !z$stationary) 
    stop("Object cannot be from a Gumbel and/or a nonstationary fit!")
  data <- as.matrix(z$data)
  theta <- as.numeric(z$par.ests)
  opt <- match.arg(opt)
  sol <- c(theta[3], theta[1])
  gevrLikShape <- function(a, sh) {
    if (a[2] <= 0) {
      #out <- .Machine$double.xmax
    }else {
      out <- dgevr(data, loc = a[2], scale = sh, shape = a[1], 
                   log.d = TRUE)
      out <- -sum(out)
      if (out == Inf) 
        out <- .Machine$double.xmax
    }
    out
  }
  cutoff <- qchisq(conf, 1)
  prof <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    2 * (lmax - lci) - cutoff
  }
  proff <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    lci
  }
  prof <- Vectorize(prof)
  proff <- Vectorize(proff)
  suppressWarnings(out1 <- uniroot(prof, c(theta[2] - 1e-06, theta[2]), extendInt = "downX"))
  suppressWarnings(out2 <- uniroot(prof, c(theta[2], theta[2] + 1e-06), extendInt = "upX"))
  CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  if(plot) {
    prof1 <- function(sh){proff(sh)}
    curve(prof1, from = CI[1]*limi_1, to = CI[2]*limi_2, n = nub, xlab = "Scale", ylab = "")
    title(main = "Profile Likelihood Function")
    abline(v = theta[2], col = "blue")
    abline(v = CI[1], col = "blue")
    abline(v = CI[2], col = "blue")
    abline(h = sum(dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE))-cutoff/2, lty = 2)
  }
  out <- list(theta[2], CI, conf)
  names(out) <- c("Estimate", "CI", "ConfLevel")
  out
}
gevrProfShape = function(z, conf = 0.95, plot = TRUE, limi_1 = 0.93, limi_2 = 1.08, opt = c("Nelder-Mead"), nub = 100){
  
  if (z$gumbel | !z$stationary) 
    stop("Object cannot be from a Gumbel and/or a nonstationary fit!")
  data <- as.matrix(z$data)
  theta <- as.numeric(z$par.ests)
  opt <- match.arg(opt)
  sol <- c(theta[1], theta[2])
  gevrLikShape <- function(a, sh) {
    if (a[2] <= 0) {
      out <- .Machine$double.xmax 
    }
    else {
      out <- dgevr(data, loc = a[1], scale = a[2], shape = sh, 
                   log.d = TRUE)
      out <- -sum(out)
      if (out == Inf) 
        out <- .Machine$double.xmax
    }
    out
  }
  cutoff <- qchisq(conf, 1)
  prof <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    2 * (lmax - lci) - cutoff
  }
  proff <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    lci
  }
  prof <- Vectorize(prof)
  proff <- Vectorize(proff)
  suppressWarnings(out1 <- uniroot(prof, c(theta[3] - 1e-06, theta[3]), extendInt = "downX"))
  suppressWarnings(out2 <- uniroot(prof, c(theta[3], theta[3] +  1e-06), extendInt = "upX"))
  CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  if (plot) {
    prof1 <- function(sh){proff(sh)}
    curve(prof1, from = CI[1]*limi_1, to = CI[2]*limi_2, n = nub, xlab = "Shape", ylab = "")
    title(main = "Profile Likelihood Function")
    abline(v = theta[3], col = "blue")
    abline(v = CI[1], col = "blue")
    abline(v = CI[2], col = "blue")
    abline(h = sum(dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE))-cutoff/2, lty = 2)
  }
  out <- list(theta[3], CI, conf)
  names(out) <- c("Estimate", "CI", "ConfLevel")
  out
}
gevrRl = function(z, period, conf = 0.95, nub = 100, method = c("delta", "profile"), plot = TRUE, limi_1 = 0.93, limi_2 = 1.08, opt = c("Nelder-Mead")){
  
  if (!z$stationary) 
    stop("Return levels can only be produced for the stationary model!")
  method <- match.arg(method)
  data <- as.matrix(z$data)
  theta <- z$par.ests
  cov <- z$varcov
  m <- -log1p(-(1/period))
  if (!z$gumbel) 
    est <- theta[1] - (theta[2]/theta[3]) * (-expm1(-theta[3] * 
                                                      log(m)))
  else est <- theta[1] - theta[2] * log(m)
  est <- as.numeric(est)
  if (method == "delta") {
    if (!z$gumbel) {
      del <- matrix(ncol = 1, nrow = 3)
      del[1, 1] <- 1
      del[2, 1] <- -((theta[3])^(-1)) * (-expm1(-theta[3] * 
                                                  log(m)))
      del[3, 1] <- ((theta[2]) * (theta[3]^(-2)) * (-expm1(-theta[3] * 
                                                             log(m)))) - ((theta[2]) * ((theta[3])^(-1)) * 
                                                                            ((1 + expm1(-theta[3] * log(m))) * log(m)))
      se <- as.numeric(sqrt(t(del) %*% cov %*% del))
    }
    else {
      se <- as.numeric(sqrt(cov[1, 1] - ((cov[2, 2] + cov[2, 
                                                          1]) * log(m)) + (cov[2, 2] * (log(m))^2)))
    }
    alpha <- (1 - conf)/2
    lower <- est - qnorm(1 - alpha) * se
    upper <- est + qnorm(1 - alpha) * se
    CI <- c(lower, upper)
  }
  else {
    opt <- match.arg(opt)
    if (!z$gumbel) 
      sol <- c(theta[2], theta[3])
    else sol <- c(theta[2])
    gevrLik <- function(a, xp) {
      loc <- xp + (a[1]/a[2]) * (-expm1(-a[2] * log(m)))
      if (a[1] <= 0) {
        out <- .Machine$double.xmax
      }
      else {
        out <- dgevr(data, loc = loc, scale = a[1], shape = a[2], 
                     log.d = TRUE)
        out <- -sum(out)
        if (out == Inf) 
          out <- .Machine$double.xmax
      }
      out
    }
    gumLik <- function(a, xp) {
      loc <- xp + a * log(m)
      if (a <= 0) {
        out <- .Machine$double.xmax
      }
      else {
        out <- dgevr(data, loc = loc, scale = a, shape = 0, 
                     log.d = TRUE)
        out <- -sum(out)
        if (out == Inf) 
          out <- .Machine$double.xmax
      }
      out
    }
    cutoff <- qchisq(conf, 1)
    prof <- function(xp) {
      if (!z$gumbel) {
        lmax <- dgevr(data, theta[1], theta[2], theta[3], 
                      log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gevrLik, method = opt, xp = xp)
      }
      else {
        lmax <- dgevr(data, theta[1], theta[2], 0, log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gumLik, method = opt, xp = xp)
      }
      sol <- yes$par
      lci <- -yes$value
      2 * (lmax - lci) - cutoff
    }
    prof <- Vectorize(prof)
    proff <- function(xp) {
      if (!z$gumbel) {
        lmax <- dgevr(data, theta[1], theta[2], theta[3], 
                      log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gevrLik, method = opt, xp = xp)
      }
      else {
        lmax <- dgevr(data, theta[1], theta[2], 0, log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gumLik, method = opt, xp = xp)
      }
      sol <- yes$par
      lci <- -yes$value
      lci
    }
    proff <- Vectorize(proff)
    suppressWarnings(out1 <- uniroot(prof, c(est - 1e-06, est), extendInt = "downX"))
    suppressWarnings(out2 <- uniroot(prof, c(est, est + 1e-06), extendInt = "upX"))
    CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
    if (plot) {
      prof1 <- function(xp) {proff(xp)}
      suppressWarnings(curve(prof1, from = CI[1]*limi_1, to = CI[2]*limi_2, n = nub, xlab = "Return Level", ylab = ""))
      title(main = "Profile Likelihood Function")
      abline(v = est, col = "blue")
      abline(v = CI[1], col = "blue")
      abline(v = CI[2], col = "blue")
      if(length(fit$par.ests) == 2){
        
        abline(h = sum(dgevr(as.matrix(fit$data), theta[1], theta[2], 0, log.d = TRUE))-cutoff/2, lty = 2)
        
      }else{
        
        abline(h = sum(dgevr(as.matrix(fit$data), theta[1], theta[2], theta[3], log.d = TRUE))-cutoff/2, lty = 2)
        
      }
      #abline(h = sum(dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE))-cutoff/2, lty = 2)
    }
  }
  out <- list(est, CI, period, conf)
  names(out) <- c("Estimate", "CI", "Period", "ConfLevel")
  out
}
plot_prof = function(data, name = "", ys = "", ks, limi_1 = 0.9, limi_2 = 1.1){
  
  plot(data$r, data$Est, main = name, xlab = "k", ylab = ys, xlim = c(1, ks), ylim = c(min(data$Lower)*limi_1, max(data$Upper)*limi_2))
  polygon(c(rev(data$r), data$r), c(rev(data$Lower), data$Upper), col = 'grey80', border = NA)
  points(data$r, data$Est, pch = 19, col = 'black')
  lines(data$r, data$Est, lty = 'solid', col = 'black')
  lines(data$r, data$Lower, lty = 'dashed', col = 'red')
  lines(data$r, data$Upper, lty = 'dashed', col = 'red')
  
}
rlarg.diag = function(z, n = z$r){
  
  z2 = z
  z2$data = z$data[, 1]
  
  if(z$trans){
    
    for (i in 1:n) {
      rlarg.pp(c(0, 1, 0), z$data[, 1:z$r], i)
      rlarg.qq(c(0, 1, 0), z$data[, 1:z$r], i)
    }
    
  }else{
    
    par(mfcol = c(2, 2))
    gev.diag(z2)
    
    for (i in 1:n) {
      rlarg.pp(z$mle, z$data, i)
      rlarg.qq(z$mle, z$data, i)
    }
    
    par(mfcol = c(1, 1))
    
  }
  
}
rlarg.pp = function(a, dat, k){
  
  da = dat[!is.na(dat[, k]), k]
  plot((1:length(da))/length(da), rlargf(a, sort(da), k), xlab = "", ylab = "")
  title(paste("k=", k, sep = ""), cex = 0.7)
  abline(0, 1, col = 4)
  
}
rlarg.qq = function(a, dat, k){
  
  da = dat[!is.na(dat[, k]), k]
  plot(rlargq(a, 1 - (1:length(da)/(length(da) + 1)), k, da), sort(da), xlab = "", ylab = "")
  title(paste("k=", k, sep = ""), cex = 0.7)
  abline(0, 1, col = 4)
  
}
rlargq = function(a, p, k, dat){
  
  res = NULL
  for(i in 1:length(p)) {
    inter = c(min(dat) - 1e+06, max(dat) + 1e+06)
    res[i] = uniroot(rlargq2, inter, a = a, kk = k, p = p[i])$root
  }
  res
  
}
rlargq2 = function(x, a, kk, p){
  
  res = rlargf(a, x, kk) - (1 - p)
  res
  
}
rlargf = function(a, z, k){
  
  eps = 10^(-6)
  res = NULL
  
  if(abs(a[3]) < eps){
    
    tau = exp(-(z - a[1])/a[2])
    
  }else{
    
    tau = (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3])
    
  }
  
  for(i in 1:length(tau)){
    
    if(is.na(tau[i])){ 
      
      res[i] = 1
      
    }else{
      
      res[i] = exp(-tau[i])*sum(tau[i]^(0:(k - 1))/gamma(1:(k)))
      
    }
    
  }
  res
  
}
gumb.diag = function(z, gumbel = FALSE){
  
  n <- length(z$data)
  x <- (1:n)/(n + 1)
  if (z$trans) {
    oldpar <- par(mfrow = c(1, 2))
    plot(x, exp(-exp(-sort(z$data))), xlab = "Empirical", 
         ylab = "Model")
    abline(0, 1, col = 4)
    title("Residual Probability Plot")
    plot(-log(-log(x)), sort(z$data), ylab = "Empirical", 
         xlab = "Model")
    abline(0, 1, col = 4)
    title("Residual Quantile Plot (Gumbel Scale)")
  }
  else {
    
    if(!gumbel){
      
      oldpar <- par(mfrow = c(2, 2))
      gev.pp(z$mle, z$data)
      gev.qq(z$mle, z$data)
      gev.rl(z$mle, z$cov, z$data)
      gev.his(z$mle, z$data)
      
    }else{
      
      oldpar <- par(mfrow = c(2, 2))
      gev.pp(c(z$par.ests,0), z$data[,1])
      gev.qq(c(z$par.ests,0), z$data[,1])
      gev.rl(c(z$par.ests,0), dat = fitter$data[,1], fit = fitter)
      gev.his(c(z$par.ests,0), z$data[,1])
      
    }
    
  }
  par(oldpar)
  invisible()
}
gev.rl = function(a, dat, ty, fit  = NA){
  
  up = c()
  down = c()
  
  f = c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
  q = gevq(a, 1 - f)
  for(r in 1:length(f)){
    
    up = c(up, gevrRl(fit, period = 1/(1-f[r]), method = "profile", plot = FALSE)$CI[2])
    down = c(down, gevrRl(fit, period = 1/(1-f[r]), method = "profile", plot = FALSE)$CI[1])
    
  }
  
  plot(-1/log(f), q, log = "x", type = "n", xlim = c(0.1, 1000), xlab = "Return Period", ylab = "Return Level")
  title("Return Level Plot")
  lines(-1/log(f), q)
  lines(-1/log(f), up, col = 4)
  lines(-1/log(f), down, col = 4)
  points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat))
  
  # plot(1/(1-f), q, log = "x", type = "n", xlab = "Return Period", ylab = "Return Level", ylim = c(min(q,up, down)*0.9, max(q, up, down)*1.1))
  # title("Return Level Plot")
  # lines(1/(1-f), q)
  # lines(1/(1-f), up, col = 4)
  # lines(1/(1-f), down, col = 4)
  # points(1/(1 - (1:length(dat))/(length(dat) + 1)), sort(dat))
  
  tt = gevrRl(fit, period = ty, method = "profile", plot = FALSE)
  
  segments(ty, 0, ty, tt$Estimate, col = 'midnightblue', lty = 6)
  segments(0.0000001, tt$Estimate, ty, tt$Estimate , col ='midnightblue', lty = 6)
  
}
gpd.prof = function (z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100, mrt = "BFGS") {
  cat("If routine fails, try changing plotting interval", fill = TRUE)
  xdat <- z$data
  u <- z$threshold
  la <- z$rate
  v <- numeric(nint)
  x <- seq(xlow, xup, length = nint)
  m <- m * npy
  sol <- z$mle[2]
  gpd.plik <- function(a) {
    if (m != Inf) 
      sc <- (a * (xp - u))/((m * la)^a - 1)
    else sc <- (u - xp)/a
    if (abs(a) < 10^(-4)) 
      l <- length(xdat) * log(sc) + sum(xdat - u)/sc
    else {
      y <- (xdat - u)/sc
      y <- 1 + a * y
      if (any(y <= 0) || sc <= 0) 
        l <- 10^6
      else l <- length(xdat) * log(sc) + sum(log(y)) * 
        (1/a + 1)
    }
    l
  }
  for (i in 1:nint) {
    xp <- x[i]
    opt <- optim(sol, gpd.plik, method = mrt)
    sol <- opt$par
    v[i] <- opt$value
  }
  plot(x, -v, type = "l", xlab = "Return Level", ylab = "Profile Log-likelihood")
  ma <- -z$nllh
  #abline(h = ma)
  abline(h = ma - 0.5 * qchisq(conf, 1))
  invisible()
  
  return(cbind(x,-v))
}
gpd.rl = function(a, u, la, n, npy = 365, mat, dat, xdat){
  #
  # function called by gpd.diag
  # produces return level curve and 95% confidence intervals
  # for fitted gpd model
  a <- c(la, a)
  eps <- 1e-006
  a1 <- a
  a2 <- a
  a3 <- a
  a1[1] <- a[1] + eps
  a2[2] <- a[2] + eps
  a3[3] <- a[3] + eps
  jj <- seq(-1, 3.75 + log10(npy), by = 0.1)
  m <- c(1/la, 10^jj)
  q <- gpdq2(a[2:3], u, la, m)
  # d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
  # d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
  # d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
  # d <- cbind(d1, d2, d3)
  d <- t( gpd.rl.gradient( a=a, m=m))
  mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1, 2], 0, 
                  mat[2, 1], mat[2, 2]), ncol = 3)
  v <- apply(d, 1, q.form, m = mat)
  plot(m/npy, q, log = "x", type = "n", xlim = c(0.1, max(m)/npy), ylim
       = c(u, max(xdat, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1])), 
       xlab = "Return period (years)", ylab = "Return level", main = 
         "Return Level Plot")
  lines(m[q > u - 1]/npy, q[q > u - 1])
  lines(m[q > u - 1]/npy, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1], col
        = 4)
  lines(m[q > u - 1]/npy, q[q > u - 1] - 1.96 * sqrt(v)[q > u - 1], col
        = 4)
  nl <- n - length(dat) + 1
  sdat <- sort(xdat)
  points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u])	
  #	points(1/(1 - (1:n)/(n + 1))/npy, 
  #		sort(xdat))
  #	abline(h = u, col = 3)
}
sum_plot = function(data_set, title = NULL, xlabel = NULL, ylabel = NULL, type = "den", l_pos = "none"){
  
  if(type == "den"){
    
    ggplot(data_set, aes(x = `val`, y = `Month`, fill = ..x..)) +
      geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
      scale_fill_viridis(name = xlabel, option = "C", alpha = 1) + xlab(xlabel) + ylab(ylabel) +
      ggtitle(title) + 
      theme_ipsum() +
      theme(
        legend.position=l_pos,
        panel.spacing = unit(0.4, "lines"),
        strip.text.x = element_text(size = 8)
      ) + theme(plot.title=element_text(hjust=0.5), axis.title.x = element_text(hjust=0.5), axis.title.y = element_text(hjust=0.5)) 
    
  }else{
    
    ggplot(data_set, aes(x = `val`, y = `Month`,  fill = `val`)) +
      geom_density_ridges(alpha=0.6, stat="binline", bins = 25) +
      scale_fill_viridis(discrete=TRUE) +
      scale_color_viridis(discrete=TRUE) +
      theme_ipsum() +
      theme(
        legend.position=l_pos,
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +  theme(plot.title=element_text(hjust=0.5), axis.title.x = element_text(hjust=0.5), axis.title.y = element_text(hjust=0.5)) +
      ggtitle(title) + xlab(xlabel) + ylab("Assigned Probability (%)")
    
  }
  
}
box_plot = function(data_set, title = NULL, xlabel = NULL, ylabel = NULL, same = TRUE, l_pos = "none", leg_Title = NULL){
  
  if(same == TRUE){
    
    ggplot(data_set, aes(Month, val, fill = Month))+ geom_boxplot() +
      theme(legend.position = l_pos, plot.title=element_text(hjust=0.5), axis.title.x = element_text(hjust=0.5), axis.title.y = element_text(hjust=0.5)) +
      ggtitle(title) + xlab(xlabel) + ylab(ylabel) + stat_boxplot(geom ='errorbar') + geom_boxplot()
    
  }else{
    
    ggplot(data_set, aes(Month, val, fill = vel))+ geom_boxplot() +
      theme(legend.position = l_pos, plot.title=element_text(hjust=0.5), axis.title.x = element_text(hjust=0.5), axis.title.y = element_text(hjust=0.5)) +
      ggtitle(title) + xlab(xlabel) + ylab(ylabel) + scale_fill_discrete(name = leg_Title) + stat_boxplot(geom ='errorbar') + geom_boxplot()
    
  }
  
}
time2season = function(x, out.fmt="months", type="default") {
  
  # Checking that 'class(x)==Date'
  #if ( ( !( class(x) %in% c("Date", "POSIXct", "POSIXt") ) ) && TRUE ) 
  if (!( is(x, "Date") | is(x, "POSIXct") | is(x, "POSIXt") )) 
    stop("Invalid argument: 'x' must be in c('Date', 'POSIXct', 'POSIXt') !")
  
  # Checking the class of out.fmt
  if (is.na(match(out.fmt, c("seasons", "months") ) ) )
    stop("Invalid argument: 'out.fmt' must be in c('seasons', 'months')")
  
  # Checking that the user provied a valid value for 'type'   
  valid.types <- c("default", "FrenchPolynesia")    
  if (length(which(!is.na(match(type, valid.types )))) <= 0)  
    stop("Invalid argument: 'type' must be in c('default', 'FrenchPolynesia')")
  
  ####################
  months <- format(x, "%m")
  
  if (type=="default") {
    winter <- which( months %in% c("12", "01", "02") )
    spring <- which( months %in% c("03", "04", "05") )
    summer <- which( months %in% c("06", "07", "08") )
    autumm <- which( months %in% c("09", "10", "11") ) 
  } else if (type=="FrenchPolynesia") {
    winter <- which( months %in% c("12", "01", "02", "03") )
    spring <- which( months %in% c("04", "05") )
    summer <- which( months %in% c("06", "07", "08", "09") )
    autumm <- which( months %in% c("10", "11") ) 
  } # ELSE end
  
  # Creation of the output, with the same length of the 'x' input
  seasons <- rep(NA, length(x))
  
  if (out.fmt == "seasons") {
    
    seasons[winter] <- "winter"
    seasons[spring] <- "spring"
    seasons[summer] <- "summer"
    seasons[autumm] <- "autumm"
    
  } else { # out.fmt == "months"
    
    if (type=="default") {
      seasons[winter] <- "DJF"
      seasons[spring] <- "MAM"
      seasons[summer] <- "JJA"
      seasons[autumm] <- "SON"
    } else  if (type=="FrenchPolynesia") {
      seasons[winter] <- "DJFM"
      seasons[spring] <- "AM"
      seasons[summer] <- "JJAS"
      seasons[autumm] <- "ON"
    } # IF end
    
  } # IF end
  
  return(seasons)
  
}

#------------------ Done ~ Formatting Data / Providing Stats  -----------------

#To format the data
#data[,"Date"] = as.POSIXct(data[,"Date"], tz = "", format = "%d/%m/%Y %H:%M")
data[,"Date"] = as.Date(data[,"Date"], format = "%Y-%m-%d")

#Summary of data
glimpse(data)
stat.desc(data[,c(3,4)])
summary(data[,3])
summary(data[,4])

#Define the time series (as.numeric(format(data[1,1], "%j"))))
ts_1 = ts(data[,3], start = c(1969, as.numeric(format(data[1,"Date"],format = "%m"))), frequency = 12)
ts_2 = ts(data[,4], start = c(1969, as.numeric(format(data[1,"Date"],format = "%m"))), frequency = 12)

#Visual Graphics
bin = hexbin(data[,3], data[,4], xbins=50)
plot(ts_1, ylab = names(data)[3])
plot(data[,3], ylab = names(data)[3])
plot(ts_2, ylab = names(data)[4])
plot(data[,4], ylab = names(data)[4])
plot(data[,3], data[,4], xlab = names(data)[3], ylab = names(data)[4], main = paste0(names(data)[3], " vs ", names(data)[4]))
plot(bin, xlab = names(data)[3], ylab = names(data)[4], main = "Hexagonal Binning")

#All data
df = data

#Keep the months (October - March)
P1_df = data[which(month(data[,"Date"]) %in% c(10,11,12,1,2,3)),]

#Keep the months (April - September)
P2_df = data[which(month(data[,"Date"]) %in% c(4,5,6,7,8,9)),]

#Obtain summary statistics
summary(df)
summary(P1_df)
summary(P2_df)

#Grouping visuals
plot(data[,3], data[,4], xlab = "M-MDD", ylab = "SR-TWE")
ggplot(data, aes(x = MRS, y = MRD) ) +  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) +  scale_fill_viridis() +  theme(legend.position='right') + guides(fill=guide_legend(title="Density"))
p = ggplot(data, aes(x = MRS, y = MRD)) + geom_point() + theme_bw()
ggMarginal(p, type="boxplot", fill = "#0099f8")
#ggMarginal(p, type="histogram", fill = "#0099f8")

#Select the data to work with (df/P1_df/P2_df) and the variables (i.e. Date (2) and Mean Discharge (3) or Summation Total Water (4))
seldata = P1_df[,c("Date","MRD", "hydroYear")]

#------------------ Done ~ Univariate Model: Block Maxima Method (BM M) -----------------

summary(seldata)

#Idea of the distribution and outliers
hist(seldata[,2], breaks = 20, col = "green", density = 20, main = bquote(""*underline(.(paste0(paste("Histogram of", firstup(colnames(seldata)[2])))))))
boxplot(seldata[,2])
title(main = bquote(""*underline(.(paste0(paste("Histogram of", firstup(colnames(seldata)[2])))))), cex.main = 1.5)

newdata = data.frame("Dates" = as.Date(character()), "Variable" = as.numeric())
names(newdata)[2] = names(seldata)[2]
for(i in 1:length(unique(seldata$hydroYear))){
  
  newdata[i,2] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2]))]
  newdata[i,1] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),1][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2]))]
  
}

#Plot of the defined new time series
plot(seldata[,c(1,2)], main = bquote(""*underline(.(paste0(paste("Block Maxima:", firstup(colnames(seldata)[2])))))), ylab = firstup(colnames(seldata)[2]), xlab  = "Time Index")
points(newdata, col = "red")
plot(newdata, main = "", ylab = "Annual maxima", xlab  = "Time Index")
title(main = bquote(""*underline(.(paste0(paste("Block Maxima:", firstup(colnames(seldata)[2])))))), cex.main = 1.5)

#Fit the parametric GEV distribution 
#fitt = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)
fitt = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE)
fiter = extRemes::fevd(newdata[,2], method = "MLE", type = "GEV")
fitos = ismev::gev.fit(newdata[,2])

#Test the fit
gnfit(newdata[,2], dist = "gev", pr = gev.fit(newdata[,2])$mle)

#To test if the model with less parameters (Gumbel) is better
lr.test(extRemes::fevd(newdata[,2], method = "MLE", type = "Gumbel"), extRemes::fevd(newdata[,2], method = "MLE", type = "GEV"), alpha = 0.05)

summary(fiter)
confint(fitt)
ci(fiter, alpha = 0.05, type = "parameter")

#Profile log-likelihood confidence intervals 
#limits = logLik(fitt)[1] -0.5*qchisq(p = 0.95, df = 1)
#The limits are obtained by Lp(mu) > Lmax greater than 0.5 * q
limi = print(suppressWarnings(confint(profile(fitt, mesh = fitt$std.err/20), level = 0.95)))
suppressWarnings(confint(profile(fitt, mesh = fitt$std.err/20), level = 0.95))

#Profile log-likelihood plots 
par(mfrow = c(1,3))
suppressWarnings(plot(profile(fitt, mesh = fitt$std.err/55), which = names(fitt$estimate)[1], ci = c(0.95)))
suppressWarnings(plot(profile(fitt), which = names(fitt$estimate)[2], ci = c(0.95)))
suppressWarnings(plot(profile(fitt, mesh = fitt$std.err/50), which = names(fitt$estimate)[3], ci = c(0.95), xlim = c(-0.9,0.7)))
par(mfrow = c(1,1))

#To create the profile log-likelihoods separately
plot(profile(fitt, mesh = fitt$std.err/55), ci = c(0.95), which = names(fitt$estimate)[1], main = "")
title(main = bquote(""*underline(.(paste0("Profile Log-likelihood of Location")))), cex.main = 2)
abline(v = limi[1], lty=2)
abline(v = limi[4], lty=2)
plot(profile(fitt, mesh = fitt$std.err/20), ci = c(0.95), which = names(fitt$estimate)[2], main = "")
title(main = bquote(""*underline(.(paste0("Profile Log-likelihood of Scale")))), cex.main = 2)
abline(v = limi[2], lty=2)
abline(v = limi[5], lty=2)
plot(profile(fitt, mesh = fitt$std.err/50), ci = c(0.95), which = names(fitt$estimate)[3], main = "")
title(main = bquote(""*underline(.(paste0("Profile Log-likelihood of Shape")))), cex.main = 2)
abline(v = limi[3], lty=2)
abline(v = limi[6], lty=2)

#Diagnostic plots type 1
plot(fitt)
gev.diag(gev.fit(newdata[,2]))

#Diagnostic plots type 2
par(mfrow = c(2 ,2))
plot(fiter, type = c("probprob"), main = bquote(""*underline(.(paste0("Probability Plot")))), cex.main = 2)
plot(fiter, type = c("qq"), main = bquote(""*underline(.(paste0("Q-Q Plot")))), cex.main = 2)
plot(fiter, type = c("rl"), main = bquote(""*underline(.(paste0("Return level Plot")))), cex.main = 2)
plot(fiter, type = c("density"), main = bquote(""*underline(.(paste0("Density Plot")))), cex.main = 2)
par(mfrow = c(1 ,1))

#Checking requirements
#plot(logspline(newdata[,2]), add = T, col = "blue")
#curve(evd::dgev(x, loc = fitt$estimate[1], scale
#plot(logspline(y), add = T, col = "red", lty = 2)
y = revd(10000,loc = fitt$estimate[1], scale = fitt$estimate[2], shape = fitt$estimate[3])
hist(newdata[,2], probability = TRUE, xlab = "Observations", ylim = c(0,0.00055), main = bquote(""*underline(.(paste0("Density Plot")))), font.main = 1, cex.main = 2)
lines(density(newdata[,2]), col = "blue")
lines(density(y), col = "red", lty = 2)
legend("topright", legend=c("Empirical", "Fitted"), col=c("blue", "red"), lty = 1:2, cex = 1)
box(which = "plot", lty = "solid")

#Finding 0.9, 0.95 and 0.99 quantiles for the fitted GEV distribution
vec = c(0.90, 0.95, 0.99)
for(k in 1:length(vec)){print(evd::qgev(vec[k],loc = fitt$estimate[1], scale = fitt$estimate[2], shape = fitt$estimate[3])) }

#Profile log-likelihood plot for the specified return level
ret = 200
return.level(fiter, conf = 0.05, return.period = ret, do.ci = TRUE)
tt = as.numeric(return.level(fiter, conf = 0.05, return.period = ret))

#Plot Return Lev vs Return Periods (tpye = probprob/qq/qq2/density/trace)
#plot(fiter, type = "rl", rperiods = c(5, 10, 20, 30, 40, 50), main = "", pch = 16)
plot(fiter, type = "rl", main = "", pch = 16, rperiods = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 250))
title(main = bquote(""*underline(.(paste0("Return Level Curve")))), cex.main = 2)
segments(ret, 0, ret, tt, col = 'midnightblue', lty = 6)
segments(0.0000001, tt, ret, tt, col ='midnightblue', lty = 6)

#Profile log-likelihood plot for the specified return level
gev.prof(gev.fit(newdata[,2]), m = ret, 4900, 5880)
segments(tt, (-gev.fit(newdata[,2], show = FALSE)$nllh - 10), tt, -gev.fit(newdata[,2], show = FALSE)$nllh, col = 'midnightblue', lty = 6)

#title(main = bquote(""*underline(.(paste0(ret, " year Return level Profile Log-likelihood")))), cex.main = 2)
profliker(fiter, type ="return.level", return.period = ret, xrange = c(5400, 7000), main ="")
ci_rl = ci(fiter , alpha = 0.05 , type = "return.level", return.period = ret , method ="proflik", xrange = c(5400, 7000))
print(ci_rl)
abline(v = ci_rl[1])
abline(v = ci_rl[3])
abline(v = ci_rl[2])
abline(h = (-gev.fit(newdata[,2], show = FALSE)[["nllh"]] - qchisq(0.95, df=1)/2) - 0.08, lty = 2)

#------------------ Done ~ Univariate Model: Generalized Block Maxima Method (KLOS M) -----------------

#To create a dataframe of the largest order statistics
k = as.numeric(dlgInput(message = "Please enter the number of extreme observations to consider", default = 3, gui = .GUI)$res) 
dfe = data.frame(matrix(ncol = k, nrow = length(unique(seldata$hydroYear))))
dff = data.frame(matrix(ncol = 2, nrow = k*length(unique(seldata$hydroYear))))
for(i in 1:(length(unique(seldata$hydroYear)))){
  
  #Each row represents the unique years and each column represents the max 
  dfe[i,] = seldata[kit::topn(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2], k) + which(seldata$hydroYear == unique(seldata$hydroYear)[i])[1] - 1, 2]
  
  dff[(1:k)+(i-1)*k,2] =  as.numeric(t(dfe[i,]))
  dff[(1:k)+(i-1)*k,1] =  as.numeric(unique(format(seldata$Date, format = "%Y"))[i])
  
}

#Plot of the defined new time series
plot(seldata[,c(1,2)], main = bquote(""*underline(.(paste0(paste("KLOS:", firstup(colnames(seldata)[2])))))), ylab = firstup(colnames(seldata)[2]), xlab  = "Time Index")
points(x = seldata[sapply(dff[,2], function(x) which(seldata[,2]==x)[1]),1], y = seldata[sapply(dff[,2], function(x) which(seldata[,2]==x)[1]),2], col = "red")
abline(v = seq(from = seldata[1,1],length.out = length(unique(dff[,1])), by = 'year'), col = "midnightblue", lty=2)

plot(dff, xlab = "Time Periods", ylab = "Largest Observations", main = names(seldata)[2])
abline(lm(dff[,2]~dff[,1]))
summary(lm(dff[,2]~dff[,1]))

#Fitting the model for different values of k
MLE  = data.frame(matrix(ncol = 3, nrow = k))
SE   = data.frame(matrix(ncol = 3, nrow = k))
NLLH = array(0, dim = k)
for(i in 1:k){
  
  RLOS = ismev::rlarg.fit(dfe, r = i , show = FALSE)
  MLE[i, ] = RLOS$mle
  SE[i, ] = RLOS$se
  NLLH[i]  = RLOS$nllh
  
} 

#Calculation of AIC and BIC
AIC = 6 + 2*NLLH
BIC = ln(nrow(dff)) + 2*NLLH

#Plotting SE pattern for the three parameters
par(mfrow = c(1,3))
plot(SE[ ,1], type = 'l', ylab = 'SE', xlab = 'k value', main = "Standard error for location parameter", xlim = c(1,k))
plot(SE[ ,2], type = 'l', ylab = 'SE', xlab = 'k value', main = "Standard error for scale parameter", xlim = c(1,k))
plot(SE[ ,3], type = 'l', ylab = 'SE', xlab = 'k value', main = "Standard error for shape parameter", xlim = c(1,k))

#Test to obtain the optimal threshold (k)
gevrSeqTests(dfe, method = "ed")
gevrSeqTests(dfe, method = "pbscore", bootnum = 100)
gevrSeqTests(dfe, method = "multscore", bootnum = 10000)

thh = as.numeric(dlgInput(message = "Please enter the return level to consider", default = 20, gui = .GUI)$res) 
vec = c("", "_loc", "_scale", "_shape")
for(h in 1:4){
  
  assign(paste0("result", vec[h]), matrix(0, k, 4))
  
}
for(j in 1:k){
  
  fitq = gevrFit(dfe[,1:j], method = "mle")
  
  y2 = gevrRl(fitq, thh, conf = 0.95, method = "profile", plot = FALSE)
  y2_loc = gevrProfLoc(fitq, plot = FALSE)
  y2_scale = gevrProfScale(fitq, plot = FALSE)
  y2_shape = gevrProfShape(fitq, plot = FALSE)
  
  result[j, 1] = j
  result[j, 2] = y2$Estimate
  result[j, 3:4] = y2$CI 
  result_loc[j, 1] = j
  result_loc[j, 2] = y2_loc$Estimate
  result_loc[j, 3:4] = y2_loc$CI 
  result_scale[j, 1] = j
  result_scale[j, 2] = y2_scale$Estimate
  result_scale[j, 3:4] = y2_scale$CI 
  result_shape[j, 1] = j
  result_shape[j, 2] = y2_shape$Estimate
  result_shape[j, 3:4] = y2_shape$CI 
  
}

colnames(result) = c("r", "Est", "Lower", "Upper")
result = as.data.frame(result)
colnames(result_loc) = c("r", "Est", "Lower", "Upper")
result_loc = as.data.frame(result_loc)
colnames(result_scale) = c("r", "Est", "Lower", "Upper")
result_scale = as.data.frame(result_scale)
colnames(result_shape) = c("r", "Est", "Lower", "Upper")
result_shape = as.data.frame(result_shape)

#Plots of the profile likelihoods
plot_prof(result, name = "Profile Likelihood", ys = paste0(thh, " Year Return Level"), ks = k)
plot_prof(result_loc, name = "Profile Likelihood", ys = "Location", ks = k)
plot_prof(result_scale, name = "Profile Likelihood", ys = "Scale", ks = k)
plot_prof(result_shape, name = "Profile Likelihood", ys = "Shape", ks = k, limi_1 = 1.1, limi_2 = 1.1)

#To obtain the threshold
th = as.numeric(dlgInput(message = "Please enter the number of extreme observations to consider", default = 3, gui = .GUI)$res) 

#To plot the profile likelihoods for the selected threshold
#fit = gevrFit(dfe[,1:th], method = "mle", gumbel = TRUE)
fit = gevrFit(dfe[,1:th], method = "mle")

#Test the fit
#gnfit(as.numeric(unlist(dfe)), dist = "gev", pr = c(as.numeric(gevrFit(dfe[,1:k], method = "mle", gumbel = TRUE)$par.ests),0))
KLOS_BEST = rlarg.fit(dfe, r = th, show = TRUE)
gnfit(as.numeric(unlist(dfe[, 1:th])), dist = "gev", pr = KLOS_BEST$mle)

#To test if the model with less parameters (Gumbel) is better
lr.test(gevrFit(dfe[,1:k], method = "mle", gumbel = TRUE)$nllh.final, gevrFit(dfe[,1:k], method = "mle", gumbel = FALSE)$nllh.final, alpha = 0.05)

#Profile Likelihood CI
#eva::gevrDiag(fit, method = "profile")
par(mfrow = c(1,3))
gevrProfLoc(fit, nub = 30, limi_1 = 0.99, limi_2 = 1.01)
gevrProfScale(fit, nub = 20, limi_1 = 0.992, limi_2 = 1.02)
gevrProfShape(fit, nub = 100, limi_1 = 1.02, limi_2 = 1.02)
par(mfrow = c(1,1))

#Fitting the model for the selected value of k
k = th
KLOS_BEST = rlarg.fit(dfe, r = k, show = TRUE)
fit = gevrFit(dfe[,1:k], method = "mle")

#Diagnostic plots
#fitter = c(fit, trans = list(FALSE))
#gumb.diag(fitter, gumbel = TRUE)
rlarg.diag(KLOS_BEST)
par(mfrow = c(1,1))

#To plot the return period graph and obtain the profile log likelihood of the return level
gev.rl(KLOS_BEST$mle, dat = KLOS_BEST$data[,1], ty = thh, fit = fit)
gevrRl(fit, period = thh, method = "profile", limi_1 = 0.95, limi_2 = 1.1)

#------------------ Done ~ Univariate Model: Peak over Threshold Method (POT M)-----------------

#First part conditioned on wind
#extRemes::mrlplot(seldata[,2], nint = 500, main = "", xlim = c(800,5400))
mrl = POT::mrlplot(seldata[,2], main = "",  xlim = c(800,5400), ylim = c(5,1700), lwd = c(0.8, 1, 0.8), nt = 1000)
title(main = bquote(""*underline(.(paste0("Mean Residual Life Plot")))), line = 3.5, cex.main = 1.5)

#To add the third axis 
trs = seq(from = 0, to = round(max(mrl$x)), by = 100)
set = numeric(length(trs))
for(i in 1:length(set)){set[i] = length(seldata[,2][seldata[,2] > trs[i]])}
axis(side = 3, at = trs, labels = set)
mtext(text = "Number of Observations", side = 3, line = 2)

#To plot the straight line
rag = dlgInput(message = "Please enter the range of threshold", default = "2700,3200", gui = .GUI)$res
rag_1 = as.numeric(trimws(str_split_fixed(rag, ",", 2))[1]) 
rag_2 = as.numeric(trimws(str_split_fixed(rag, ",", 2))[2])
u_thre = seq(from = rag_1, to = rag_2, by = 0.5)
emp = numeric(length(u_thre))

#Evaluating the mean of X-u for all thresholds in u
for(i in 1:length(u_thre)){emp[i] = mean(seldata[,2][seldata[,2] > u_thre[i]] - u_thre[i])}
cut = data.frame("Threshold" = u_thre, "Mean_Excess" = emp)
ablineclip(lm(cut$Mean_Excess~cut$Threshold), x1 = rag_1, x2 = rag_2, col = "blue")
abline(v = rag_1, col = "blue", lty = 3)
abline(v = rag_2, col = "blue", lty = 3)

#To obtain the quantile and the total number of remaining population
#summary(new_data[,4])
#quantile(new_data[,4], probs = c(0,0.25,0.5,0.7,0.75,0.8,0.9,0.95,0.99,0.993))
paste0(round(ecdf(seldata[,2])(rag_1)*100,3), "%")
length(which(seldata[,2] > rag_1))
length(which(seldata[,2] > rag_2))

#Fitting the GPD Model Over a Range of Thresholds
ismev::gpd.fitrange(seldata[,2], (rag_1 - 100), (rag_2 + 150))
title(main = bquote(""*underline(.(paste0("Parameter Stability Plots")))), cex.main = 1.5)

#Fitting gpd using a threshold determined by the above methods
thrs = as.numeric(dlgInput(message = "Please enter threshold value", default = paste0(rag_1), gui = .GUI)$res)
fitt = mev::fit.gpd(seldata[,2], threshold = thrs)
fiter = fevd(seldata[,2], method = "MLE", type = "GP", threshold = thrs, span = length(unique(seldata$hydroYear)))
fito = evd::fpot(seldata[,2], threshold = thrs, model = "gpd", npp = 6)
gpd_fit = gpd.fit(seldata[,2], threshold = thrs, npy = 6)

#Calculation of AIC and BIC
AIC = 4 + 2*gpd_fit$nllh
BIC = ln(gpd_fit$nexc) + 2*gpd_fit$nllh

#Test the fit
gnfit(as.vector(seldata[,2]), dist = "gpd", pr = gpd_fit$mle, threshold = rag_1)
lr.test(mev::fit.gpd(seldata[,2], threshold = thrs, fpar = list("shape" = 0))$nllh, mev::fit.gpd(seldata[,2], threshold = thrs)$nllh, alpha = 0.05)
par(mfrow = c(1,1))

#Plot the remaining observations
plot(seldata[seldata[,2] > thrs,2], main = "", ylab = firstup(colnames(seldata)[1:length(colnames(seldata))][2]), xlab  = "Time Index")
title(main = bquote(""*underline(.(paste0("Excess Observations")))), cex.main = 1.5)

#Plot of the defined new time series
plot(seldata[,c(1,2)], main = bquote(""*underline(.(paste0(paste("POT:", firstup(colnames(seldata)[2])))))), ylab = firstup(colnames(seldata)[2]), xlab  = "Time Index")
points(x = seldata[sapply(seldata[seldata[,2] > thrs,2], function(x) which(seldata[,2]==x)[1]),1], y = seldata[sapply(seldata[seldata[,2] > thrs,2], function(x) which(seldata[,2]==x)[1]),2], col = "red")
abline(v = seq(from = seldata[1,1],length.out = length(unique(seldata[,1])), by = 'year'), col = "midnightblue", lty=2)
abline(h = thrs, col = "green", lty=2)

#95% Wald confidence interval
confint(fitt)
ci(fiter, confidence = 0.95, type = "parameter")

#Profile log-likelihood confidence intervals 
#limits = logLik(fitt)[1] -0.5*qchisq(p = 0.95, df = 1)
#The limits are obtained by Lp(mu) > Lmax greater than 0.5 * q
limi = print(suppressWarnings(confint(profile(fito, mesh = fito$std.err/55), level = 0.95)))
suppressWarnings(confint(profile(fito, mesh = fito$std.err/50), level = 0.95))

#Profile log-likelihood plots 
par(mfrow = c(1,2))
suppressWarnings(plot(profile(fito, mesh = fito$std.err/55), ci = c(0.95)))
par(mfrow = c(1,1))

#To create the profile log-likelihoods separately
par(mfrow = c(1,2))
plot(profile(fito, mesh = fito$std.err/55), ci = c(0.95), which = names(fito$estimate)[1], main = "")
title(main = bquote(""*underline(.(paste0("Profile Log-likelihood of Scale")))), cex.main = 1.1)
abline(v = limi[1])
abline(v = limi[3])
plot(profile(fito, mesh = fito$std.err/55), ci = c(0.95), which = names(fito$estimate)[2], main = "")
title(main = bquote(""*underline(.(paste0("Profile Log-likelihood of Shape")))), cex.main = 1.1)
abline(v = limi[2])
abline(v = limi[4])
par(mfrow = c(1,1))

#Diagnostic plots type 2
par(mfrow = c(2 ,2))
plot(fiter, type = c("probprob"), main = bquote(""*underline(.(paste0("Probability Plot")))), cex.main = 1.5)
plot(fiter, type = c("qq"), main = bquote(""*underline(.(paste0("Q-Q Plot")))), cex.main = 1.5)
plot(fiter, type = c("rl"), main = bquote(""*underline(.(paste0("Return level Plot")))), cex.main = 1.5)
plot(fiter, type = c("density"), main = bquote(""*underline(.(paste0("Density Plot")))), cex.main = 1.5)
par(mfrow = c(1 ,1))

#Diagnostic plots
plot(fitt)
gpd.diag(gpd_fit)
plot(fiter)

par(mfrow = c(1,1))
#plot(logspline(new_data[ ,4][new_data[ ,4] > thrs]), add = T, col = "blue")
#plot(logspline(y + thrs), add = T, col = "red", lty = 2)
#curve(evd::dgpd(x, loc = 0, scale = fitt$estimate[1], shape = fitt$estimate[2], log = FALSE), add=TRUE)
y = extRemes::revd(10000, scale = gpd_fit$mle[1], shape = gpd_fit$mle[2], type = "GP")
hist(seldata[ ,2][seldata[ ,2] > thrs], ylim = c(0,0.14), probability = TRUE, xlab = "Observations", main = bquote(""*underline(.(paste0("Density Plot")))), font.main = 1, cex.main = 1.5)
lines(density(seldata[ ,2][seldata[ ,2] > thrs]), col = "blue")
lines(density(y + thrs), col = "red")
legend("topright", legend=c("Empirical", "Fitted"), col=c("blue", "red"), lty = 1:2, cex = 1)
box(which = "plot", lty = "solid")

#Finding 0.9, 0.95 and 0.99 quantiles for the fitted GEV distribution
vec = c(0.90, 0.95, 0.99)
for(k in 1:length(vec)){print(evd::qgpd(vec[k], scale = gpd_fit$mle[1], shape = gpd_fit$mle[2]))}

#The return level for the specified return level
rev = 20
jj = rev*gpd_fit$npy
tt = gpdq2(gpd_fit$mle, gpd_fit$threshold, gpd_fit$rate, m = jj)

#Plot Return Lev vs Return Periods (tpye = probprob/qq/qq2/density/trace)
#plot(fiter, type = "rl", main = bquote(""*underline(.(paste0("POT Method")))), pch = 16)
gpd.rl(gpd_fit$mle, gpd_fit$threshold, gpd_fit$rate, gpd_fit$n, gpd_fit$npy, gpd_fit$cov, gpd_fit$data, gpd_fit$xdata)
segments(rev, 0, rev, tt, col = 'midnightblue',lty = 6)
segments(0.00001, tt, rev, tt, col = 'midnightblue', lty = 6)

#Profile log-likelihood plot for the specified return level
x = gpd.prof(gpd_fit, m = rev, 5050, 5800, npy = 6, nint = 100, mrt = "Nelder-Mead")
title(main = bquote(""*underline(.(paste0("Profile Likelihood Function")))), cex.main = 1.4)
abline(v = tt, col = 'midnightblue', lty = 2)

fy = x[-which(x[,1]>=tt),]
ly = x[which(x[,1]>=tt),]
fy = as.data.frame(fy)
ly = as.data.frame(ly)
names(fy) = c("rt","nllh")
names(ly) = c("rt","nllh")

fy[,2] = fy[,2]-(-gpd_fit$nllh - 0.5 * qchisq(0.95, 1))
ly[,2] = ly[,2]-(-gpd_fit$nllh - 0.5 * qchisq(0.95, 1))

fy = fy[order(abs(fy$nllh)),]
ly = ly[order(abs(ly$nllh)),]

abline(v = ly[1,1], col = 'midnightblue', lty = 2)
abline(v = fy[1,1], col = 'midnightblue', lty = 2)

print(ly[1,1])
print(fy[1,1])
