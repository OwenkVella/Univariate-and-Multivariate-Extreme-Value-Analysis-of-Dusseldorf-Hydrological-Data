#Required Packages
#install_github("cran/EcoHydRology")
#library("EcoHydRology")
#install.packages(c("tcltk", "lubridate", "dplyr", "tidyr", "pastecs", "extRemes", "svDialogs", "fExtremes", "evd", "copula", "RcppRoll"))
#install.packages(c("remotes", "zoo", "ismev", "VGAM", "SciViews", "lmomco", "copBasic", "hexbin", "scatterplot3d", "rgl", "graphics", "gumbel", "SuppDists", "openxlsx"))
#install.packages(c("ggplot2", "viridis", "ggExtra", "gridExtra", "kit", "eva", "VineCopula", "POT", "stringr", "plotrix", "mev", "logspline"))
lapply(c("tcltk", "lubridate", "dplyr", "tidyr", "pastecs", "extRemes", "svDialogs", "fExtremes", "evd", "copula", "RcppRoll"), require, character.only = TRUE)
lapply(c("ismev", "VGAM", "SciViews", "lmomco", "copBasic", "hexbin", "scatterplot3d", "rgl", "graphics", "gumbel", "SuppDists", "openxlsx"), require, character.only = TRUE)
lapply(c("remotes", "ggplot2", "viridis","ggExtra", "gridExtra", "kit", "eva", "zoo", "VineCopula", "POT", "stringr", "plotrix", "mev", "logspline"), require, character.only = TRUE)
remotes::install_github("AlexanderRitz/copR")

#To clean
rm(list = ls())
graphics.off()

#Functions
franktau = function(tau) {
  sign <- 1
  if (tau < 0) {
    sign <- -1
    tau <- -tau
  }
  if (abs(tau) > 0.99999) {
    return(Inf * sign)
  }
  thet <- function(x) {
    tau - (1 - 4 / x + 4 / x * debye(x))
  }
  rooted <- stats::uniroot(
    thet,
    lower = 0 + .Machine$double.eps ^ (1 / 2),
    upper = 5e+05,
    tol = .Machine$double.eps ^ (1 / 2)
  )
  theta <- rooted$root
  return(sign * theta)
}
debye = function(x) {
  fun <- function(x, lambda) {
    x ^ lambda / (exp(x) - 1)
  }
  up = abs(x)
  if (x == 0) {
    D = 1
  } else {
    int = stats::integrate(
      f = fun,
      lower = 0,
      upper  = up,
      lambda = 1
    )
    D = int[[1]] / up
  }
  if (x < 0) {
    D = D + up / (2)
  }
  
  return(D)
}

#Loading data (Summary_Final.csv)
#data = read.csv(file.choose())
data = read.csv(tk_choose.files(default = "", caption = "Select csv file"))

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
ggplot(data, aes(x = MRS, y = MRD) ) +  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) +  scale_fill_viridis() +  theme(legend.position='right') + guides(fill=guide_legend(title="Density"))
p = ggplot(data, aes(x = MRS, y = MRD)) + geom_point() + theme_bw()
ggMarginal(p, type="boxplot", fill = "#0099f8")
#ggMarginal(p, type="histogram", fill = "#0099f8")

#Select the data to work with (df/P1_df/P2_df) and the variables (i.e. Date (2) and Mean Discharge (3) or Summation Total Water (4))
seldata = P2_df[,c(2,3,4,1)]

#------------------ Done ~ Fitting Distribution  -----------------

summary(seldata)

#To obtain the annual maxima or the excess (need to amend accordingly)
newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
names(newdata) = names(seldata)[c(2,3)]
for(i in 1:length(unique(seldata$hydroYear))){
  
  newdata[i,1] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2]))]
  newdata[i,2] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3]))]
  
}

#Fitting the GEV/GP distribution to the columns separately
fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE, shape = 0)
fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)

#Univariate Estimates
fit_mle_1$estimate
fit_mle_2$estimate

#Probability transformation of the component-wise maxima sample
newdata_2 = newdata
newdata_2[,1] = pevd(newdata_2[,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0)
newdata_2[,2] = pevd(newdata_2[,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0)

#Computes the Kendall's tau to obtain the parameter of the Gumbel copula
ken = cor(newdata_2[,1], newdata_2[,2], method = "kendall")
par_clay = 2/((1/ken) - 1)
par_fran = franktau(ken)

#IFM fit
fit.gamma = fitCopula(frankCopula(), newdata_2, method = "ml", start = par_fran)
fit.gamma

#stats
fit = fit.gamma
AIC(fit)
BIC(fit)
confint(fit) 
coef(fit, SE = TRUE)

#Profile likelihood for the Gumbel copula parameter estimate c("Clayton", "Frank")
opt = "Frank"
cop = onacopulaL(opt, list(fit.gamma@estimate,1:dim(seldata)[2]))
n = 2000
U = rnacopula(n,cop)
efm = emle(U, cop)
pfm = profile(efm)
ci = confint(pfm, level=0.95)
logL = function(x) -efm@minuslogl(x)
logL. = Vectorize(logL)
I = c(cop@copula@iTau(0.15), cop@copula@iTau(0.22))
curve(logL., from=I[1], to=I[2], xlab=quote(theta),
      ylab="log-likelihood",
      main=paste("log-likelihood for", opt))
abline(v = ci)

#Goodness of fit
gofCopula(frankCopula(), newdata, N = 200, estim.method = "ml")

#Graphics: Change between frankCopula and claytonCopula
mycopula = frankCopula(param = fit.gamma@estimate, dim = 2)
wireframe2(mycopula, pCopula) 
contour(mycopula, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")
scatterplot3d(newdata_2[,1], newdata_2[,2], pCopula(cbind(newdata_2[,1], newdata_2[,2]), mycopula), color = "red", main = "CDF", xlab = "u1", ylab = "u2", zlab = "C", pch = 21)

#To plot the copula in 3d (pt parameter is activated after the contour section is done)
funt = function(u_1,u_2){pCopula(cbind(u_1, u_2), mycopula)}
pt = FALSE
x_1 = seq(0, 1, by = 0.03)
x_2 = x_1
y = expand.grid(x_1, x_2)
y_1 = cbind(y, funt(y[,1], y[,2]), "grey")
y_2 = cbind(newdata_2, funt(newdata_2[,1], newdata_2[,2]), "red")
if(pt){
  
  y_4 = cbind(newdata_4[,1:2], funt(newdata_4[,1], newdata_4[,2]), "yellow")
  colnames(y_2) = colnames(y_1)
  colnames(y_4) = colnames(y_1)
  df = rbind(y_1, y_2, y_4)
  df = data.frame(x = df[,1], y = df[,2], z = df[,3], colt = df[,4])
  with(df, plot3d(x, y, z, col = colt,  main = "CDF", xlab = "u1", ylab = "u2", zlab = "C", pch = 21, size = 5))
  
}else{
  
  colnames(y_2) = colnames(y_1)
  df = rbind(y_1, y_2)
  df = data.frame(x = df[,1], y = df[,2], z = df[,3], colt = df[,4])
  with(df, plot3d(x, y, z, col = colt,  main = "", xlab = "MRD", ylab = "MRS", zlab = "C", pch = 21, size = 5))
  
}

#Return period
FR = function(u, v, para){
  
  return((-1/para)*(log(1+(((exp(-para*u)-1)*(exp(-para*v) - 1))/(exp(-para) - 1)))))
  
}
"Fcop" = function(u,v, ...) { asCOP(u,v, f = FR, ...) }
T = 200
q = T2prob(T)
T_dual = lmomco::prob2T(duCOP(q,q, cop = Fcop, para = fit.gamma@estimate)) #And case
T_coop = lmomco::prob2T(COP(q,q, cop = Fcop, para = fit.gamma@estimate)) #Or Case
q_dual =  T2prob(T_dual)
q_coop =  T2prob(T_coop)

#To obtain the required q to acquire the exact T period for the OR case
Rt_coop = function(q, alpha, rt){
  
  #return(exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
  #return((max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
  return((-1/alpha)*(log(1+(((exp(-alpha*q)-1)*(exp(-alpha*q) - 1))/(exp(-alpha) - 1)))) - 1 + 1/rt)
  
}
rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.gamma@estimate, rt = T)$root
lmomco::prob2T(COP(rt_coop, rt_coop, cop = Fcop, para = fit.gamma@estimate))

#To obtain the required q to acquire the exact T period for the OR case
Rt_dual = function(q, alpha, rt){
  
  #return(2*q - exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
  #return(2*q - (max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
  return(2*q - (-1/alpha)*(log(1+(((exp(-alpha*q)-1)*(exp(-alpha*q) - 1))/(exp(-alpha) - 1)))) - 1 + 1/rt)
  
}
rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.gamma@estimate, rt = T)$root
lmomco::prob2T(duCOP(rt_dual, rt_dual, cop = Fcop, para = fit.gamma@estimate))

#------------------ Done ~ Plotting Isolines / Design Event -- OR case scenario -----------------

contour(mycopula, pCopula, xlim = c(0, 1), ylim = c(0, 1), main = "Contour plot", nlevels = 15, xlab = "MRD", ylab = "MRS")
points(newdata_2[,1], newdata_2[,2], col = "black")

#Define the Clayton/Frank and the corresponding generator (change the return of the functions)
gumbel = function(u_1, u_2, alpha){
  
  #return(exp(-((-log(u_1))^(alpha) + (-log(u_2))^(alpha))^(1/alpha)))
  #return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
  return((-1/alpha)*(log(1+(((exp(-alpha*u_1)-1)*(exp(-alpha*u_2) - 1))/(exp(-alpha) - 1)))))
  
}
invphigumbel = function (t, alpha = 1){
  
  #exp(-t^(1/alpha))
  #(1+alpha*t)^(-1/alpha)
  (-1/alpha)*(log(1+exp(-t)*(exp(-alpha)-1)))
  
} 
phigumbel = function (t, alpha = 1){
  
  #(-log(t))^alpha
  #(1/alpha)*((t^(-alpha)) - 1)
  -log((exp(-alpha*t)-1)/(exp(-alpha)-1))
  
}

#Obtain elements of the danger zone according to a level z = rt_coop
z = rt_coop
x = seq(0, 1, by = 0.001)
y = rep(0, length(x))
for(i in 1:length(x)){y[i] = invphigumbel(phigumbel(z, alpha = fit.gamma@estimate) - phigumbel(x[i], alpha = fit.gamma@estimate), fit.gamma@estimate)}
x = x[which(y<1&y>0)]
y = y[which(y<1&y>0)]
mat = na.omit(cbind(x, y, gumbel(x,y, alpha = fit.gamma@estimate)))
if(!identical(which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)), integer(0))){
  
  mat = mat[-which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)),]
  
}
mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]

#Graphics of the contour line
plot(newdata_2[,1], newdata_2[,2], xlab = names(newdata_2)[1], ylab = names(newdata_2)[2], main = paste0("MRD", " vs ", "MRS"))
points(mat[,1], mat[,2], type = "l")

#Most-Likely Design
newdata_4 = cbind(mat, rep(0, length(mat[,1])))
for(j in 1:length(mat[,1])){
  
  newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2],0)
  
} 
colnames(newdata_4)[3:4] = c("c", "ML")

pt = TRUE
newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV")
qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV")

plot(newdata_4[,4], type = "l", main = paste("Most-Likely Design Function for Return Period", T, "Years"), ylab = "ML Weight Function")
segments(lty = "dotted", y0 = -10, x1 = which(newdata_4[,4] == max(newdata_4[,4])), x0 = which(newdata_4[,4] == max(newdata_4[,4])), y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])
segments(lty = "dotted", y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4], x0 = which(newdata_4[,4] == max(newdata_4[,4])), x1 = -10, y0 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])

#To calculate the confidence intervals of the RP estimate
CI = function(data, ts, zquant = NA){
  
  newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
  names(newdata) = names(data)[c(2,3)]
  for(i in 1:length(unique(data$hydroYear))){
    
    newdata[i,1] = data[which(data$hydroYear == unique(data$hydroYear)[i]),2][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),2] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),2]))[1]]
    newdata[i,2] = data[which(data$hydroYear == unique(data$hydroYear)[i]),3][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),3] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),3]))[1]]
    
  }
  
  #Fitting the GEV/GP distribution to the columns separately
  fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE, shape = 0)
  fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)
  
  #Probability transformation of the component-wise maxima sample
  newdata_2 = newdata
  newdata_2[,1] = revd(length(newdata_2[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
  newdata_2[,2] = revd(length(newdata_2[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
  newdata_2[,1] = pevd(newdata_2[,1], fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
  newdata_2[,2] = pevd(newdata_2[,2], fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
  
  #IFM fit
  fit.gamma = fitCopula(frankCopula(), newdata_2, method = "ml")
  mycopula = frankCopula(param = fit.gamma@estimate, dim = 2)
  
  #Return period GHcop, GLcop, HR, tEV
  FR = function(u, v, para){
    
    return((-1/para)*(log(1+(((exp(-para*u)-1)*(exp(-para*v) - 1))/(exp(-para) - 1)))))
    
  }
  "Fcop" = function(u,v, ...) { asCOP(u,v, f = FR, ...) }
  T = ts
  q = T2prob(T)
  T_dual = lmomco::prob2T(duCOP(q,q, cop = Fcop, para = fit.gamma@estimate)) #And case
  T_coop = lmomco::prob2T(COP(q,q, cop = Fcop, para = fit.gamma@estimate)) #Or Case
  q_dual =  T2prob(T_dual)
  q_coop =  T2prob(T_coop)
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_coop = function(q, alpha, rt){
    
    #return(exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
    #return((max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    return((-1/alpha)*(log(1+(((exp(-alpha*q)-1)*(exp(-alpha*q) - 1))/(exp(-alpha) - 1)))) - 1 + 1/rt)
    
  }
  rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.gamma@estimate, rt = T)$root
  lmomco::prob2T(COP(rt_coop, rt_coop, cop = Fcop, para = fit.gamma@estimate))
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_dual = function(q, alpha, rt){
    
    #return(2*q - exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
    #return(2*q - (max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    return(2*q - (-1/alpha)*(log(1+(((exp(-alpha*q)-1)*(exp(-alpha*q) - 1))/(exp(-alpha) - 1)))) - 1 + 1/rt)
    
  }
  rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.gamma@estimate, rt = T)$root
  lmomco::prob2T(duCOP(rt_dual, rt_dual, cop = Fcop, para = fit.gamma@estimate))
  
  #Define the Clayton/Frank and the corresponding generator (change the return of the functions)
  gumbel = function(u_1, u_2, alpha){
    
    #return(exp(-((-log(u_1))^(alpha) + (-log(u_2))^(alpha))^(1/alpha)))
    #return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
    return((-1/alpha)*(log(1+(((exp(-alpha*q)-1)*(exp(-alpha*q) - 1))/(exp(-alpha) - 1)))))
    
  }
  invphigumbel = function (t, alpha = 1){
    
    #exp(-t^(1/alpha))
    #(1+alpha*t)^(-1/alpha)
    (-1/alpha)*(log(1+exp(-t)*(exp(-alpha)-1)))
    
  } 
  phigumbel = function (t, alpha = 1){
    
    #(-log(t))^alpha
    #(1/alpha)*((t^(-alpha)) - 1)
    -log((exp(-alpha*t)-1)/(exp(-alpha)-1))
    
  }
  
  #Obtain elements of the danger zone according to a level z = 0.08/0.23
  if(is.na(zquant) == TRUE){
    
    z = rt_coop
    
  }else{
    
    z = zquant
    
  }
  x = seq(0, 1, by = 0.001)
  y = rep(0, length(x))
  for(i in 1:length(x)){y[i] = invphigumbel(phigumbel(z, alpha = fit.gamma@estimate) - phigumbel(x[i], alpha = fit.gamma@estimate), fit.gamma@estimate)}
  x = x[which(y<1&y>0)]
  y = y[which(y<1&y>0)]
  mat = na.omit(cbind(x, y, gumbel(x,y, alpha = fit.gamma@estimate)))
  if(!identical(which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)), integer(0))){
    
    mat = mat[-which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)),]
    
  }
  mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]
  
  #Most-Likely Design
  newdata_4 = cbind(mat, rep(0, length(mat[,1])))
  for(j in 1:length(mat[,1])){
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2],0)
    
  } 
  colnames(newdata_4)[3:4] = c("c", "ML")
  
  pt = TRUE
  newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
  dd=qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV")
  rttt=qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV")
  
  return(c(dd,rttt))
  
}
rep = 1000
zquantt = 0.974
mat_CI = matrix(nrow = rep, ncol = 2)
for(jj in 1:rep){
  
  tryCatch({
    
    print(jj)
    
    boot_CI = CI(seldata, ts = T, zquant = zquantt)
    mat_CI[jj,1] = as.numeric(boot_CI[1])
    mat_CI[jj,2] = as.numeric(boot_CI[2])
    
  }, error=function(e){})
  
}
range(na.omit(mat_CI[,1]))
range(na.omit(mat_CI[,2]))

#------------------ Done ~ Plotting Isolines / Design Event -- AND case scenario -----------------

cc = seq(0, 1, by = 0.001)
mat1 = matrix(0, nrow = length(cc), ncol = length(cc), byrow=TRUE)
for(o in 1:length(cc)){
  
  for(r in 1:length(cc)){
    
    mat1[o,r] = cc[o] + cc[r] - gumbel(cc[o], cc[r], alpha = fit.gamma@estimate)
    
  }
  
}
contour(cc,cc,mat1, nlevels = 15, lty = 1, xlim = c(0, 1), ylim = c(0, 1), main = "Dual Contour plot", xlab = "MRD", ylab = "MRS")
points(newdata_2[,1], newdata_2[,2], col = "black")

#AND case scenario
#Obtain elements of the danger zone according to a level 
z = rt_dual
x = seq(0, 1, by = 0.001)
y = rep(0, length(x))
dual = function(u, v, s, alpha){
  
  - u - v + gumbel(u, v, alpha) + s
  
}
for(t in 1:length(y)){
  
  tryCatch({
    
    y[t] = uniroot(dual, interval = c(0, 1), v = x[t] , s = z, alpha = fit.gamma@estimate)$root
    
  }, error=function(e){})
  
}

mat = na.omit(cbind(x, y, x + y - gumbel(x, y, alpha = fit.gamma@estimate)))

#Graphics of the contour line
plot(newdata_2[,1], newdata_2[,2], xlab = names(newdata_2)[1], ylab = names(newdata_2)[2], main = paste0(names(newdata_2)[1], " vs ", names(newdata_2)[2]))
points(mat[,1], mat[,2], type = "l")

#Most-Likely Design
mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]
newdata_4 = cbind(mat, rep(0, length(mat[,1])))
for(j in 1:length(mat[,1])){
  
  newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], 0)
  
} 
colnames(newdata_4)[3:4] = c("c", "ML")

pt = TRUE
newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV")
qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV")

plot(newdata_4[,4], type = "l", main = "Most-Likely Design")
segments(lty = "dotted", y0 = -10, x1 = which(newdata_4[,4] == max(newdata_4[,4])), x0 = which(newdata_4[,4] == max(newdata_4[,4])), y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])
segments(lty = "dotted", y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4], x0 = which(newdata_4[,4] == max(newdata_4[,4])), x1 = -10, y0 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])

#To calculate the confidence intervals of the RP estimate
CI = function(data, ts, zquant = NA){
  
  newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
  names(newdata) = names(data)[c(2,3)]
  for(i in 1:length(unique(data$hydroYear))){
    
    newdata[i,1] = data[which(data$hydroYear == unique(data$hydroYear)[i]),2][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),2] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),2]))[1]]
    newdata[i,2] = data[which(data$hydroYear == unique(data$hydroYear)[i]),3][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),3] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),3]))[1]]
    
  }
  
  #Fitting the GEV distribution to the columns separately
  fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE)
  fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE)
  
  #Probability transformation of the component-wise maxima sample
  newdata_2 = newdata
  newdata_2[,1] = revd(length(newdata_2[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3])
  newdata_2[,2] = revd(length(newdata_2[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3])
  newdata_2[,1] = pevd(newdata_2[,1], fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3])
  newdata_2[,2] = pevd(newdata_2[,2], fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3])
  
  #IFM fit
  fit.beta = fitCopula(frankCopula(), newdata_2, method = "ml")
  mycopula = frankCopula(param = fit.gamma@estimate, dim = 2)
  
  #Return period GHcop, GLcop, HR, tEV
  #q = 1-1/500
  z = rt_dual
  q = T2prob(ts)
  T_dual = lmomco::prob2T(duCOP(q,q, cop = GHcop, para = fit.beta@estimate)) #And case
  T_coop = lmomco::prob2T(COP(q,q, cop = GHcop, para = fit.beta@estimate)) #Or Case
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_coop = function(q, alpha, rt){
    
    #return(invphigumbel(phigumbel(u_1,alpha) + phigumbel(u_2,alpha), alpha))
    return(exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
    
  }
  rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.beta@estimate, rt = ts)$root
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_dual = function(q, alpha, rt){
    
    #return(invphigumbel(phigumbel(u_1,alpha) + phigumbel(u_2,alpha), alpha))
    return(2*q - exp(-((-log(q))^(alpha) + (-log(q))^(alpha))^(1/alpha)) - 1 + 1/rt)
    
  }
  rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.beta@estimate, rt = ts)$root
  
  #Define the Gumbel and the corresponding generator
  gumbel = function(u_1, u_2, alpha){
    
    #return(exp(-((-log(u_1))^(alpha) + (-log(u_2))^(alpha))^(1/alpha)))
    #return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
    return((-1/alpha)*(log(1+(((exp(-alpha*u_1)-1)*(exp(-alpha*u_2) - 1))/(exp(-alpha) - 1)))))
    
  }
  invphigumbel = function (t, alpha = 1){
    
    #exp(-t^(1/alpha))
    #(1+alpha*t)^(-1/alpha)
    (-1/alpha)*(log(1+exp(-t)*(exp(-alpha)-1)))
    
  } 
  phigumbel = function (t, alpha = 1){
    
    #(-log(t))^alpha
    #(1/alpha)*((t^(-alpha)) - 1)
    -log((exp(-alpha*t)-1)/(exp(-alpha)-1))
    
  }
  
  #Obtain elements of the danger zone according to a level z = 0.08/0.23
  z = rt_coop
  x = seq(0, 1, by = 0.001)
  y = rep(0, length(x))
  for(i in 1:length(x)){y[i] = invphigumbel(phigumbel(z, alpha = fit.beta@estimate) - phigumbel(x[i], alpha = fit.beta@estimate), fit.beta@estimate)}
  
  mat = na.omit(cbind(x, y, gumbel(x,y, alpha = fit.beta@estimate)))
  if(!identical(which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)), integer(0))){
    
    mat = mat[-which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)),]
    
  }
  
  #Most-Likely Design
  mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]
  newdata_4 = cbind(mat, rep(0, length(mat[,1])))
  for(j in 1:length(mat[,1])){
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],fit_mle_1$estimate[3])*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2], fit_mle_1$estimate[3])
    
  } 
  colnames(newdata_4)[3:4] = c("c", "ML")
  
  pt = TRUE
  newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
  dd=qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV")
  rttt=qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV")
  
  return(c(dd,rttt))
  
}
rep = 1000
zquantt = 0.974
mat_CI = matrix(nrow = rep, ncol = 2)
for(jj in 1:rep){
  
  tryCatch({
    
    print(jj)
    
    boot_CI = CI(seldata, ts = T)
    mat_CI[jj,1] = as.numeric(boot_CI[1])
    mat_CI[jj,2] = as.numeric(boot_CI[2])
    
  }, error=function(e){})
  
}
range(na.omit(mat_CI[,1]))
range(na.omit(mat_CI[,2]))
