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
seldata = P1_df[,c(2,3,4,1)]
point = "P1"

#------------------ Done ~ Fitting Distribution  -----------------

summary(seldata)

#To obtain the annual maxima or the excess (need to amend accordingly)
newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
names(newdata) = names(seldata)[c(2,3)]
for(i in 1:length(unique(seldata$hydroYear))){
  
  newdata[i,1] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2]))]
  newdata[i,2] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3]))]
  
}

#Showing component-wise maxima from disjoint blocks in red
#plot(new_data[,c(7,17)], col = "black")
plot(seldata[, c(2,3)], col = "black", xlim = 0.9*c(min(newdata[,1], seldata[,2]), 1.1*max(newdata[,1], seldata[,2])), ylim = 0.9*c(min(newdata[,2], seldata[,3]), 1.1*max(newdata[,2], seldata[,3])))
points(newdata, col = "red")

hk = rbind(cbind(seldata[, c(2,3)],"G" = 1), cbind(newdata, "G" = 2), cbind(intersect(newdata, seldata[,c(2,3)]),"G" = 3))
hk$G = as.factor(hk$G)
p1 = ggplot(hk, aes(x = MRS, y = MRD, colour = G)) + geom_point() + theme_bw()+ theme(legend.position = "none") + scale_color_manual(values = c("black","red", "green"))
ggMarginal(p1, type="boxplot", fill = "#0099f8") 

#To plot a matrix of scatter plots
upper.panel = function(x, y){
  points(x,y, pch = 19)
}
panel.hist = function(x, ...){
  usr = par("usr"); on.exit(par("usr"))
  par(usr = c(usr[1:2], 0, 1.5) )
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$counts; y = y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor = function(x, y){
  usr = par("usr"); on.exit(par("usr"))
  par(usr = c(0, 1, 0, 1))
  r = round(kendall.tau(x, y), digits=2)
  t = round(cor(x, y, method = "pearson"), digits=2)
  txt = bquote(tau == .(round(r,2)))
  txt_1 = bquote(rho == .(round(t,2)))
  cex.cor = 0.8/strwidth(txt)
  text(0.5, 0.6, txt_1, cex = 2.5)
  text(0.5, 0.4, txt, cex = 2.5)
}
pairs(newdata, lower.panel = panel.cor,  upper.panel = upper.panel,diag.panel=panel.hist)

#Fitting the GEV distribution to the columns separately
if(point == "P1"){
  
  fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE)
  fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE)
  
}else{
  
  fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE, shape = 0)
  fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)
  
}

#Univariate Estimates
fit_mle_1$estimate
fit_mle_2$estimate

#Probability transformation of the component-wise maxima sample
newdata_2 = newdata
if(point == "P1"){
  
  newdata_2[,1] = pevd(newdata_2[,1], fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3])
  newdata_2[,2] = pevd(newdata_2[,2], fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3]) 
  
}else{
  
  newdata_2[,1] = pevd(newdata_2[,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0)
  newdata_2[,2] = pevd(newdata_2[,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0)
  
}

#Computes the Kendall's tau to obtain the parameter of the Clayton copula
ken = cor(newdata_2[,1], newdata_2[,2], method = "kendall")
par_clay = 2/((1/ken) - 1)

#IFM fit
fit.beta = fitCopula(claytonCopula(), newdata_2, method = "ml")
fit.beta

#Stats
AIC(fit.beta)
BIC(fit.beta)
confint(fit.beta) 
coef(fit.beta, SE = TRUE)

#Profile likelihood for the Clayton copula parameter estimate
opt = "Clayton"
cop = onacopulaL(opt, list(fit.beta@estimate,1:dim(seldata)[2]))
n = 2000
U = rnacopula(n,cop)
efm = emle(U, cop)
pfm = profile(efm)
ci = confint(pfm, level=0.95)
logL = function(x) -efm@minuslogl(x)
logL. = Vectorize(logL)
I = c(cop@copula@iTau(0.22), cop@copula@iTau(0.3))
curve(logL., from=I[1], to=I[2], xlab=quote(theta),
      ylab="log-likelihood",
      main=paste("log-likelihood for", opt))
abline(v = ci)

#Goodness of fit
gofCopula(claytonCopula(), newdata, N = 200, estim.method = "ml")

#Graphics: Change between frankCopula and claytonCopula
mycopula = claytonCopula(param = fit.beta@estimate, dim = 2)
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
T = 20
q = T2prob(T)
T_dual = lmomco::prob2T(duCOP(q,q, cop = CLcop, para = fit.beta@estimate)) #And case
T_coop = lmomco::prob2T(COP(q,q, cop = CLcop, para = fit.beta@estimate)) #Or Case
q_dual =  T2prob(T_dual)
q_coop =  T2prob(T_coop)

#To obtain the required q to acquire the exact T period for the OR case
Rt_coop = function(q, alpha, rt){
  
  return((max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
  
}
rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.beta@estimate, rt = T)$root
lmomco::prob2T(COP(rt_coop, rt_coop, cop = CLcop, para = fit.beta@estimate))

#To obtain the required q to acquire the exact T period for the AND case
Rt_dual = function(q, alpha, rt){
  
  return(2*q - (max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
  
}
rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.beta@estimate, rt = T)$root
lmomco::prob2T(duCOP(rt_dual, rt_dual, cop = CLcop, para = fit.beta@estimate))

#------------------ Done ~ Plotting Isolines / Design Event -- OR case scenario -----------------

contour(mycopula, pCopula, xlim = c(0, 1), ylim = c(0, 1), main = "Contour plot", nlevels = 15, xlab = "MRD", ylab = "MRS")
points(newdata_2[,1], newdata_2[,2], col = "black")

#Define the Clayton and the corresponding generator (change the return of the functions)
clayton = function(u_1, u_2, alpha){
  
  return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
  
}
invphiclayton = function (t, alpha = 1){
  
  (1+alpha*t)^(-1/alpha)
  
} 
phiclayton = function (t, alpha = 1){
  
  (1/alpha)*((t^(-alpha)) - 1)
  
}

#Obtain elements of the danger zone according to a level z = rt_coop
z = rt_coop
x = seq(0, 1, by = 0.001)
y = rep(0, length(x))
for(i in 1:length(x)){y[i] = invphiclayton(phiclayton(z, alpha = fit.beta@estimate) - phiclayton(x[i], alpha = fit.beta@estimate), fit.beta@estimate)}
x = x[which(y<1&y>0)]
y = y[which(y<1&y>0)]
mat = na.omit(cbind(x, y, clayton(x,y, alpha = fit.beta@estimate)))
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
  
  if(point == "P1"){
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],fit_mle_1$estimate[3])*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], fit_mle_2$estimate[3])
    
  }else{
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2],0)
    
  }
  
} 
colnames(newdata_4)[3:4] = c("c", "ML")

pt = TRUE
newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
if(point == "P1"){
  
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"))
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"))
  
}else{
  
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"))
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"))
  
}

plot(newdata_4[,4], type = "l", main = paste("Most-Likely Design Function for Return Period", T, "Years"), ylab = "ML Weight Function")
segments(lty = "dotted", y0 = -10, x1 = which(newdata_4[,4] == max(newdata_4[,4])), x0 = which(newdata_4[,4] == max(newdata_4[,4])), y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])
segments(lty = "dotted", y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4], x0 = which(newdata_4[,4] == max(newdata_4[,4])), x1 = -10, y0 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])

#To calculate the confidence intervals of the RP estimate
CI = function(data, ts){
  
  newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
  names(newdata) = names(data)[c(2,3)]
  for(i in 1:length(unique(data$hydroYear))){
    
    newdata[i,1] = data[which(data$hydroYear == unique(data$hydroYear)[i]),2][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),2] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),2]))[1]]
    newdata[i,2] = data[which(data$hydroYear == unique(data$hydroYear)[i]),3][which(data[which(data$hydroYear == unique(data$hydroYear)[i]),3] == max(data[which(data$hydroYear == unique(data$hydroYear)[i]),3]))[1]]
    
  }
  
  #Fitting the GEV distribution to the columns separately
  if(point == "P1"){
    
    fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE)
    fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE)
    
  }else{
    
    fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE, shape = 0)
    fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)
    
  }
  
  #Probability transformation of the component-wise maxima sample
  newdata_2 = matrix(nrow = dim(newdata)[1], ncol = 2)
  if(point == "P1"){
    
    newdata_2[,1] = pevd(revd(length(newdata[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3])
    newdata_2[,2] = pevd(revd(length(newdata[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3])
    
  }else{
    
    newdata_2[,1] = revd(length(newdata_2[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
    newdata_2[,2] = revd(length(newdata_2[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
    newdata_2[,1] = pevd(newdata_2[,1], fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
    newdata_2[,2] = pevd(newdata_2[,2], fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
    
  }
  
  #IFM fit
  fit.beta = fitCopula(claytonCopula(), newdata_2, method = "ml")
  mycopula = claytonCopula(param = fit.beta@estimate, dim = 2)
  
  #Return period 
  #q = 1-1/500
  q = T2prob(ts)
  T_dual = lmomco::prob2T(duCOP(q,q, cop = CLcop, para = fit.beta@estimate)) #AND case
  T_coop = lmomco::prob2T(COP(q,q, cop = CLcop, para = fit.beta@estimate)) #OR Case
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_coop = function(q, alpha, rt){
    
    return((max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    
  }
  rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.beta@estimate, rt = ts)$root
  
  #To obtain the required q to acquire the exact T period for the AND case
  Rt_dual = function(q, alpha, rt){
    
    return(2*q - (max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    
  }
  rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.beta@estimate, rt = ts)$root
  
  #Define the Clayton and the corresponding generator (change the return of the functions)
  clayton = function(u_1, u_2, alpha){
    
    return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
    
  }
  invphiclayton = function (t, alpha = 1){
    
    (1+alpha*t)^(-1/alpha)
    
  } 
  phiclayton = function (t, alpha = 1){
    
    (1/alpha)*((t^(-alpha)) - 1)
    
  }
  
  #Obtain elements of the danger zone according to a level z = rt_coop
  z = rt_coop
  x = seq(0, 1, by = 0.00001)
  y = rep(0, length(x))
  for(i in 1:length(x)){y[i] = invphiclayton(phiclayton(z, alpha = fit.beta@estimate) - phiclayton(x[i], alpha = fit.beta@estimate), fit.beta@estimate)}
  x = x[which(y<1&y>0)]
  y = y[which(y<1&y>0)]
  mat = na.omit(cbind(x, y, clayton(x,y, alpha = fit.beta@estimate)))
  if(!identical(which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)), integer(0))){
    
    mat = mat[-which((mat[,1] > 1 | mat[,1] < 0) | (mat[,2] > 1 | mat[,2] < 0)),]
    
  }
  mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]
  
  #Most-Likely Design
  newdata_4 = cbind(mat, rep(0, length(mat[,1])))
  for(j in 1:length(mat[,1])){
    
    if(point == "P1"){
      
      newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],fit_mle_1$estimate[3])*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], fit_mle_2$estimate[3])
      
    }else{
      
      newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2],0)
      
    }
    
  } 
  colnames(newdata_4)[3:4] = c("c", "ML")
  
  newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
  if(point == "P1"){
    
    dd = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV")
    rttt = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV")
    
  }else{
    
    dd = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV")
    rttt = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV")
    
  }
  
  return(c(dd,rttt))
  
}
rep = 1000
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

#------------------ Done ~ Plotting Isolines / Design Event -- AND case scenario -----------------

cc = seq(0, 1, by = 0.001)
mat1 = matrix(0, nrow = length(cc), ncol = length(cc), byrow=TRUE)
for(o in 1:length(cc)){
  
  for(r in 1:length(cc)){
    
    mat1[o,r] = cc[o] + cc[r] - clayton(cc[o], cc[r], alpha = fit.beta@estimate)
    
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
  
  - u - v + clayton(u, v, alpha) + s
  
}
for(t in 1:length(y)){
  
  tryCatch({
    
    y[t] = uniroot(dual, interval = c(0, 1), v = x[t] , s = z, alpha = fit.beta@estimate)$root
    
  }, error=function(e){})
  
}
x = x[which(y<1&y>0)]
y = y[which(y<1&y>0)]
mat = na.omit(cbind(x, y, x + y - clayton(x, y, alpha = fit.beta@estimate)))
mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]

#Graphics of the contour line
plot(newdata_2[,1], newdata_2[,2], xlab = names(newdata_2)[1], ylab = names(newdata_2)[2], main = paste0("MRD", " vs ", "MRS"))
points(mat[,1], mat[,2], type = "l")

#Most-Likely Design
newdata_4 = cbind(mat, rep(0, length(mat[,1])))
for(j in 1:length(mat[,1])){
  
  if(point == "P1"){
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],fit_mle_1$estimate[3])*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], fit_mle_2$estimate[3])
    
  }else{
    
    newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], 0)
    
  }
  
} 
colnames(newdata_4)[3:4] = c("c", "ML")

pt = TRUE
newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
if(point == "P1"){
  
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"))
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"))
  
}else{
  
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"))
  print(qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"))
  
}

plot(newdata_4[,4], type = "l", main = "Most-Likely Design")
segments(lty = "dotted", y0 = -10, x1 = which(newdata_4[,4] == max(newdata_4[,4])), x0 = which(newdata_4[,4] == max(newdata_4[,4])), y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])
segments(lty = "dotted", y1 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4], x0 = which(newdata_4[,4] == max(newdata_4[,4])), x1 = -10, y0 = newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),4])

#To calculate the confidence intervals of the RP estimate
CI = function(data, ts){
  
  #To obtain the annual maxima or the excess (need to amend accordingly)
  newdata = data.frame("V1" = as.numeric(), "V2" = as.numeric())
  names(newdata) = names(seldata)[c(2,3)]
  for(i in 1:length(unique(seldata$hydroYear))){
    
    newdata[i,1] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),2]))]
    newdata[i,2] = seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3][which(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3] == max(seldata[which(seldata$hydroYear == unique(seldata$hydroYear)[i]),3]))]
    
  }
  
  #Fitting the GEV distribution to the columns separately
  if(point == "P1"){
    
    fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE)
    fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE)
    
  }else{
    
    fit_mle_1 = evd::fgev(x = newdata[,1], std.err = TRUE, corr = TRUE, shape = 0)
    fit_mle_2 = evd::fgev(x = newdata[,2], std.err = TRUE, corr = TRUE, shape = 0)
    
  }
  
  #Probability transformation of the component-wise maxima sample
  newdata_2 = matrix(nrow = dim(newdata)[1], ncol = 2)
  if(point == "P1"){
    
    newdata_2[,1] = pevd(revd(length(newdata[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],fit_mle_1$estimate[3])
    newdata_2[,2] = pevd(revd(length(newdata[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],fit_mle_2$estimate[3])
    
  }else{
    
    newdata_2[,1] = revd(length(newdata_2[,1]), fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
    newdata_2[,2] = revd(length(newdata_2[,2]), fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
    newdata_2[,1] = pevd(newdata_2[,1], fit_mle_1$estimate[1],fit_mle_1$estimate[2],0)
    newdata_2[,2] = pevd(newdata_2[,2], fit_mle_2$estimate[1],fit_mle_2$estimate[2],0)
    
  }
  
  #IFM fit
  fit.beta = fitCopula(claytonCopula(), newdata_2, method = "ml")
  fit.beta
  
  #Return period
  T = ts
  q = T2prob(T)
  T_dual = lmomco::prob2T(duCOP(q,q, cop = CLcop, para = fit.beta@estimate)) #AND case
  T_coop = lmomco::prob2T(COP(q,q, cop = CLcop, para = fit.beta@estimate)) #OR Case
  q_dual = T2prob(T_dual)
  q_coop = T2prob(T_coop)
  
  #To obtain the required q to acquire the exact T period for the OR case
  Rt_coop = function(q, alpha, rt){
    
    return((max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    
  }
  rt_coop = uniroot(Rt_coop, interval = c(0,1), alpha = fit.beta@estimate, rt = T)$root
  lmomco::prob2T(COP(rt_coop, rt_coop, cop = CLcop, para = fit.beta@estimate))
  
  #To obtain the required q to acquire the exact T period for the AND case
  Rt_dual = function(q, alpha, rt){
    
    return(2*q - (max((q^(-alpha))+(q^(-alpha))-1,0))^(-1/alpha) - 1 + 1/rt)
    
  }
  rt_dual = uniroot(Rt_dual, interval = c(0,1), alpha = fit.beta@estimate, rt = T)$root
  lmomco::prob2T(duCOP(rt_dual, rt_dual, cop = CLcop, para = fit.beta@estimate))
  
  #Define the Clayton and the corresponding generator (change the return of the functions)
  clayton = function(u_1, u_2, alpha){
    
    return((max((u_1^(-alpha))+(u_2^(-alpha))-1,0))^(-1/alpha))
    
  }
  invphiclayton = function (t, alpha = 1){
    
    (1+alpha*t)^(-1/alpha)
    
  } 
  phiclayton = function (t, alpha = 1){
    
    (1/alpha)*((t^(-alpha)) - 1)
    
  }
  
  #AND case scenario: Obtain elements of the danger zone according to a level 
  z = rt_dual
  x = seq(0, 1, by = 0.001)
  y = rep(0, length(x))
  dual = function(u, v, s, alpha){
    
    - u - v + clayton(u, v, alpha) + s
    
  }
  for(t in 1:length(y)){
    
    tryCatch({
      
      y[t] = uniroot(dual, interval = c(0, 1), v = x[t] , s = z, alpha = fit.beta@estimate)$root
      
    }, error=function(e){})
    
  }
  x = x[which(y<1&y>0)]
  y = y[which(y<1&y>0)]
  mat = na.omit(cbind(x, y, x + y - clayton(x, y, alpha = fit.beta@estimate)))
  mat = mat[-which(mat[,1] == 0 | mat[,1] == 1 | mat[,2] == 0 | mat[,2] == 1),]
  
  #Most-Likely Design
  newdata_4 = cbind(mat, rep(0, length(mat[,1])))
  for(j in 1:length(mat[,1])){
    
    if(point == "P1"){
      
      newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],fit_mle_1$estimate[3])*devd(qevd(as.data.frame(mat)[j,2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2], fit_mle_2$estimate[3])
      
    }else{
      
      newdata_4[j,4] = dCopula(cbind(mat[j,1], mat[j,2]), mycopula)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV"), fit_mle_1$estimate[1], fit_mle_1$estimate[2],0)*devd(qevd(as.data.frame(mat)[j,1], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV"), fit_mle_2$estimate[1], fit_mle_2$estimate[2],0)
      
    }
    
  } 
  colnames(newdata_4)[3:4] = c("c", "ML")
  
  newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),]
  if(point == "P1"){
    
    dd = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = fit_mle_1$estimate[3], type = "GEV")
    rttt = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = fit_mle_2$estimate[3], type = "GEV")
    
  }else{
    
    dd = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),1], loc = fit_mle_1$estimate[1], scale = fit_mle_1$estimate[2], shape = 0, type = "GEV")
    rttt = qevd(newdata_4[which(newdata_4[,4] == max(newdata_4[,4])),2], loc = fit_mle_2$estimate[1], scale = fit_mle_2$estimate[2], shape = 0, type = "GEV")
    
  }
  
  return(c(dd,rttt))
  
}
rep = 1000
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
