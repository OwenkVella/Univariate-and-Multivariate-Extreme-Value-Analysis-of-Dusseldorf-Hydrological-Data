#install.packages("evd")
#install.packages("latex2exp")
library("evd")
library("latex2exp")

#To remove variables/ functions and plots
rm(list = ls())
graphics.off()

####-GEV-####

#Plotting margins
layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), heights = c(8, 1))
par(mai=rep(0.5, 4))

#Generating a sequence
x = seq(-5, 5, by = 0.1)
y = seq(-1, 8, by = 0.1)

#Extreme Value Densities
Gumbel_density  = evd::dgev(x, loc = 0, scale = 1)
Frechet_density = evd::dgev(x, shape = 0.5, loc = 0, scale = 1)
Weibull_density = evd::dgev(x, shape = -0.5, loc = 0, scale = 1)

#Extreme Value Distributions
Gumbel_dist  = evd::pgev(x, loc = 0, scale = 1)
Frechet_dist = evd::pgev(x, shape = 0.5, loc = 0, scale = 1)
Weibull_dist = evd::pgev(x, shape = -0.5, loc = 0, scale = 1)

#Extreme Value densities plots
plot(c(x,x,x), c(Gumbel_density,Frechet_density, Weibull_density), type = 'n', xlab = "", ylab = "", las = 1)
lines(x, Gumbel_density,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(x, Frechet_density, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(x, Weibull_density, type = 'l', lty = 3, col = 'blue', lwd = 2)

#Extreme Value distribution plots
plot(c(x,x,x), c(Gumbel_dist,Frechet_dist, Weibull_dist), type = 'n', xlab = "", ylab = "", las = 1)
lines(x, Gumbel_dist,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(x, Frechet_dist, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(x, Weibull_dist, type = 'l', lty = 3, col = 'blue', lwd = 2)

#Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x = "center", ncol = 3.2, legend = c(TeX(r'($\xi = 0$)'),TeX(r'($\xi = 0.5$)'),TeX(r'($\xi = -0.5$)')), fill = c("green","red","blue"), pt.cex = 1, cex = 1.2)

####-GGEV-####

#Plotting margins
layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), heights = c(8, 1))
par(mai=rep(0.5, 4))

#Selecting shape parameter/cut-off
alpha = 0.5
k = 3

#Generalized Extreme Order Statistics Densities
Gumbel_density_k = exp(-exp(-x))*(exp(-k*x))/factorial(k-1)
Frechet_density_k = ifelse(1 + x*alpha > 0, (exp(-(1 + x*alpha)^(-1/alpha))*(1 + x*alpha)^(-((1/alpha)*k+1)))/factorial(k-1),0)
Weibull_density_k = ifelse(1 - x*alpha > 0, (exp(-(1 - x*alpha)^(1/alpha))*(1 - x*alpha)^(-((-1/alpha)*k+1)))/factorial(k-1),0)

#Extreme Value Distributions (GEV)
Gumbel_dist  = exp(-exp(-x))
Frechet_dist = ifelse(1 + x*alpha > 0, exp(-(1 + x*alpha)^(-1/alpha)),0)
Weibull_dist = ifelse(1 - x*alpha > 0, exp(-(1 - x*alpha)^(1/alpha)),1)

#Obtaining the generalized order statistics distribution part
Gumb_r = c()
Frec_r = c()
Weib_r = c()
tvec_r = matrix(0, nrow = length(x), ncol = 3, byrow = TRUE)
colnames(tvec_r) = c("Gumbel_k", "Frechet_k", "Weibull_k")
for(f in 1:length(x)){
  
  for(j in 0:(k-1)) {
    
    term_Gumb_r = (exp(-j*x[f]))/(factorial(j))
    Gumb_r = c(Gumb_r,term_Gumb_r)
    
    term_Frec_r = ifelse(1 + x[f]*alpha > 0,((1 + x[f]*alpha)^(-(j/alpha))),0)/(factorial(j))
    Frec_r = c(Frec_r,term_Frec_r)
    
    term_Weib_r = ifelse(1 - x[f]*alpha > 0,((1 - x[f]*alpha)^(j/alpha))/(factorial(j)),1)
    Weib_r = c(Weib_r,term_Weib_r)
    
  }
  
  tvec_r[f,] = c(sum(Gumb_r),sum(Frec_r),ifelse(x[f]==0,sum(Frec_r),ifelse(sum(Weib_r) == k,1,sum(Weib_r))))
  Gumb_r = c()
  Frec_r = c()
  Weib_r = c()
  
}

#Generalized Extreme Order Statistics Distributions
Gumbel_dist_k  = Gumbel_dist*tvec_r[,1]
Frechet_dist_k = Frechet_dist*tvec_r[,2]
Weibull_dist_k = Weibull_dist*tvec_r[,3]

#Extreme Value densities plots
plot(c(x,x,x), c(Gumbel_density_k,Frechet_density_k, Weibull_density_k), type = 'n', xlab = "", ylab = "", las = 1)
lines(x, Gumbel_density_k,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(x, Frechet_density_k, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(x, Weibull_density_k, type = 'l', lty = 3, col = 'blue', lwd = 2)

#Extreme Value distribution plots
plot(c(x,x,x), c(Gumbel_dist_k,Frechet_dist_k, Weibull_dist_k), type = 'n', xlab = "", ylab = "", las = 1)
lines(x, Gumbel_dist_k,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(x, Frechet_dist_k, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(x, Weibull_dist_k, type = 'l', lty = 3, col = 'blue', lwd = 2)

#Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x = "center", ncol = 3.2, legend = c(TeX(r'($\xi = 0$)'),TeX(r'($\xi = 0.5$)'),TeX(r'($\xi = -0.5$)')), fill = c("green","red","blue"), pt.cex = 1, cex = 1.2)

####-GPD-####

#Plotting margins
layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), heights = c(8, 1))
par(mai = rep(0.5, 4))

#GPD densities
Expo_density  = evd::dgpd(y, loc = 0, scale = 1)
Pareto_density = evd::dgpd(y, shape = 0.5, loc = 0, scale = 1)
Pareto_type_II_density = evd::dgpd(y, shape = -0.5, loc = 0, scale = 1)

#GPD distributions
Expo_dist           = evd::pgpd(y, loc = 0, scale = 1)
Pareto_dist         = evd::pgpd(y, shape = 0.5, loc = 0, scale = 1)
Pareto_type_II_dist = evd::pgpd(y, shape = -0.5, loc = 0, scale = 1)

#GPD densities plot
plot(c(y,y,y), c(Expo_density,Pareto_density, Pareto_type_II_density), type = 'n', xlab = "", ylab = "", las = 1)
lines(y, Expo_density,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(y, Pareto_density, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(y, Pareto_type_II_density, type = 'l', lty = 3, col = 'blue', lwd = 2)

#GPD densities plot
plot(c(y,y,y), c(Expo_dist, Pareto_dist, Pareto_type_II_dist), type = 'n', xlab = "", ylab = "", las = 1)
lines(y, Expo_dist,  type = 'l', lty = 1, col = 'green', lwd = 2)
lines(y, Pareto_dist, type = 'l', lty = 2, col = 'red', lwd = 2)
lines(y, Pareto_type_II_dist, type = 'l', lty = 3, col = 'blue', lwd = 2)

#Legend
par(mai = c(0,0,0,0))
plot.new()
legend(x = "center", ncol = 3.2, legend = c(TeX(r'($\xi = 0$)'),TeX(r'($\xi = 0.5$)'),TeX(r'($\xi = -0.5$)')), fill = c("green","red","blue"), pt.cex = 1, cex = 1.2)
