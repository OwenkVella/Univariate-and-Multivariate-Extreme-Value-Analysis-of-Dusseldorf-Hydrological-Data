#Required Packages
#install.packages(c("lmtest", "tcltk", "forecast", "dplyr", "rdwd", "zoo", "pastecs", "ggplot2", "VGAM", "utils", "pracma", "RcppRoll", "twosamples", "extRemes", "copula", "tidyr", "lubridate", "Kendall", "tseries"))
lapply(c("lmtest", "tcltk", "forecast", "dplyr", "rdwd", "zoo", "pastecs", "ggplot2", "VGAM", "utils", "pracma", "RcppRoll", "twosamples", "extRemes", "copula", "tidyr", "lubridate", "Kendall", "tseries"), require, character.only = TRUE)

#To clean
rm(list = ls())
graphics.off()

#Set the directory
setwd(tk_choose.dir(default = "", caption = "Select directory"))

#Load the data from the ftp server
#First step: Station data
statID = 1078
link_2 = selectDWD(id = statID, res = "daily", var = "more_precip", per = "historical")
file_2 = dataDWD(link_2, read = FALSE)
clim_2 = readDWD(file_2, varnames = FALSE)
clim_2[, "MESS_DATUM"] = as.Date(clim_2[,"MESS_DATUM"], "%Y-%m-%d")
link_3 = selectDWD(id = statID, res = "daily", var = "kl", per = "historical")
file_3 = dataDWD(link_3, read = FALSE)
clim_3 = readDWD(file_3, varnames = FALSE)
clim_3[, "MESS_DATUM"] = as.Date(clim_3[,"MESS_DATUM"], "%Y-%m-%d")

#Second step: Extract Data from Grid data
link_1 = paste0("https://opendata.dwd.de/climate_environment/CDC/derived_germany/soil/daily/historical/derived_germany_soil_daily_historical_",statID,".txt.gz")
download.file(link_1,paste0("derived_germany_soil_daily_historical_",statID,".txt.gz"))
clim_1 = read.table(gzfile(paste0("derived_germany_soil_daily_historical_",statID,".txt.gz")),sep=";")
colnames(clim_1) = as.character(clim_1[1,])
clim_1 = clim_1[-1,]
clim_1[,"Datum"] = paste0(substr(clim_1[,"Datum"], 1, 4),"-",substr(clim_1[,"Datum"], 5, 6),"-",substr(clim_1[,"Datum"], 7, 8))
clim_1[,"Datum"] = as.Date(clim_1[,"Datum"], "%Y-%m-%d")

#Third step: Formatting 
dw_data_2 = clim_2[,c(which(names(clim_2) %in% c("MESS_DATUM", "RS", "NSH_TAG", "SH_TAG")))]
names(dw_data_2) = c("Date", "Precipitation", "Total Snow Depth", "New Snow Depth")
dw_data_3 = clim_3[,c(which(names(clim_3) %in% c("MESS_DATUM", "FM","TMK")))]
names(dw_data_3) = c("Date", "Wind Speed", "Temperature")
dw_data_1 = clim_1[,c(which(names(clim_1) %in% c("Datum", "VPGB")))]
names(dw_data_1) = c("Date", "Potential Evapotranspiration")
dw_data_1[,"Potential Evapotranspiration"] = as.numeric(dw_data_1[,"Potential Evapotranspiration"])

#Loading and formatting the data (Extended Discharge and Evapotranspiration (PET))
Extended_DR_data = read.csv(tk_choose.files(default = "", caption = "Select csv file"))
names(Extended_DR_data)[1] = "Date"
names(Extended_DR_data)[2] = "Water Level"
names(Extended_DR_data)[3] = "River Discharge"
Extended_DR_data[,"Date"] = as.Date(Extended_DR_data[,"Date"], format = "%d/%m/%Y")

PET = read.csv(tk_choose.files(default = "", caption = "Select csv file"))
names(PET)[1] = "Date"
names(PET)[2] = "Potential Evapotranspiration"
PET[,"Date"] = as.Date(PET[,"Date"], format = "%d/%m/%Y")

#Combining data set
int_date = dw_data_2[1,"Date"]
fin_date = "2021-12-31"
Eva = rbind(PET[which(PET[,"Date"] == int_date):(which(PET[,"Date"] == dw_data_1[1,"Date"])-1),],dw_data_1[1:which(dw_data_1[,"Date"] == fin_date),])
fn_dw = cbind(dw_data_2,Eva[,2],Extended_DR_data[which(Extended_DR_data[,"Date"] == int_date):which(Extended_DR_data[,"Date"] == fin_date),c(2,3)], dw_data_3[which(dw_data_3[,"Date"] == int_date):which(dw_data_3[,"Date"] == fin_date),c(2,3)])
names(fn_dw)[5:9] = c(names(Eva)[2],names(Extended_DR_data)[2:3],names(dw_data_3)[2:3])
rownames(fn_dw) = 1:length(rownames(fn_dw))
fn_dw = mutate(fn_dw,`Wind Speed` = na.approx(`Wind Speed`), `Total Snow Depth` = na.approx(`Total Snow Depth`), `New Snow Depth` = na.approx(`New Snow Depth`))
#fn_dw[,c("Total Snow Depth","New Snow Depth")] = fn_dw[,c("Total Snow Depth","New Snow Depth")]*0.01
#fn_dw[,c("Precipitation", "Potential Evapotranspiration")] = fn_dw[,c("Precipitation", "Potential Evapotranspiration")]*0.001

#### Repeat for each time series to get the plots and the summary statistics ####

#To Plot each series 
#hist(sec_data[,2])
sec_data = fn_dw[,c("Date", "Potential Evapotranspiration")]
x = zoo(sec_data[,2], sec_data[,1])
plot(x, xaxt ="n", xlab = "Date", ylab = "", col = "black")
title(ylab = bquote("Potential Evapotranspiration ( " ~ mm ~ ")"), line = 2.5)
axis(side = 1, at = time(x)[seq(1, length(time(x)), length.out=5)], labels = format(time(x)[seq(1, length(time(x)), length.out=5)], "%d-%m-%Y"),  cex.axis = 0.7)

#To extract summary statistics
summary(sec_data[,2])
stat.desc(sec_data[,2])

#### Analysis on the potential evapotranspiration estimates ####
set_1 = PET[which(PET[,"Date"] == int_date):(which(PET[,"Date"] == dw_data_1[1,"Date"])-1),]
set_2 = dw_data_1[1:which(dw_data_1[,"Date"] == fin_date),]

sew_1 = mutate(set_1, Group = 1)
sew_2 = mutate(set_2, Group = 2)
sew = union(sew_1,sew_2)
sew$Group = factor(sew$Group)

ggplot(sew, aes(x = Date, y = `Potential Evapotranspiration`, group = Group)) + geom_line(aes(color = Group)) + theme_bw() + scale_x_date(date_breaks = "5 year", date_labels =  "%Y") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom", legend.title = element_blank()) + scale_colour_discrete(labels=c(bquote("PE"[HS]), bquote("PE"[AMBAV]))) + labs(y =  bquote("Potential Evapotranspiration ( " ~ mm ~ ")"))

int_com_date = max(range(PET[,"Date"])[1],range(dw_data_1[,"Date"])[1])
fin_com_date = min(range(PET[,"Date"])[2],range(dw_data_1[,"Date"])[2])
dat = mutate(dw_data_1[(which(dw_data_1[,"Date"] == int_com_date)):(which(dw_data_1[,"Date"] == fin_com_date)),], Group = 1)
dat = union(dat, mutate(PET[(which(PET[,"Date"] == int_com_date)):(which(PET[,"Date"] == fin_com_date)),], Group = 2))
dat$Group = factor(dat$Group)

ggplot(dat, aes(x = Date, y = `Potential Evapotranspiration`, group = Group)) + geom_line(aes(color = Group)) + theme_bw() + scale_x_date(date_breaks = "2 year", date_labels =  "%Y") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom", legend.title = element_blank()) + scale_colour_discrete(labels=c(bquote("PE"[AMBAV]), bquote("PE"[HS]))) + labs(y =  bquote("Potential Evapotranspiration ( " ~ mm ~ ")"))

#To plot a matrix of scatter plots
upper.panel = function(x, y){
  points(x,y, pch = 19)
  linear_regression = lm(x~y)
  linear_regression_line = abline(linear_regression, col = "red", lwd = 2)
  cf = round(coef(linear_regression), 2) 
  eq = paste0(names(ser)[1]," = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " ",names(ser)[2] , " ")
  mtext(eq, 3, line=-2)
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
  #r = round(kendall.tau(x, y), digits=2)
  t = round(cor(x, y, method = "pearson"), digits=2)
  #s = round(taildep(x, y, u = 0.9, type = "chi"), digits=2)
  #txt = bquote(tau == .(round(r,2)))
  txt_1 = bquote(rho == .(round(t,2)))
  #txt_2 = bquote(lambda[u] == .(round(s,2)))
  #cex.cor = 0.8/strwidth(txt)
  #text(0.5, 0.7, txt_1, cex = 2.5)
  text(0.5, 0.5, txt_1, cex = 2.5)
  #text(0.5, 0.5, txt, cex = 2.5)
  #text(0.5, 0.3, txt_2, cex = 2.5)
}
ser = as.data.frame(cbind(dat[which(dat$Group == 1), 2], dat[which(dat$Group == 2), 2]))
names(ser) = c("AMBAV", "HS")
pairs(ser, lower.panel = panel.cor,  upper.panel = upper.panel,diag.panel=panel.hist)

#### Target sequences ####
td_data = fn_dw
td_data = td_data %>% mutate(SD = NA, rho = NA, SMWE = NA, PS = NA, CI = NA)
for(h in 1:dim(td_data)[1]){
  
  if(h == 1){td_data[1,"SD"] = 0}else{td_data[h,"SD"] = td_data[h-1,"Total Snow Depth"] + td_data[h,"New Snow Depth"] - td_data[h,"Total Snow Depth"]}
  td_data[h,"rho"] = (90 + 130*(sqrt(td_data[h,"Total Snow Depth"])))*(1.5 + 0.17*(pracma::nthroot(td_data[h,"Temperature"], 3)))*(1 + 0.1*(sqrt(td_data[h,"Wind Speed"])))
  td_data[h,"SMWE"] = td_data[h,"SD"]*td_data[h,"rho"]/100
  td_data[h,"PS"] = td_data[h,"Precipitation"] - td_data[h,"Potential Evapotranspiration"]
  td_data[h,"CI"] = max((td_data[h,"SMWE"] + td_data[h,"PS"]),0)
    
}

sec_data = td_data[,c("Date", "River Discharge")]
x = zoo(sec_data[,2], sec_data[,1])
plot(x, xaxt ="n", xlab = "Date", ylab = "", col = "black")
title(ylab = bquote("T"[RD] ~ " (" ~ m^3*s^-1 ~ ")"), line = 2.5)
axis(side = 1, at = time(x)[seq(1, length(time(x)), length.out=5)], labels = format(time(x)[seq(1, length(time(x)), length.out=5)], "%d-%m-%Y"),  cex.axis = 0.7)

#Extreme Points
Fl = sec_data[which(sec_data[,2] %in% sort(sec_data[,2],decreasing=T)[1:5]),]
Fl = Fl[order(Fl[,2], decreasing = T),]
points(Fl, col = "red")
td_data[which(td_data[,"Date"] %in% Fl[,"Date"]),c("Date", "CI", "River Discharge")]

#To extract summary statistics
summary(sec_data[,2])
stat.desc(sec_data[,2])

#To plot particular hydrographs
hydrograph = function (input = matrix(ncol = 2, nrow = 2), streamflow = input[,2], timeSeries = input[, 1], streamflow2 = NULL, precip = NULL, begin = 1, endindex = length(streamflow), P.units = "", S.units = "", S1.col = "black", S2.col = "red", stream.label = "Streamflow", flow.label = "Precipitation", streamflow3=NULL,streamflow4=NULL, precip2=NULL){
  
  if (is.null(streamflow2) & (ncol(input) > 3)) 
    streamflow2 <- input[, 4]
  if (is.null(precip) & (ncol(input) > 2)) {
    precip <- input[, 2]
    streamflow <- input[, 3]
  }
  if (!is.null(precip))  {
    par(mar = c(3, 5, 1, 4))
    barplot(precip[begin:endindex], yaxt = "n", space = NULL, 
            ylim = rev(c(0, 4 * max(na.omit(precip[begin:endindex])))), 
            xaxt = "n")
    axis(side = 3, pos = 0, tck = 0,xaxt = "n")
    axis(side = 4, at = seq(0, floor(max(na.omit(precip[begin:endindex])) + 
                                       1), length = (1 + ifelse(floor(max(na.omit(precip[begin:endindex])) + 
                                                                        1) < 10, floor(max(na.omit(precip[begin:endindex])) + 1), 
                                                                4))), labels = as.integer(seq(0, floor(max(na.omit(precip[begin:endindex])) + 
                                                                                                         1), length = (1 + ifelse(floor(max(na.omit(precip[begin:endindex])) + 
                                                                                                                                          1) < 10, floor(max(na.omit(precip[begin:endindex])) + 1), 
                                                                                                                                  4)))))
    if (P.units=="") {
      mtext(paste(flow.label, P.units), 4, line = 2, cex = 0.9, adj = 1)
    } else  mtext(paste(flow.label, " (", P.units, ")", sep=""), 4, line = 2, cex = 0.9, adj = 1)
    par(new = T)
  }
  if (!is.null(precip2)){
    barplot(precip2[begin:endindex], yaxt = "n", space = NULL, col="pink", 
            ylim = rev(c(0, 4 * max(na.omit(precip[begin:endindex])))), 
            xaxt = "n")
    par(new = T)
  }
  
  plot(streamflow[begin:endindex], col = S1.col, type = "l", 
       lwd = 1, ylab = stream.label, xaxt = "n", xlab = "date", 
       ylim = c(0, 1.2 * max(na.omit(streamflow[begin:endindex]), na.omit(streamflow2[begin:endindex]))), 
       axes = FALSE)
  #mtext (expression(paste("                              ", " (" , m^3/s, ")", sep="")), 2,3)
  if (S.units=="m3/s" | S.units=="m3s"){
    mtext (expression(paste(" (" , m^3/s, ")", sep="")), 2,1.5)
  } else if (S.units=="ft3/s" | S.units=="ft3s") {
    mtext (expression(paste(" (" , ft^3/s, ")", sep="")), 2,1.5)
  } else if (S.units!="") mtext (paste(" (" , S.units, ")", sep=""), 2,1.5)
  lines(streamflow2[begin:endindex], col = S2.col, lwd = 1, 
        lty = 2, xaxt = "n")
  if (!is.null(streamflow3)){
    lines(streamflow3[begin:endindex], col = "blue", lwd = 1,   ##potential for more streamflows
          lty = 3, xaxt = "n")
  }
  if (!is.null(streamflow4)){
    lines(streamflow4[begin:endindex], col = "green", lwd = 1,   ##potential for more streamflows
          lty = 4, xaxt = "n")
  }
  axis(side = 1, at = seq(1, (endindex - begin + 1), length = 14), 
       pos = 0, labels = format(timeSeries[seq(begin, endindex, 
                                               length = 14)], "%d-%m-%y"))
  
  axis(side = 2, pos = 0)
}
ut = 1995
lag = 15
pt = td_data[which(as.numeric(format(td_data[,"Date"], format = "%Y")) == ut),c("Date", "CI", "River Discharge")]
hydrograph(pt, stream.label = names(pt)[3], flow.label = "Total Water")
title(xlab = "Date", line = 1)
abline(v = which(pt[,3] == max(pt[,3])), lty = 6)
abline(v = which(pt[,3] == max(pt[,3])) - (lag+1) + which(pt[(which(pt[,3] == max(pt[,3]))-lag):(which(pt[,3] == max(pt[,3]))+lag),2] == max(pt[(which(pt[,3] == max(pt[,3]))-lag):(which(pt[,3] == max(pt[,3]))+lag),2])), lty = 6)

#Create a data.frame with discharge and rolling sum of total water for selected rolling period
rolday = 7
RS = c(rep(0,rolday-1), roll_sum(x = td_data[,"CI"], n = rolday, by = 1))
fd_data = cbind(td_data,RS)[7:nrow(td_data),]

sec_data = fd_data[,c("Date", "RS")]
x = zoo(sec_data[,2], sec_data[,1])
plot(x, xaxt ="n", xlab = "Date", ylab = "", col = "black")
title(ylab = bquote("T"[RS] ~ " (" ~ mm ~ ")"), line = 2.5)
axis(side = 1, at = time(x)[seq(1, length(time(x)), length.out=5)], labels = format(time(x)[seq(1, length(time(x)), length.out=5)], "%d-%m-%Y"),  cex.axis = 0.7)

#To extract summary statistics
summary(sec_data[,2])
stat.desc(sec_data[,2])

#To plot a matrix of scatter plots
upper.panel = function(x, y){
  points(x,y, pch = 19)}
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
  t = round(cor(x, y, method = "spearman"), digits=2)
  s = round(taildep(x, y, u = 0.9, type = "chi"), digits=2)
  txt = bquote(tau == .(round(r,2)))
  txt_1 = bquote(rho == .(round(t,2)))
  txt_2 = bquote(lambda[u] == .(round(s,2)))
  #cex.cor = 0.8/strwidth(txt)
  text(0.5, 0.7, txt_1, cex = 2.5)
  text(0.5, 0.5, txt, cex = 2.5)
  text(0.5, 0.3, txt_2, cex = 2.5)
}
pairs(fd_data[,c("River Discharge", "RS")], lower.panel = panel.cor,  upper.panel = upper.panel,diag.panel=panel.hist)

#To Group data
Kendall_Lag = function(Data, Lags = seq(-6, 6, 1), PLOT = TRUE, GAP = 0.1){
  
  Lag = function(x, k){
    
    if(k > 0){
      
      return(c(rep(NA,k), x)[1:length(x)])
      
    }else{
      
      return(c(x[(-k + 1):length(x)], rep(NA, -k)))
      
    }
  }
  
  correlation = function(Data_Cor,lag){
    
    for(i in 2:(ncol(Data_Cor))){
      
      Data_Cor[,i] = Lag(Data_Cor[,i], - lag[i-1])
      
    }
    
    Data_Cor = na.omit(Data_Cor)
    return(corKendall(Data_Cor[2:ncol(Data_Cor)], checkNA = FALSE)[which(lower.tri(cor(Data_Cor[2:ncol(Data_Cor)])) == T)])
    
  }
  
  n = ncol(Data) - 1
  
  if(n == 3){
    
    Var1_Var2 = numeric(length(Lags))
    Var2_Var3 = numeric(length(Lags))
    Var1_Var3 = numeric(length(Lags))
    
    for(i in 1:length(Lags)){
      
      Var1_Var2[i] = correlation(Data_Cor = Data, c(Lags[i],0,0))[1]
      Var2_Var3[i] = correlation(Data_Cor = Data, c(0,0,Lags[i]))[2]
      Var1_Var3[i] = correlation(Data_Cor = Data, c(0,0,Lags[i]))[3]
      
    }
    
    correlation.test = function(Data_Cor,lag){
      
      for(i in 2:(ncol(Data_Cor))){
        
        Data_Cor[,i] = Lag(Data_Cor[,i],-lag[i-1])
        
      }
      
      Data_Cor = na.omit(Data_Cor)
      return(cor.test.default(Data_Cor[,2],Data_Cor[,3],method = "kendall")$p.value)
      
    }
    
    Var1_Var2_Test = numeric(length(Lags))
    Var2_Var3_Test = numeric(length(Lags))
    Var1_Var3_Test = numeric(length(Lags))
    
    for(i in 1:length(Lags)){
      
      Var1_Var2_Test[i] = correlation.test(Data_Cor = Data[,-4],c(Lags[i],0))
      Var2_Var3_Test[i] = correlation.test(Data_Cor = Data[,-3],c(0,Lags[i]))
      Var1_Var3_Test[i] = correlation.test(Data_Cor = Data[,-2],c(0,Lags[i]))
      
    }
    
    if(PLOT==TRUE){
      
      yx = max(c(Var1_Var2,Var2_Var3,Var1_Var3))-min(c(Var1_Var2,Var2_Var3,Var1_Var3))
      plot(Lags,Var1_Var2, ylim = c(min(c(Var1_Var2,Var2_Var3,Var1_Var3))-GAP*yx,max(c(Var1_Var2,Var2_Var3,Var1_Var3)))+GAP*yx,type='l',xlab="Lag",ylab=expression(paste("Kendall's "*tau*' coefficient')),cex.lab=1,cex.axis=0.8,lwd=2.5)
      abline(h = 0, lty = 2)
      lines(Lags, Var2_Var3, col=2)
      lines(Lags, Var1_Var3, col=3)
      
      points(Lags,Var1_Var2, pch = ifelse(Var1_Var2_Test < 0.05, 16, 21), bg = ifelse(Var1_Var2_Test < 0.05,1, "White"), cex=1.5)
      points(Lags,Var2_Var3, pch = ifelse(Var2_Var3_Test < 0.05, 16, 21), col = 2, bg=ifelse(Var2_Var3_Test < 0.05,2, "White"), cex=1.5)
      points(Lags,Var1_Var3, pch = ifelse(Var1_Var3_Test < 0.05, 16, 21), col = 3, bg=ifelse(Var1_Var3_Test < 0.05,3, "White"), cex=1.5)
      legend("topleft", c(paste(colnames(Data)[2], " & ", colnames(Data)[3], sep=""),
                          paste(colnames(Data)[2], " & ", colnames(Data)[4], sep=""),
                          paste(colnames(Data)[3], " & ", colnames(Data)[4], sep="")),
             bty = "n", lwd = 2.5, col = c(1,2,3))
    }
    
    Value = list()
    Value[[paste(names(Data)[2], '_', names(Data)[3], sep="")]]= Var1_Var2
    Value[[paste(names(Data)[3], '_', names(Data)[4], sep="")]]= Var2_Var3
    Value[[paste(names(Data)[2], '_', names(Data)[4], sep="")]]= Var1_Var3
    
    Test = list()
    Test[[paste(names(Data)[2], '_', names(Data)[3], '_Test',sep="")]]= Var1_Var2_Test
    Test[[paste(names(Data)[3], '_', names(Data)[4], '_Test',sep="")]]= Var2_Var3_Test
    Test[[paste(names(Data)[2], '_', names(Data)[4], '_Test',sep="")]]= Var1_Var3_Test
  }
  
  if(n == 2){
    
    Var1_Var2 = numeric(length(Lags))
    
    for (i in 1:length(Lags)) {
      
      Var1_Var2[i] = correlation(Data_Cor = Data, c(Lags[i], 0))[1]
      
    }
    
    correlation.test = function(Data_Cor, lag) {
      
      for (i in 2:(ncol(Data_Cor))) {
        
        Data_Cor[, i] = Lag(Data_Cor[, i], -lag[i - 1])
        
      }
      
      Data_Cor = na.omit(Data_Cor)
      return(cor.test.default(Data_Cor[, 2], Data_Cor[, 3], method = "kendall")$p.value)
      
    }
    
    Var1_Var2_Test = numeric(length(Lags))
    
    for (i in 1:length(Lags)) {
      
      Var1_Var2_Test[i] = correlation.test(Data_Cor = Data[,-4], c(Lags[i], 0))
      
    }
    
    if (PLOT == TRUE) {
      
      yx = max(c(Var1_Var2)) - min(c(Var1_Var2))
      
      plot(Lags, Var1_Var2, ylim = c(min(c(Var1_Var2)) - GAP * yx, max(c(Var1_Var2)) + GAP * yx), type = "l", xlab = "Lag",
           ylab = expression(paste("Kendall's " * tau * " coefficient")),
           cex.lab = 1.45, cex.axis = 1.45, lwd = 2.5)
      abline(h = 0, lty = 2)
      points(Lags, Var1_Var2, pch = ifelse(Var1_Var2_Test <  0.05, 16, 21), bg = ifelse(Var1_Var2_Test < 0.05,1, "White"), cex = 1.5)
      legend("topleft", c(paste(colnames(Data)[2], " & ", colnames(Data)[3],sep = "")), bty = "n", lwd = 2.5, col = c(1, 2, 3))
      
    }
    
    Value = list()
    Value[[paste(names(Data)[2], '_', names(Data)[3], sep="")]]= Var1_Var2
    
    Test = list()
    Test[[paste(names(Data)[2], '_', names(Data)[3], '_Test', sep="")]]= Var1_Var2_Test
    
  }
  
  res = list("Value" = Value, "Test" = Test)
  return(res)
  
}
Pearson_Lag = function(Data, Lags = seq(-6, 6, 1), PLOT = TRUE, GAP = 0.1){
  
  Lag = function(x, k){
    
    if(k > 0){
      
      return(c(rep(NA,k), x)[1:length(x)])
      
    }else{
      
      return(c(x[(-k + 1):length(x)], rep(NA, -k)))
      
    }
  }
  
  correlation = function(Data_Cor,lag){
    
    for(i in 2:(ncol(Data_Cor))){
      
      Data_Cor[,i] = Lag(Data_Cor[,i], - lag[i-1])
      
    }
    
    Data_Cor = na.omit(Data_Cor)
    return(cor(Data_Cor[2:ncol(Data_Cor)], method = "pearson")[which(lower.tri(cor(Data_Cor[2:ncol(Data_Cor)])) == T)])
    
  }
  
  n = ncol(Data) - 1
  
  if(n == 3){
    
    Var1_Var2 = numeric(length(Lags))
    Var2_Var3 = numeric(length(Lags))
    Var1_Var3 = numeric(length(Lags))
    
    for(i in 1:length(Lags)){
      
      Var1_Var2[i] = correlation(Data_Cor = Data, c(Lags[i],0,0))[1]
      Var2_Var3[i] = correlation(Data_Cor = Data, c(0,0,Lags[i]))[2]
      Var1_Var3[i] = correlation(Data_Cor = Data, c(0,0,Lags[i]))[3]
      
    }
    
    correlation.test = function(Data_Cor,lag){
      
      for(i in 2:(ncol(Data_Cor))){
        
        Data_Cor[,i] = Lag(Data_Cor[,i],-lag[i-1])
        
      }
      
      Data_Cor = na.omit(Data_Cor)
      return(cor.test.default(Data_Cor[,2],Data_Cor[,3],method = "pearson")$p.value)
      
    }
    
    Var1_Var2_Test = numeric(length(Lags))
    Var2_Var3_Test = numeric(length(Lags))
    Var1_Var3_Test = numeric(length(Lags))
    
    for(i in 1:length(Lags)){
      
      Var1_Var2_Test[i] = correlation.test(Data_Cor = Data[,-4],c(Lags[i],0))
      Var2_Var3_Test[i] = correlation.test(Data_Cor = Data[,-3],c(0,Lags[i]))
      Var1_Var3_Test[i] = correlation.test(Data_Cor = Data[,-2],c(0,Lags[i]))
      
    }
    
    if(PLOT==TRUE){
      
      yx = max(c(Var1_Var2,Var2_Var3,Var1_Var3))-min(c(Var1_Var2,Var2_Var3,Var1_Var3))
      plot(Lags,Var1_Var2, ylim = c(min(c(Var1_Var2,Var2_Var3,Var1_Var3))-GAP*yx,max(c(Var1_Var2,Var2_Var3,Var1_Var3)))+GAP*yx,type='l',xlab="Lag",ylab=expression(paste("Pearson's coefficient")),cex.lab=1,cex.axis=1,lwd=2.5)
      abline(h = 0, lty = 2)
      lines(Lags, Var2_Var3, col=2)
      lines(Lags, Var1_Var3, col=3)
      
      points(Lags,Var1_Var2, pch = ifelse(Var1_Var2_Test < 0.05, 16, 21), bg = ifelse(Var1_Var2_Test < 0.05,1, "White"), cex=1.5)
      points(Lags,Var2_Var3, pch = ifelse(Var2_Var3_Test < 0.05, 16, 21), col = 2, bg=ifelse(Var2_Var3_Test < 0.05,2, "White"), cex=1.5)
      points(Lags,Var1_Var3, pch = ifelse(Var1_Var3_Test < 0.05, 16, 21), col = 3, bg=ifelse(Var1_Var3_Test < 0.05,3, "White"), cex=1.5)
      legend("topleft", c(paste(colnames(Data)[2], " & ", colnames(Data)[3], sep=""),
                          paste(colnames(Data)[2], " & ", colnames(Data)[4], sep=""),
                          paste(colnames(Data)[3], " & ", colnames(Data)[4], sep="")),
             bty = "n", lwd = 2.5, col = c(1,2,3))
    }
    
    Value = list()
    Value[[paste(names(Data)[2], '_', names(Data)[3], sep="")]]= Var1_Var2
    Value[[paste(names(Data)[3], '_', names(Data)[4], sep="")]]= Var2_Var3
    Value[[paste(names(Data)[2], '_', names(Data)[4], sep="")]]= Var1_Var3
    
    Test = list()
    Test[[paste(names(Data)[2], '_', names(Data)[3], '_Test',sep="")]]= Var1_Var2_Test
    Test[[paste(names(Data)[3], '_', names(Data)[4], '_Test',sep="")]]= Var2_Var3_Test
    Test[[paste(names(Data)[2], '_', names(Data)[4], '_Test',sep="")]]= Var1_Var3_Test
  }
  
  if(n == 2){
    
    Var1_Var2 = numeric(length(Lags))
    
    for (i in 1:length(Lags)) {
      
      Var1_Var2[i] = correlation(Data_Cor = Data, c(Lags[i], 0))[1]
      
    }
    
    correlation.test = function(Data_Cor, lag) {
      
      for (i in 2:(ncol(Data_Cor))) {
        
        Data_Cor[, i] = Lag(Data_Cor[, i], -lag[i - 1])
        
      }
      
      Data_Cor = na.omit(Data_Cor)
      return(cor.test.default(Data_Cor[, 2], Data_Cor[, 3], method = "pearson")$p.value)
      
    }
    
    Var1_Var2_Test = numeric(length(Lags))
    
    for (i in 1:length(Lags)) {
      
      Var1_Var2_Test[i] = correlation.test(Data_Cor = Data[,-4], c(Lags[i], 0))
      
    }
    
    if (PLOT == TRUE) {
      
      yx = max(c(Var1_Var2)) - min(c(Var1_Var2))
      
      plot(Lags, Var1_Var2, ylim = c(min(c(Var1_Var2)) - GAP * yx, max(c(Var1_Var2)) + GAP * yx), type = "l", xlab = "Lag",
           ylab = expression(paste("Pearson's coefficient")),
           cex.lab = 1.45, cex.axis = 1.45, lwd = 2.5)
      abline(h = 0, lty = 2)
      points(Lags, Var1_Var2, pch = ifelse(Var1_Var2_Test <  0.05, 16, 21), bg = ifelse(Var1_Var2_Test < 0.05,1, "White"), cex = 1.5)
      legend("topleft", c(paste(colnames(Data)[2], " & ", colnames(Data)[3],sep = "")), bty = "n", lwd = 2.5, col = c(1, 2, 3))
      
    }
    
    Value = list()
    Value[[paste(names(Data)[2], '_', names(Data)[3], sep="")]]= Var1_Var2
    
    Test = list()
    Test[[paste(names(Data)[2], '_', names(Data)[3], '_Test', sep="")]]= Var1_Var2_Test
    
  }
  
  res = list("Value" = Value, "Test" = Test)
  return(res)
  
}
cor.test.default = function(x, y, alternative = c("two.sided", "less", "greater"), method = c("pearson", "kendall", "spearman"), exact = NULL,conf.level = 0.95, continuity = FALSE, ...){
  
  alternative = match.arg(alternative)
  method = match.arg(method)
  DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  if(length(x) != length(y)) stop("'x' and 'y' must have the same length")
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  if(!is.numeric(y)) stop("'y' must be a numeric vector")
  
  OK = complete.cases(x, y)
  x = x[OK]
  y = y[OK]
  n = length(x)
  
  NVAL = 0
  conf.int = FALSE
  
  if(method == "pearson") {
    if(n < 3L)
      stop("not enough finite observations")
    method <- "Pearson's product-moment correlation"
    names(NVAL) <- "correlation"
    r <- cor(x, y)
    df <- n - 2L
    ESTIMATE <- c(cor = r)
    PARAMETER <- c(df = df)
    STATISTIC <- c(t = sqrt(df) * r / sqrt(1 - r^2))
    if(n > 3) { ## confidence int.
      if(!missing(conf.level) &&
         (length(conf.level) != 1 || !is.finite(conf.level) ||
          conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
      conf.int <- TRUE
      z <- atanh(r)
      sigma <- 1 / sqrt(n - 3)
      cint <-
        switch(alternative,
               less = c(-Inf, z + sigma * qnorm(conf.level)),
               greater = c(z - sigma * qnorm(conf.level), Inf),
               two.sided = z +
                 c(-1, 1) * sigma * qnorm((1 + conf.level) / 2))
      cint <- tanh(cint)
      attr(cint, "conf.level") <- conf.level
    }
    PVAL <- switch(alternative,
                   "less" = pt(STATISTIC, df),
                   "greater" = pt(STATISTIC, df, lower.tail=FALSE),
                   "two.sided" = 2 * min(pt(STATISTIC, df),
                                         pt(STATISTIC, df, lower.tail=FALSE)))
  }
  else {
    
    if(n < 2) stop("not enough finite observations")
    
    PARAMETER = NULL
    TIES = (min(length(unique(x)), length(unique(y))) < n)
    
    if(method == "kendall") {
      method = "Kendall's rank correlation tau"
      names(NVAL) = "tau"
      r = kendall.tau(x, y)
      ESTIMATE = c(tau = r)
      
      if(!is.finite(ESTIMATE)){  # all x or all y the same
        
        ESTIMATE[] <- NA
        STATISTIC <- c(T = NA)
        PVAL <- NA
        
      }else {
        if(is.null(exact))
          exact <- (n < 50)
        if(exact && !TIES){
          q <- round((r + 1) * n * (n - 1) / 4)
          STATISTIC <- c(T = q)
          pkendall <- function(q,n){SuppDists::pKendall(q, n)}
          PVAL <-
            switch(alternative,
                   "two.sided" = {
                     if(q > n * (n - 1) / 4)
                       p <- 1 - pkendall(q - 1, n)
                     else
                       p <- pkendall(q, n)
                     min(2 * p, 1)
                   },
                   "greater" = 1 - pkendall(q - 1, n),
                   "less" = pkendall(q, n))
        } else {
          
          xties = table(x[duplicated(x)]) + 1
          yties = table(y[duplicated(y)]) + 1
          T0 = n * (n - 1)/2
          T1 = sum(xties * (xties - 1))/2
          T2 = sum(yties * (yties - 1))/2
          S = r * sqrt((T0 - T1) * (T0 - T2))
          v0 = n * (n - 1) * (2 * n + 5)
          vt = sum(xties * (xties - 1)*(2 * xties + 5))
          vu = sum(yties * (yties - 1)*(2 * yties + 5))
          v1 = sum(xties * (xties - 1))*sum(yties*(yties - 1))
          v2 = sum(xties * (xties - 1)*(xties - 2))*sum(yties * (yties - 1) * (yties - 2))
          
          var_S = (v0 - vt - vu) / 18 +
            v1 / (2 * n * (n - 1)) +
            v2 / (9 * n * (n - 1) * (n - 2))
          
          if(exact && TIES)
            warning("Cannot compute exact p-value with ties")
          if (continuity) S <- sign(S) * (abs(S) - 1)
          STATISTIC <- c(z = S / sqrt(var_S))
          PVAL <- switch(alternative,
                         "less" = pnorm(STATISTIC),
                         "greater" = pnorm(STATISTIC, lower.tail=FALSE),
                         "two.sided" = 2 * min(pnorm(STATISTIC),
                                               pnorm(STATISTIC, lower.tail=FALSE)))
        }
      }
    } else {
      method <- "Spearman's rank correlation rho"
      if (is.null(exact))
        exact <- TRUE
      names(NVAL) <- "rho"
      r <- cor(rank(x), rank(y))
      ESTIMATE <- c(rho = r)
      if(!is.finite(ESTIMATE)) {  # all x or all y the same
        ESTIMATE[] <- NA
        STATISTIC <- c(S = NA)
        PVAL <- NA
      }
      else {
        ## Use the test statistic S = sum(rank(x) - rank(y))^2
        ## and AS 89 for obtaining better p-values than via the
        ## simple normal approximation.
        ## In the case of no ties, S = (1-rho) * (n^3-n)/6.
        pspearman <- function(q, n, lower.tail = TRUE) {
          if(n <= 1290 && exact) # n*(n^2 - 1) does not overflow
            .Call(C_pRho, round(q) + 2*lower.tail, n, lower.tail)
          else { # for large n: asymptotic t_{n-2}
            den <- (n*(n^2-1))/6 # careful for overflow
            ## Kendall et all (1939) p. 260
            if (continuity) den <- den + 1
            r <- 1 - q/den
            pt(r / sqrt((1 - r^2)/(n-2)), df = n-2,
               lower.tail = !lower.tail)
          }
        }
        q <- (n^3 - n) * (1 - r) / 6
        STATISTIC <- c(S = q)
        if(TIES && exact){
          exact <- FALSE
          warning("Cannot compute exact p-value with ties")
        }
        PVAL <-
          switch(alternative,
                 "two.sided" = {
                   p <- if(q > (n^3 - n) / 6)
                     pspearman(q, n, lower.tail = FALSE)
                   else
                     pspearman(q, n, lower.tail = TRUE)
                   min(2 * p, 1)
                 },
                 "greater" = pspearman(q, n, lower.tail = TRUE),
                 "less" = pspearman(q, n, lower.tail = FALSE))
      }
    }
  }
  
  RVAL <- list(statistic = STATISTIC,
               parameter = PARAMETER,
               p.value = as.numeric(PVAL),
               estimate = ESTIMATE,
               null.value = NVAL,
               alternative = alternative,
               method = method,
               data.name = DNAME)
  if(conf.int)
    RVAL <- c(RVAL, list(conf.int = cint))
  class(RVAL) <- "htest"
  RVAL
}

gd_data = fd_data[,c("Date", "River Discharge", "RS")]
gd_data$Date = format(gd_data$Date, format = "%Y-%m") 
gd_data = gd_data %>% group_by(Date) %>% summarise(MRD = mean(`River Discharge`), MRS = mean(RS))
gd_data = as.data.frame(gd_data)
gd_data$Date = as.Date(as.yearmon(gd_data$Date))

sec_data = gd_data[,c("Date", "MRS")]
x = zoo(sec_data[,2], sec_data[,1])
plot(x, xaxt ="n", xlab = "Date", ylab = "", col = "black")
title(ylab = bquote("T"[MRS] ~ " (" ~ mm ~ ")"), line = 2.5)
axis(side = 1, at = time(x)[seq(1, length(time(x)), length.out=5)], labels = format(time(x)[seq(1, length(time(x)), length.out=5)], "%d-%m-%Y"),  cex.axis = 0.7)

#Extreme Points
Fl = sec_data[which(sec_data[,2] %in% sort(sec_data[,2],decreasing=T)[1:5]),]
Fl = Fl[order(Fl[,2], decreasing = T),]
points(Fl, col = "red")
gd_data[which(gd_data[,"Date"] %in% Fl[,"Date"]),c("Date", "MRS", "MRD")]

#To extract summary statistics
summary(sec_data[,2])
stat.desc(sec_data[,2])

pairs(gd_data[,c("MRD", "MRS")], lower.panel = panel.cor,  upper.panel = upper.panel,diag.panel=panel.hist)
Kendall_Lag(gd_data[,c("Date", "MRD", "MRS")], Lag = seq(-6,6,1))
Pearson_Lag(gd_data[,c("Date", "MRD", "MRS")], Lag = seq(-6,6,1))

td_data[,c("Date", "River Discharge")] %>% mutate(Month = format(Date, format = "%m")) %>% group_by(Month) %>% summarise(mean = mean(`River Discharge`)) %>%
  mutate(month = month.name[as.numeric(Month)]) %>% ggplot(aes(x = month, y = mean, fill = "#00abff")) + labs(x = "Month", y = bquote("Daily Mean River Discharge" ~ " (" ~ m^3*s^-1 ~ ")")) +
  geom_col() + theme_bw() + scale_x_discrete(limits = month.name) + theme(legend.position="none") + scale_fill_manual(values = "#00BE67") 

dew = td_data[,c("Date", "River Discharge")] %>% mutate(Year = format(Date, format = "%Y")) %>% group_by(Year) %>% summarise(max = max(`River Discharge`), Date = "")
dew = as.data.frame(dew)
for(h in 1:nrow(dew)){
  
  dew[h,3] = td_data[which(format(td_data[,"Date"],format = "%Y") == as.numeric(dew[h,1])),c("Date","River Discharge")][which(td_data[which(format(td_data[,"Date"],format = "%Y") == as.numeric(dew[h,1])),"River Discharge"] == as.numeric(dew[h,2])),"Date"]
  
}
dew[,"Date"] = as.Date(as.numeric(dew[,"Date"]))
dew = dew %>% mutate(Month = format(Date,format = "%m"))
100*length(which(dew[,"Month"] %in% c("11", "12", "01", "02", "03")))/nrow(dew)
100*length(which(dew[,"Month"] %in% c("11", "12", "01", "02", "03")))/nrow(dew)
100*length(which(dew[,"Month"] %in% c("05", "06", "07", "08", "09")))/nrow(dew)

#To change to water year
data = gd_data
if(length(which(data[,1] == paste0(unique(year(data$Date))[1],"-10-01") | data[,1] == paste0(unique(year(data$Date))[2],"-10-01")))<2){
  
  breaks = seq(data[which(data[,1] == paste0(unique(year(data$Date))[1],"-10-01") | data[,1] == paste0(unique(year(data$Date))[2],"-10-01"))[1],1], length = (length(unique(year(data$Date)))), by="year")
  data$hydroYear = cut(data$Date, breaks, labels = unique(year(data$Date))[2]:(unique(year(data$Date))[length(unique(year(data$Date)))]))
  
}else{
  
  breaks = seq(data[which(data[,1] == paste0(unique(year(data$Date))[1],"-10-01") | data[,1] == paste0(unique(year(data$Date))[2],"-10-01"))[1],1], length = (length(unique(year(data$Date))) + 1), by="year")
  data$hydroYear = cut(data$Date, breaks, labels = unique(year(data$Date))[1]:(unique(year(data$Date))[length(unique(year(data$Date)))]))
  
}
new_data = data %>% dplyr::select(names(data)[length(names(data))], everything())
new_data = na.omit(new_data)
for(i in 1:length(unique(new_data$hydroYear))){
  
  if(length(which(new_data$hydroYear == unique(new_data$hydroYear)[i])) < 12){
    
    new_data = new_data[-which(new_data$hydroYear == unique(new_data$hydroYear)[i]),]
    
  }
  
}

sec_data = new_data[,c("Date", "MRS")]
x = zoo(sec_data[,2], sec_data[,1])
plot(x, xaxt ="n", xlab = "Date", ylab = "", col = "black")
title(ylab = bquote("T[MRS] ( " ~ m^3*s^-1 ~ ")"), line = 2.5)
axis(side = 1, at = time(x)[seq(1, length(time(x)), length.out=5)], labels = format(time(x)[seq(1, length(time(x)), length.out=5)], "%d-%m-%Y"),  cex.axis = 0.7)

#Extreme Points
Fl = sec_data[which(sec_data[,2] %in% sort(sec_data[,2],decreasing=T)[1:5]),]
Fl = Fl[order(Fl[,2], decreasing = T),]
points(Fl, col = "red")
new_data[which(new_data[,"Date"] %in% Fl[,"Date"]),c("Date", "MRS", "MRD")]

#Mean Line
#eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
reg = lm(MRS ~ Date, data = new_data)
coeff = coefficients(reg)
abline(reg, col="blue")

#To extract summary statistics
summary(sec_data[,2])
stat.desc(sec_data[,2])

#To Save
write.csv(new_data, paste0(getwd(),"\\Summary_Final.csv"), row.names=FALSE)

#Assumption Models Autocorrelation function vs Periodogram 
new_data$hydroYear = factor(new_data$hydroYear)
new_data %>% group_by(hydroYear) %>%  mutate(month = month.name[as.numeric(format(Date, format = "%m"))]) %>%
  ggplot(aes(month,MRD)) + geom_boxplot(fill = "#2BE381") + theme_bw() + scale_x_discrete(limits = month.name) + labs(x = "Month")

new_data %>% group_by(hydroYear) %>%  mutate(month = month.name[as.numeric(format(Date, format = "%m"))]) %>%
  ggplot(aes(month,MRS)) + geom_boxplot(fill = "#2BE381") + theme_bw() + scale_x_discrete(limits = month.name) + labs(x = "Month")

#Hypothesis test (KPSS test) H0 = trend stationary vs H1 unit root
tseries::pp.test(new_data$MRD, alternative = "stationary")
tseries::pp.test(new_data$MRS, alternative = "stationary")

#Trend Test (MannKendall trend test) -> Hypothesis H0: There is no discernible pattern in the data & H1: There is a trend in the data
MannKendall(new_data$MRD)
MannKendall(new_data$MRS)
SeasonalMannKendall(ts(new_data$MRD, start = c("1969",10), freq = 12))
SeasonalMannKendall(ts(new_data$MRS, start = c("1969",10), freq = 12))

tmod = lm(formula = MRD ~ Date, data = new_data)
dwtest(formula = tmod,  alternative = "two.sided")

acf(new_data$MRS)
acf(new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),"MRS"], lag.max = 20)
acf(new_data[which(format(new_data$Date, format = "%m") %in% c("04","05", "06", "07", "08", "09")),"MRS"], lag.max = 20)

acf(new_data$MRD)
acf(new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),"MRD"], lag.max = 20)
acf(new_data[which(format(new_data$Date, format = "%m") %in% c("04","05", "06", "07", "08", "09")),"MRD"], lag.max = 20)

s = spectrum(new_data$MRS)
abline(v = s$freq[which(s$spec == max(s$spec))], col = "blue", lty = 2)
1/s$freq[which(s$spec == max(s$spec))]

s_1 = spectrum(new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),"MRS"])
abline(v = s_1$freq[which(s_1$spec == max(s_1$spec))], col = "blue", lty = 2)
1/s_1$freq[which(s_1$spec == max(s_1$spec))]

s_2 = spectrum(new_data[which(format(new_data$Date, format = "%m") %in% c("04","05", "06", "07", "08", "09")),"MRS"])
abline(v = s_2$freq[which(s_2$spec == max(s_2$spec))], col = "blue", lty = 2)
1/s_2$freq[which(s_2$spec == max(s_2$spec))]

t = spectrum(new_data$MRD)
abline(v = t$freq[which(t$spec == max(t$spec))], col = "blue", lty = 2)
1/t$freq[which(t$spec == max(t$spec))]

t_1 = spectrum(new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),"MRD"])
abline(v = t_1$freq[which(t_1$spec == max(t_1$spec))], col = "blue", lty = 2)
1/t_1$freq[which(t_1$spec == max(t_1$spec))]

t_2 = spectrum(new_data[which(format(new_data$Date, format = "%m") %in% c("04","05", "06", "07", "08", "09")),"MRD"])
abline(v = t_2$freq[which(t_2$spec == max(t_2$spec))], col = "blue", lty = 2)
1/t_2$freq[which(t_2$spec == max(t_2$spec))]

#To extract summary statistics
x = new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),"MRS"]
summary(x)
stat.desc(x)

#Multivariate Analysis: Acf and periodogram
#ts = new_data
#ts = new_data[which(format(new_data$Date, format = "%m") %in% c("04","05", "06", "07", "08", "09")),]
ts = new_data[which(format(new_data$Date, format = "%m") %in% c("10", "11","12", "01", "02", "03")),]
par(mfrow = c(2,2))
acf(ts$MRD, main = "MRD", lag = 30)
ccf(ts$MRD, ts$MRS, main = "", ylab = "CCF", lag.max = 30)
title(main = "MRD & MRS")
ccf(ts$MRS, ts$MRD, main = "", ylab = "CCF", lag.max = 30)
title(main = "MRS & MRD")
acf(ts$MRS, main = "MRS", lag = 30)
par(mfrow = c(1,1))

combined = ts.intersect(ts(new_data$MRD, start= c(1969,10), frequency = 12), ts(new_data$MRS, start= c(1969,10), frequency = 12), dframe=TRUE)
spd = spec.pgram(combined, demean = TRUE, detrend = TRUE, spans = c(3,3))
dd = data.frame(freq=spd$freq, coh=spd$coh, per = 1/spd$freq)
order = dd[order(-dd$coh),]

plot(spd$freq, spd$coh, type = "l", xlab = "Frequency", ylab = "Coherence")
abline(v = order[1,1], col = "blue", lty = 2)

#### Extra Material ####

scale = 800
ut = paste("m",{supsc("3")}, "s", {supsc("-1")})
fd_data[,c("Date", "River Discharge", "Precipitation", "SMWE")] %>% mutate(SN = Precipitation + SMWE, month = format(Date, format = "%m")) %>%
  group_by(month) %>% summarise(mean = mean(`River Discharge`), mean_1 = mean(SN)*scale) %>% mutate(Month = month.name[as.numeric(month)]) %>%
  ggplot(aes(x = Month, group = 1)) + geom_col(aes(y = mean, fill = "Monthly Mean River Discharge"), width = 0.8) + geom_point(aes(y = mean_1, colour = "Mean Monthly Precipitation + Snowmelt"), size = 2.6) + geom_line(aes(y = mean_1, colour = "Mean Monthly Precipitation + Snowmelt"), size = 1.3) + scale_color_viridis(discrete = TRUE, name = "") +
  theme_bw() + labs(y = "", color = "Legend") + scale_y_continuous(label = unit_format(unit = ut), sec.axis = sec_axis(~./scale, name = "", label = unit_format(unit = "mm"))) +
  scale_fill_manual(name = "Legend:", values = c("Monthly Mean River Discharge" = "#e0a56e", "Mean Monthly Precipitation + Snowmelt" = "#0d6346"), limits = c("Monthly Mean River Discharge", "Mean Monthly Precipitation + Snowmelt")) + scale_color_manual(name = "Legend:", values = c("Monthly Mean River Discharge" = "#e0a56e", "Mean Monthly Precipitation + Snowmelt" = "#0d6346"), limits = c("Monthly Mean River Discharge","Mean Monthly Precipitation + Snowmelt")) +
  guides(color = guide_legend(override.aes = list(fill = c("#e0a56e", NA)))) + theme(legend.position = "bottom") + scale_x_discrete(limits = month.name)

P3_df = new_data[which(month(new_data[,2]) %in% c(4)),]
ggplot(mutate(new_data, vel = factor(ifelse(year(Date)>2000,"2000-2020","1969-1999")), Month = factor(month.abb[month(Date)],levels=month.abb), val = MRD), aes(x=Date, y=MRD, group=vel, color=vel)) + geom_line() + scale_color_viridis(discrete = TRUE, name="") + ggtitle("MRD") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
ggplot(mutate(P3_df, vel = factor(ifelse(year(Date)>2000,"2000-2020","1969-1999")), Month = factor(month.abb[month(Date)],levels=month.abb), val = MRD), aes(x=Date, y=MRD, group=Month, color=Month)) + geom_line() + scale_color_viridis(discrete = TRUE, name="") + ggtitle("MRD") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))


