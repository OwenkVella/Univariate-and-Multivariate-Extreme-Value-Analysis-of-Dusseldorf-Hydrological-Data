#Required Packages
#install.packages(c("tcltk", "dplyr", "writexl"))
lapply(c("tcltk", "dplyr", "writexl"), require, character.only = TRUE)

#To clean
rm(list = ls())
graphics.off()

#Loading and formatting the data
data = read.csv(tk_choose.files(default = "", caption = "Select csv file"))
names(data)[1] = "Date"
data[,"Date"] = as.Date(data[,1], format = "%d/%m/%Y")
new_data = data[1:(which(is.na(data$Mean_Daily_Discharge))[1]-1),]
new_data[,"Water_Level"] = new_data[,"Water_Level"]*0.01

#Scatter plot of Discharge vs Water Level
plot(new_data$Water_Level, new_data$Mean_Daily_Discharge, xlab = "Water Level", ylab = "River Discharge")
points(new_data$Water_Level[which(format(new_data$Date, format = "%Y") == "2000")[1]:length(new_data$Water_Level)], 
       new_data$Mean_Daily_Discharge[which(format(new_data$Date, format = "%Y") == "2000")[1]:length(new_data$Mean_Daily_Discharge)], col = "red")

#Factoring the data
factory = 20
datas = mutate(new_data, Year = as.numeric(format(Date, format = "%Y")), Duration = (Year - min(Year)), Group = ifelse(Duration%%factory == 0, (Duration%/%(factory)),(Duration%/%(factory) + 1)))
datas[which(datas$Year == min(datas$Year)),"Group"] = 1
datas$Group = factor(datas$Group)
datas = mutate(datas, Legend = "1")
for(i in 1:length(datas$Group)){
  
  fili = datas$Year[which(datas$Group == datas$Group[i])]
  
  datas[i,"Legend"] = paste0("[", min(fili),", ", max(fili), "]")
  
}
datas$Legend = factor(datas$Legend)

#Changing the palette and plotting
palette("default")
palette(c("magenta3", palette()[-which(palette() == "black")]))
palette(c("gold2", palette()[-which(palette() == "#28E2E5")]))
palette(c("royalblue", palette()[-which(palette() == "#2297E6")]))

plot(datas$Water_Level, datas$Mean_Daily_Discharge, xlab = "Water Level (m)", ylab = "", col = datas$Legend)
legend(x = "topleft",levels(datas$Legend),col=1:length(levels(datas$Legend)),pch=1)
title(ylab = bquote("River Discharge (" ~ m^3 ~ s^-1 ~ ")"), line = 2.5)

#Plot the relationship of Water level vs Discharge for the last group
new_datas = datas[which(datas$Group == tail(levels(datas$Group),1)),]
new_datas$Water_Level = new_datas$Water_Level/0.01

model = stats::nls(Mean_Daily_Discharge~C*(Water_Level + a)^(B), data=new_datas, start = list(C = 1, a = 1, B = 1))
summary(model)
nlstools::confint2(model, level = 0.95, method = "asymptotic")
cor(new_datas[,3], predict(model))

new_datas$Water_Level = new_datas$Water_Level*0.01
plot(new_datas[,c(2,3)], xlab = "Water Level (m)", ylab = "")
lines(new_datas[,2], predict(model), col = "red")
title(ylab = bquote("River Discharge (" ~ m^3 ~ s^-1 ~ ")"), line = 2.5)

#Extract river discharge data for the period 2020-2021
finall_data = data
for(j in 1:length(which(is.na(data$Mean_Daily_Discharge)))){
  
  finall_data[which(is.na(data$Mean_Daily_Discharge))[j], "Mean_Daily_Discharge"] = as.numeric(coef(model)[1])*(finall_data[which(is.na(data$Mean_Daily_Discharge))[j], "Water_Level"] + as.numeric(coef(model)[2]))^(as.numeric(coef(model)[3]))
  
}
finall_data$Water_Level = finall_data$Water_Level*0.01
plot(finall_data$Water_Level, finall_data$Mean_Daily_Discharge, xlab = "Water Level (m)", ylab = "")
points(finall_data[which(is.na(data$Mean_Daily_Discharge)),"Water_Level"], finall_data[which(is.na(data$Mean_Daily_Discharge)),"Mean_Daily_Discharge"], col = "red")
title(ylab = bquote("River Discharge (" ~ m^3 ~ s^-1 ~ ")"), line = 2.5)

#To save the data
#write_xlsx(finall_data, paths)
paths = paste0(gsub("/","\\\\", tk_choose.dir(default = "", caption = "Select directory")),"\\Extended_DR_Dataset.csv")
write.csv(finall_data, paths, row.names=FALSE)

