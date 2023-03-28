#Required Packages
#install.packages(c("grid", "tcltk", "ncdf4", "raster", "rgdal", "stringr", "lubridate", "dplyr", "ggplot2", "viridis", "scales"))
lapply(c("grid", "tcltk", "ncdf4", "raster", "rgdal", "stringr", "lubridate", "dplyr", "ggplot2", "viridis", "scales"), require, character.only = TRUE)

#To clean
rm(list = ls())
graphics.off()

#### NetCDF format ####
setwd(tk_choose.dir(default = "", caption = "Select directory"))
name = file.choose()

nc_data = nc_open(name)

sink('FLUX.txt')
print(nc_data)
sink()

#Workings Not important
#nc_data[["dim"]][["time"]][["vals"]]
#nc_data[["dim"]][["time"]][["units"]]
#lat = ncatt_get(nc_data, attributes(nc_data$var)$names[1])
#lon = ncatt_get(nc_data, attributes(nc_data$var)$names[2])
#t = ncvar_get(nc_data, "time")
#lon = ncvar_get(nc_data, "lon")
#lat = ncvar_get(nc_data, "lat")
#pet = ncatt_get(nc_data, attributes(nc_data$var)$names[3])

#Relative dimensions
t = as.Date(ncvar_get(nc_data, attributes(nc_data$dim)$names[1]), origin = as.Date(trimws(str_remove_all(nc_data[["dim"]][["time"]][["units"]], "[a-z]"))))
lat = ncvar_get(nc_data, attributes(nc_data$var)$names[1])
lon = ncvar_get(nc_data, attributes(nc_data$var)$names[2])

#Selected coordinates in the vicinity of Station D?sseldorf:
#Station D?sseldorf: lat = 51.2960 and lon = 6.7686
sel_x = 25 
sel_y = 106
x_index = as.numeric(which(lat == lat[sel_x,sel_y]))
y_index = as.numeric(which(lon == lon[sel_x,sel_y]))

pet = ncvar_get(nc_data, attributes(nc_data$var)$names[3], start = c(sel_x, sel_y, 1), count = c(length(x_index), length(y_index), -1))
data = data.frame("Date" = t, "Evapotranspiration" = as.vector(pet))
write.csv(data, "PET.csv", row.names = F)

#Graphical plots
pt = data[which(format(data[,1], format = "%Y") == 1980),]
plot(pt, type = "l", xlab = "Date", ylab = "Observations", main = paste(str_to_title(paste0(nc_data[["var"]][["pet"]][["longname"]])), "in", format(pt[1,1], format = "%Y")))
plot(t, pet, type = "l", xlab = "Date", ylab = "Observations", main = str_to_title(paste0(nc_data[["var"]][["pet"]][["longname"]])))

#Average Visualisation
data %>% mutate(Month = factor(month.abb[month(Date)], levels=month.abb)) %>%
  group_by(Month) %>%
  summarise(mean = mean(Evapotranspiration)) %>%
  ggplot(aes(x = Month, y = mean, group = 1)) + geom_point(aes(y = mean, colour = "Daily Potential Evapotranspiration (mm)"), size = 2.6) + geom_line(aes(y = mean, colour = "Daily Potential Evapotranspiration (mm)"), size = 1.3) + scale_color_viridis(discrete = TRUE, name = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + labs(y = "") + scale_y_continuous(label = unit_format(unit = "mm"), breaks = seq(0, 6, 0.4)) + theme(legend.position = "bottom") +
  scale_fill_manual(name = "Legend:", values = c("Daily Potential Evapotranspiration (mm)" = "#d2290f"), limits = c("Daily Potential Evapotranspiration (mm)")) + scale_color_manual(name = "Legend:", values = c("Daily Potential Evapotranspiration (mm)" = "#d2290f"), limits = c("Daily Potential Evapotranspiration (mm)")) +
  guides(color = guide_legend(override.aes = list(fill = c(NA)))) + theme(legend.position = "bottom")

