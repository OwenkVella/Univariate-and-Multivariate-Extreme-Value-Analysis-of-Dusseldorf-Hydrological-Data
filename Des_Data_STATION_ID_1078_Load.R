#Required Packages
#install.packages(c("tcltk", "lubridate", "dplyr", "ggplot2", "viridis", "scales"))
lapply(c("tcltk", "lubridate", "dplyr", "ggplot2", "viridis", "scales"), require, character.only = TRUE)

#To clean
rm(list = ls())
graphics.off()

#Loading data
#data = read.csv(file.choose())
data = read.csv(tk_choose.files(default = "", caption = "Select csv file"))

#------------------ Done ~ Formatting Data -----------------

#To format the data
data[,2] = as.Date(data[,2], format = "%d/%m/%Y")
names(data) = gsub("."," ", names(data), fixed = TRUE)
data = mutate(data, Month = factor(month.abb[month(Date)], levels=month.abb), Year = as.numeric(format(data[,2], format = "%Y")), Key = paste0(Month, "_", Year))
data = dplyr::select(data, 1:2, (length(names(data)) - 2):(length(names(data))), 3:(length(names(data)) - 2))

#------------------ Done ~ Plotting Data / Providing Stats  -----------------

#Graph parameters for daily mean
scale = 6
scale_1 = 15
names(data)[7] = "value"
names(data)[8] = "value_1"
names(data)[6] = "value_2"
names(data)[9] = "value_3"
names(data)[10] = "value_4"

#Daily Amount of Precipitation (mm) vs Daily Mean Temperature (intToUtf8(176) C)
ns = paste0("Daily Mean Temperature (",intToUtf8(176),"C)")
data %>%
  group_by(Month) %>%
  summarise(mean = mean(value), mean_1 = mean(value_1)/scale) %>%
  ggplot(aes(x = Month, group = 1)) + geom_col(aes(y = mean, fill = "Daily Amount of Precipitation (mm)"), width = 0.8) + geom_point(aes(y = mean_1, colour = ns), size = 2.6) + geom_line(aes(y = mean_1, colour = ns), size = 1.3) + scale_color_viridis(discrete = TRUE, name = "") +
  ggtitle(paste0("Climate chart -  D",intToUtf8(252), "sseldorf")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + labs(y = "", color = "Legend") + scale_y_continuous(label = unit_format(unit = "mm"), breaks = seq(0, 3, 0.5), sec.axis = sec_axis(~.*scale, name = "", label = unit_format(unit = paste0(intToUtf8(176), "C")), breaks = seq(0, 30, 3))) +
  scale_fill_manual(name = "Legend:", values = c("Daily Amount of Precipitation (mm)" = "#0d6346", ns = "#d2290f"), limits = c("Daily Amount of Precipitation (mm)", ns)) + scale_color_manual(name = "Legend:", values = c("Daily Amount of Precipitation (mm)" = "#0d6346",ns = "#d2290f"), limits = c("Daily Amount of Precipitation (mm)",ns)) +
  guides(color = guide_legend(override.aes = list(fill = c("#0d6346", NA)))) + theme(legend.position = "bottom")

#Sunshine Duration Daily Total (h) vs Daily Mean of the Relative Humidity (%) vs Daily mean wind speed (m/s)
# data %>%
#   group_by(Month) %>%
#   summarise(mean = mean(na.omit(value_2)), mean_1 = mean(na.omit(value_3))/scale_1) %>%
#   ggplot(aes(x = Month, group = 1)) + geom_col(aes(y = mean, fill = "Sunshine Duration Daily Total (h)"), width = 0.8) + geom_text(aes(y = round(mean,2), label = round(mean,2)), position=position_dodge(width = 0.9), vjust = 1.5, colour = "white") + geom_point(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 2.6) + geom_line(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 1.3) + scale_color_viridis(discrete = TRUE, name = "") +
#   ggtitle(paste0("Climate chart -  D",intToUtf8(252), "sseldorf")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + labs(y = "", color = "Legend") + scale_y_continuous(label = unit_format(unit = "h"), breaks = seq(0, 16, 0.5), sec.axis = sec_axis(~.*scale_1, name = "", label = unit_format(unit = "%"), breaks = seq(0, 105, 8))) +
#   scale_fill_manual(name = "Legend:", values = c("Sunshine Duration Daily Total (h)" = "#2f09ad", "Daily Mean of the Relative Humidity (%)" = "#b80d0d"), limits = c("Sunshine Duration Daily Total (h)", "Daily Mean of the Relative Humidity (%)")) + scale_color_manual(name = "Legend:", values = c("Sunshine Duration Daily Total (h)" = "#2f09ad", "Daily Mean of the Relative Humidity (%)" = "#b80d0d"), limits = c("Sunshine Duration Daily Total (h)", "Daily Mean of the Relative Humidity (%)")) +
#   guides(color = guide_legend(override.aes = list(fill = c("#2f09ad", NA)))) + theme(legend.position = "bottom")
data %>%
  group_by(Month) %>%
  summarise(mean = mean(na.omit(value_2)), mean_1 = mean(na.omit(value_3))/scale_1, mean_2 = mean(na.omit(value_4))) %>%
  ggplot(aes(x = Month, group = 1)) + geom_point(aes(y = mean_2, colour = "Daily Mean Wind Speed (m/s)"), size = 2.6) + geom_line(aes(y = mean_2, colour = "Daily Mean Wind Speed (m/s)"), size = 1.3) + 
  geom_point(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 2.6) + geom_line(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 1.3) +
  geom_col(aes(y = mean, fill = "Sunshine Duration Daily Total (h)"), width = 0.8, position = position_dodge(2)) + geom_text(aes(y = round(mean,2), label =paste(round(mean,2), "h")), position=position_dodge(width = 0.9), vjust = 1.5, colour = "white") +
  geom_point(aes(y = mean_2, colour = "Daily Mean Wind Speed (m/s)"), size = 2.6) + geom_line(aes(y = mean_2, colour = "Daily Mean Wind Speed (m/s)"), size = 1.3) + 
  geom_point(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 2.6) + geom_line(aes(y = mean_1, colour = "Daily Mean of the Relative Humidity (%)"), size = 1.3) +
  ggtitle(paste0("Climate chart -  D",intToUtf8(252), "sseldorf")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + labs(y = "", color = "Legend: ") + theme(legend.position = "bottom") + 
  scale_y_continuous(label = unit_format(unit = "m/s"), breaks = seq(0, 8, 0.5), sec.axis = sec_axis(~.*scale_1, name = "", label = unit_format(unit = "%"), breaks = seq(0, 100, 9))) +
  scale_fill_manual(name = "Legend:", values = c("Sunshine Duration Daily Total (h)" = "#2f09ad", "Daily Mean of the Relative Humidity (%)" = "#b80d0d", "Daily Mean Wind Speed (m/s)" = "#1da4a8"), limits = c("Sunshine Duration Daily Total (h)", "Daily Mean of the Relative Humidity (%)", "Daily Mean Wind Speed (m/s)")) + scale_color_manual(name = "Legend:", values = c("Sunshine Duration Daily Total (h)" = "#2f09ad", "Daily Mean of the Relative Humidity (%)" = "#b80d0d", "Daily Mean Wind Speed (m/s)" = "#1da4a8"), limits = c("Sunshine Duration Daily Total (h)", "Daily Mean of the Relative Humidity (%)", "Daily Mean Wind Speed (m/s)")) +
  guides(color = guide_legend(override.aes = list(fill = c("#2f09ad", NA, NA))))
  
