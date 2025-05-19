spec <- read.csv('daily_SPEC_2021.csv')
sample <- spec[1:2000,]
length(unique(spec$Longitude))
#545 unique positions
sort(table(spec$Parameter.Name))
length(unique(spec$Date.Local))
#use Sulfate PM2.5 LC
spec_sulfate <- spec[spec$Parameter.Name == 'Sulfate PM2.5 LC',]
length(unique(spec_sulfate$Longitude))
#299 unique positions
spec_surfate_order <- spec_sulfate[order(spec_sulfate$Latitude),]
sort(table(spec_surfate_order$Longitude))
table(spec_sulfate$Date.Local)
write.csv(spec_surfate_order,'surfate2021.csv')
spec_sulfate_24 <- spec_surfate_order[spec_surfate_order$Sample.Duration == '24 HOUR',]
table(spec_sulfate_24$Date.Local)
spec_sulfate_24[spec_sulfate_24$Date.Local == '2021-12-26',]
table(spec_sulfate_24$Latitude)
write.csv(spec_sulfate_24,'sulfate2021_24.csv')
sort(table(spec_sulfate_24$Latitude))
spec_sulfate_24_3col <- spec_sulfate_24[,grepl("Latitude|Longitude|Arithmetic.Mean", colnames(spec_sulfate_24))] #提取列
write.csv(spec_sulfate_24_3col,'spec_sulfate_24_3col.csv')
library(fda)
fdata(spec_sulfate_24_3col[1:50,3], argvals = NULL, rangeval = NULL, names = NULL, fdata2d = FALSE)

#load and substract X - the characterisitics of counties
pop2021 <- read.csv('ACSST5Y2021.S0101-Data.csv',header = T, skip = 1)
pop2021_sub <- pop2021[,c(1,2,3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,127,131,
                          135,139,143)]
com2021 <- read.csv('ACSDT5Y2021.B28003-Data.csv',header = T, skip = 1)
com2021_sub <- com2021[,c(1,2,3,7,11,15,19,23)]
library(tidyverse)
spec_sulfate_24$name<-c(' County, ')
spec_sulfate_24$countyName<-str_c(spec_sulfate_24$County.Name,spec_sulfate_24$name)
spec_sulfate_24$fullname <- str_c(spec_sulfate_24$countyName,spec_sulfate_24$State.Name)
spec_sulfate_24$code <- str_c(spec_sulfate_24$State.Code,spec_sulfate_24$County.Code,spec_sulfate_24$Site.Num)
spec_sulfate_24$sitename <- str_c(spec_sulfate_24$fullname,spec_sulfate_24$Site.Num)
fullname <- data.frame(unique(spec_sulfate_24$fullname))
a <- merge(pop2021_sub, fullname, by.x = "Geographic.Area.Name", by.y = "unique.spec_sulfate_24.fullname.")
unique(spec_sulfate_24$sitename)
Y_f_299 <- read.csv('Y_f_299.csv',header = F)
Y_f_299$Latitude <- unique(spec_sulfate_24$Latitude)
Y_f_299$Longitude <- unique(spec_sulfate_24$Longitude)
Y_f_299$sitename <- unique(spec_sulfate_24$sitename)
Y_f_299$fullname <- gsub("\\d","",Y_f_299$sitename)
Y_f_264 <- Y_f_299[!duplicated(Y_f_299$fullname),]
XY <- merge(Y_f_264,pop2021_sub,by.x = 'fullname',by.y = "Geographic.Area.Name")
fulldata <- merge(XY,com2021_sub,by.x = 'fullname',by.y = "Geographic.Area.Name")
write.csv(fulldata,'fulldata.csv')

library(stringr)
aqidata <- read.csv('daily_aqi_by_county_2021.csv')
aqidata$fullname <- str_c(aqidata$State.Name,aqidata$county.Name)
uscounties <- read.csv('uscounties.csv')
uscounties$fullname <- str_c(uscounties$state_name,uscounties$county)
write.csv(aqidata,'aqidata.csv')
Y_f_1002 <- read.csv('Y_f_1002.csv',header = F)
Y_f_1002$fullname <- unique(aqidata$fullname)
fulldata <- merge(Y_f_1002,uscounties,by = 'fullname')
fulldata <- merge(fulldata,pop2021_sub,by.x = 'fullname',by.y = "Geographic.Area.Name")
fulldata$matchname <- str_c(fulldata$county_full,', ',fulldata$state_name)
XY <- merge(fulldata,com2021_sub,by.x = 'matchname',by.y = "Geographic.Area.Name")
XY <- merge(XY,pop2021_sub,by.x = 'matchname',by.y = "Geographic.Area.Name")
write.csv(XY,'fulldata_aqi.csv')

fulldata_aqi <- read.csv('fulldata_aqi.csv')
idx <- read.csv('idx_2.csv',header = F)
idx <- unlist(idx)
fulldata_aqi$idx <- idx
library(cluster)
library(factoextra)
fviz_nbclust(data.frame(fulldata_aqi$V90), kmeans, method = "wss")
kmeans_result <- kmeans(fulldata_aqi$V90, 5)
fulldata_aqi$idx <- kmeans_result$cluster



##extract origin data to plot
group <- unique(aqidata$fullname)
day1 <- data.frame(matrix(0,1002,11))
for (i in 1:length(group)){
  day1[i,] <- subset(aqidata,fullname == group[i] & Date == '2021-01-07')[1,]
}
day1 <- na.omit(day1)
day1full <- merge(day1,uscounties,by.x = 'X11', by.y = 'fullname')
library(cluster)
library(factoextra)
fviz_nbclust(data.frame(day1full$X6), kmeans, method = "wss")
kmeans_result <- kmeans(day1full$X6, 5)
day1full$idx <- kmeans_result$cluster


#install.packages('fadcluster')
library(fdacluster)
out1 <- fdakmeans(
  simulated30$x,
  simulated30$y,
  seeds = c(1, 21),
  n_clusters = 2,
  centroid_type = "mean",
  warping_class = "affine",
  metric = "pearson", 
  cluster_on_phase = FALSE
)
plot(out1, type = "amplitude")

b <- c(1:50)
B1x <- matrix(1:50,960,50,byrow = T)
B1y <- read.csv('B1.csv',header = F)
B1y <- as.matrix(B1y)
B1y <- array(B1y, c(960, 1, 50))
out2 <- fdakmeans(
  B1x,
  B1y,
  seeds = c(1,11,21,31),
  n_clusters = 4,
  centroid_type = "mean",
  warping_class = "srsf",
  metric = "pearson", 
  cluster_on_phase = FALSE
)
center <- matrix(out2$center_curves,4,50)
write.csv(center,'center.csv')
