library(echarts4r)
library(colorspace)
#XY <- read.csv('fulldata_aqi.csv')
XY <- day1full
XY$new <- XY$X6
XY$new <- XY$idx
XY |>
  e_charts(x = lng) |>
  e_geo(
    roam = TRUE,
    map = "world",
    boundingCoords = list(
      c(-170, 0),
      c(0, 70)
    )
  ) |>
  e_scatter(
    serie = lat,
    scale = NULL, # 去掉尺度变换
    size = new, 
    symbol = 'circle',
    coord_system = "geo"
  ) |>
  
  e_visual_map(
    serie = new, 
    type = 'piecewise', 
    left = 0, bottom = 0, 
    text  = c("high", "low"), 
    textStyle = list(color = "white"),
    inRange = list(
      #color = c('lightgreen','yellow', 'red','purple'), 
      #color = c('sienna','palegreen4','cadetblue'), 
      #color = hcl.colors(3,palette=hcl.pals()[15]), 
      #color = c("#1e5670","#4198b9","#6bb3c0","#91cfc9","#cde8f3"),
      #color = rainbow_hcl(5),
      color = hcl.colors(3,palette="Set2"),
      colorAlpha = 1, 
      symbolSize = c(12, 12) 
      
    )
  ) |>
  e_theme(name = "chalk")



fulldata_aqi$new <- fulldata_aqi$idx
# 绘图
fulldata_aqi |>
  e_charts(x = lng) |>
  e_geo(
    roam = TRUE,
    map = "world",
    boundingCoords = list(
      c(-170, 0),
      c(0, 70)
    )
  ) |>
  e_scatter(
    serie = lat,
    scale = NULL,
    size = new, 
    symbol = 'pin',
    coord_system = "geo"
  ) |>
  
  e_visual_map(
    serie = new, 
    left = 0, bottom = 0, 
    text  = c("high", "low"),
    textStyle = list(color = "white"),
    inRange = list(
      color = hcl.colors(5,palette=hcl.pals()[15]), 
      colorAlpha = 1, 
      symbolSize = c(9, 9) 
      
    )
  ) |>
  e_theme(name = "chalk")
 
library(plotly)
plotly::plot_ly(
  data = day1full,
  type = "choropleth",
  locations = ~county_fips,
  locationmode = "USA-states", 
  colorscale = "Viridis", 
  z = ~X6
) |>
  plotly::colorbar(title = "population") |> 
  plotly::layout(
    geo = list(scope = "usa"),
    title = "Population in 1974"
  ) |>
  plotly::config(displayModeBar = FALSE)

library(usmap)
day1us <- merge(aqidata,uscounties,by='fullname')
day1us <- data.frame(cbind(day1us$county_fips,day1us$AQI))
colnames(day1us) <- c('fips','values')
plot_usmap(data = day1us)
group <- unique(aqidata$fullname)
day1 <- data.frame(matrix(0,2000,11))
for (i in 1:length(group)){
  day1[i,] <- subset(aqidata,fullname == group[i] & Date == '2021-01-01')[1,]
}
day1 <- na.omit(day1)
day1full <- merge(day1,uscounties,by.x = 'X11', by.y = 'fullname')




plot(day1$X6)
plot(aqidata[1:280,6])
library(xts)
aqi2 <- aqidata[aqidata$State.Code==20&aqidata$County.Code==57,]
aqits1 <- xts(aqidata[1:280,6], order.by=as.Date(aqidata[1:280,5], "%Y-%m-%d"))
aqits2 <- xts(aqi2[,6], order.by=as.Date(aqi2[,5], "%Y-%m-%d"))
plot(aqits2)
plot(aqits1, main = "Dividend date and amount")
lines(aqits2, col = "orange", lwd = 2)
axis(side = 4, at = pretty(data$citigroup), col = "orange")



library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
incidents <-read.table("https://data.seattle.gov/api/views/3k2p-39jp/rows.csv?accessType=DOWNLOAD", head=TRUE, sep=",", fill=TRUE, stringsAsFactors=F) 

#If the above data set is unavailable please use this code
df= fread('https://raw.githubusercontent.com/lgellis/MiscTutorial/master/ggmap/i2Sample.csv', stringsAsFactors = FALSE)
#incidents <- df
incidents <- read.csv('i2Sample.csv')

#Assign color variables
col1 = "#d8e1cf" 
col2 = "#438484"

#Peek at the data set and attach the column names
head(incidents)
attach(incidents)
str(incidents)
incidents$ymd <-mdy_hms(Event.Clearance.Date)
incidents$month <- month(incidents$ymd)
incidents$year <- year(incidents$ymd)
incidents$wday <- wday(incidents$ymd)
incidents$hour <- hour(incidents$ymd)
attach(incidents)
head(incidents)
groupSummary <- ddply(incidents, c( "Event.Clearance.Group", "hour"), summarise,
                      N    = length(ymd)
)
groupSummary <- ddply(aqidata[1:50000,], c( "fullname", "Date"), summarise,
                      N    = length(aqidata)
)
groupSummarysub <- groupSummary[1:10000,]
groupsub = sample(group,20)
groupsub1 <- groupsub[c(1,2,4,5,7,8,10,11,17,18,20)]
groupsub1 <- groupsub[c(18,9,1,7,13,6,11,8,3,15)]
aqidatasub = aqidata[aqidata$fullname %in% groupsub1,]
aqidatasub$Date <- as.Date(aqidatasub$Date)
#overall summary
ggplot(aqidatasub, aes(Date,fullname)) + geom_tile(aes(fill = AQI),colour = "white") +
  scale_fill_gradient2(limits=c(0, 320),mid = 'chartreuse',high = 'firebrick1',midpoint = 0)
  #scale_fill_gradient(low = 'springgreen4',high = 'firebrick') +  
  guides(fill=guide_legend(title="AQI Index")) +
  labs(x = "Date", y = "County") +
  theme_bw() + theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()+
          scale_x_date(date_breaks = "1 month"))
 
Sys.setlocale('LC_TIME','English')
