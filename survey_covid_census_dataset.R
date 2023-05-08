#Load Data
AR <- read.csv("https://www.dropbox.com/s/rua5k2havwh7qp5/AR.csv?dl=1")
CT <- read.csv("https://www.dropbox.com/s/wcuop8ukry8qxgs/CT.csv?dl=1")
FL <- read.csv("https://www.dropbox.com/s/m5z3y0p8atu3ml2/FL.csv?dl=1")
IA <- read.csv("https://www.dropbox.com/s/8bdrsg29fgkfvfk/IA.csv?dl=1")
MO <- read.csv("https://www.dropbox.com/s/e4cvme3s3f0853w/MO.csv?dl=1")
NC <- read.csv("https://www.dropbox.com/s/q8syfrkwxjuz793/NC.csv?dl=1")
SC <- read.csv("https://www.dropbox.com/s/2hahe4o8x4rm3ya/SC.csv?dl=1")
TX <- read.csv("https://www.dropbox.com/s/5w39v9d8qvbr9yf/TX.csv?dl=1")
UT <- read.csv("https://www.dropbox.com/s/v606rrufsspzkmh/UT.csv?dl=1")
WY <- read.csv("https://www.dropbox.com/s/gcg3ztqr1is7n8i/WY.csv?dl=1")

####Linking survey data to shapefiles####
library(tigris)
library(tidyverse)
library(sf)

# Multiple States (first half)
df1 <- rbind(AR, CT, FL, IA, MO)
df1 <- filter(df1, Finished == TRUE & Q4 == "Yes, I give my consent to participate in this study.")
df1.pts <- st_as_sf(df1, coords=c("LocationLongitude", "LocationLatitude"))
df1_count <- tigris::counties(state = c("Arkansas", "Connecticut", "Florida", "Iowa", "Missouri"))
st_crs(df1.pts) <- 4269
pts_4_df1 <- st_intersection(df1.pts, df1_count)
head(pts_4_df1)
ggplot() +
  geom_sf(data=df1_count) +
  geom_sf(data=pts_4_df1) 

# Multiple States (second half)
df2 <- rbind(NC, SC, TX, UT, WY)
df2 <- filter(df2, Finished == TRUE & Q4 == "Yes, I give my consent to participate in this study.")
df2.pts <- st_as_sf(df2, coords=c("LocationLongitude", "LocationLatitude"))
df2_count <- tigris::counties(state = c("North Carolina", "South Carolina", "Texas", "Utah", "Wyoming"))
st_crs(df2.pts) <- 4269
pts_4_df2 <- st_intersection(df2.pts, df2_count)
head(pts_4_df2)
ggplot() +
  geom_sf(data=df2_count) +
  geom_sf(data=pts_4_df2) 

#All states
pts_4_all <- rbind(pts_4_df1, pts_4_df2)
all_count <- rbind(df1_count, df2_count)
ggplot() +
  geom_sf(data=all_count) +
  geom_sf(data=pts_4_all) 

#Make unique fips identifier
pts_4_all$fips <- paste(pts_4_all$STATEFP, pts_4_all$COUNTYFP)
pts_4_all$fips <- gsub(" ", "", paste(pts_4_all$STATEFP, pts_4_all$COUNTYF))

####Linking survey shapefiles to covid data####
#CovidDat=read.csv("COVID-19_Case_Surveillance_Public_Use_Data_with_Geography (3).csv",header=TRUE)
#CDC data seems too large for both R and Excel so using nyt data for now

library (readr)
urlfile="https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"
CovidDat<-read_csv(url(urlfile))

#select survey states
CovidDat <- CovidDat %>% filter(state %in% c("Arkansas", "Connecticut", "Florida", "Iowa", "Missouri", "North Carolina", "South Carolina", "Texas", "Utah", "Wyoming"))

#discard dates before cases were tracked and after survey scope.
CovidSmoothed <- CovidDat %>% filter(date > as.Date("2020-03-01") & date < as.Date("2020-05-31"))

#Aggregate total cases by county
CaseTotals=aggregate(list(CovidSmoothed$cases),
                     by=list(CovidSmoothed$fips),
                     FUN=max)
colnames(CaseTotals)<-c("fips","cases")

#pad fips code for uniform match between datasets
library(stringr)
CaseTotals$fips <- str_pad(CaseTotals$fips, 5, pad = "0")

#Merge data
mergedDat=merge(pts_4_all, CaseTotals, by="fips")

####Including census data####
#select survey states
census2020 <- read.csv("https://www.dropbox.com/s/3gqydymogj1ts08/CountyCensus2020.csv?dl=1")

census2020 <- census2020 %>%
  filter(STNAME %in% c("Arkansas", "Connecticut", "Florida", "Iowa", "Missouri", "North Carolina", "South Carolina", "Texas", "Utah", "Wyoming"))


#Make unique fips identifier
census2020$STATE2 <- str_pad(census2020$STATE, 2, pad = "0")
census2020$COUNTY2 <- str_pad(census2020$COUNTY, 3, pad = "0")
census2020$fips <- paste(census2020$STATE2, census2020$COUNTY2)
census2020$fips <- gsub(" ", "", paste(census2020$STATE2, census2020$COUNTY2))

#Merge data
allDat=merge(mergedDat, census2020[,c("fips", "POPESTIMATE2020")], by="fips")

#I stuck with merge because it was a better indicator if something went wrong, but I can switch to left join if need be

####Calculate population density, add categorical and categorical groupings####
allDat$area <- allDat$ALAND + allDat$AWATER
allDat$area <- allDat$area / 1000000 #convert to square km
allDat$POPdensity <- allDat$POPESTIMATE2020 / allDat$area

#categorically define urban and rural counties
allDat$URcategory <- cut(allDat$POPdensity, breaks=c(0,500,Inf), include.lowest=TRUE, labels=c("rural","urban"))
summary(allDat)
head(allDat)
####Calculate covid density####
allDat$COVIDdensity <- allDat$cases / allDat$POPESTIMATE2020

#save merged data 
setwd("C:/Users/akaz1/Dropbox/COVID Project/Manuscripts/UrbanRural Analysis/Analysis")
write.csv(allDat,'allDat.csv')
