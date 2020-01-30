devtools::install_github("USGS-R/EflowStats")

#need to build site x environment dataframe
head(SiteID) #these are my sites
#want to find oldest and youngest year at each site
Site_data<-AxL %>% group_by(Species, Site) %>%
  summarize(maxYear=max(Year),
            minYear=min(Year)) %>% 
  left_join(SiteID) %>%
  dplyr::select()

library(dataRetrieval)

#use Eflow to pull out some hydrological variables -----
library(EflowStats)
dailyQ <- readNWISdata(siteNumber=SiteID$`USGS gage Identifier`[1],
                       parameterCd="00060",
                       startDate="",
                       endDate="")
dailyQClean <- validate_data(dailyQ[c("Date","X_00060_00003")],yearType="water")



#Get drainage area
siteInfo <- readNWISsite(siteNumber = "04085427")
drainageArea <- siteInfo$drain_area_va

#Get peak flows
peakFlows <- readNWISpeak(siteNumber = "04085427",
                          startDate = "2000-10-01",
                          endDate = "2012-9-30")

#Check data for completeness

dailyQClean <- validate_data(dailyQ[c("Date", "X_00060_00003")], yearType="water")

#Get flood recurence threshold
floodThresh <- get_peakThreshold(dailyQClean[c("date","discharge")],
                                peakFlows[c("peak_dt","peak_va")])

#Calculate all hit stats
calc_allHITOut <- calc_allHIT(dailyQClean,
                             drainArea=drainageArea,
                              floodThreshold=floodThresh)

# figure out water temperature -----



install.packages("nhdplusTools")
library(nhdplusTools)
library(sf)

start_point <- st_sfc(st_point(c(-89.362239, 43.090266)), crs = 4269)
start_comid <- discover_nhdplus_id(start_point)

flowline <- navigate_nldi(list(featureSource = "comid", 
                               featureID = start_comid), 
                          mode = "upstreamTributaries", 
                          data_source = "")

subset_file <- tempfile(fileext = ".gpkg")
subset <- subset_nhdplus(comids = flowline$nhdplus_comid,
                         output_file = subset_file,
                         nhdplus_data = "download", 
                         return_data = TRUE)
flowline <- subset$NHDFlowline_Network
catchment <- subset$CatchmentSP
waterbody <- subset$NHDWaterbody


#
