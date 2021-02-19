library(lubridate);library(dataRetrieval);library(EflowStats)

#need to build site x environment dataframe
head(SiteID) #these are my sites
#want to find oldest and youngest year at each site
Site.data<-AxL %>% group_by(Species, Site) %>%
  rename(SiteID=Site)%>%
  summarize(maxYear=max(Year),
            minYear=min(Year)) %>% 
  left_join(SiteID, by="SiteID") %>%
  dplyr::select(Species,Site.Agg, HUC.8,HUC.12,USGS.Gage, maxYear, minYear, Lat.cor, Long.cor) %>%
  filter(!duplicated(Lat.cor))%>%
  rowwise() %>%
  mutate(d.av.begin=min(whatNWISdata(siteNumber=USGS.Gage, parameterCd="00060", service="dv")$begin_date),
         prob.year=minYear-year(d.av.begin)) %>%
  arrange(prob.year)
#whatNWISdata(siteNumber="05586100", parameterCd="00060", service="dv")

#use Eflow to pull out some hydrological variables -----
get_dis_stats<- function(site,species, gage, maxyear, shellyear, gageyear, hit_var) {
  dailyQ<-readNWISdata(siteNumber=gage,
                      parameterCd="00060",
                      startDate=paste(max(shellyear, gageyear+1),"10-01", sep="-"),
                      endDate=paste(maxyear,"09-30", sep="-"),
                      service="dv") %>%
  mutate(Date=as.Date(dateTime))
  dailyQClean <- validate_data(dailyQ[c("Date","X_00060_00003")],
                               yearType="water") 
  siteInfo <- readNWISsite(siteNumber = gage)
  drainageArea <- siteInfo$drain_area_va #get drainage area

  #Get peak flows
  peakFlows <- readNWISpeak(siteNumber = gage,
                            startDate=paste(max(shellyear, gageyear+1),"10-01", sep="-"),
                            endDate=paste(maxyear,"09-30", sep="-"))
  #Get flood recurence threshold
  floodThresh <- get_peakThreshold(dailyQClean[c("date","discharge")],
                                  peakFlows[c("peak_dt","peak_va")])
  #Calculate all hit stats
  #get a lot of warnings, but its just a slight function language change,
  #functionality is the same (changed . to _)
  calc_allHITOut <- calc_allHIT(dailyQClean,
                             drainArea=drainageArea,
                              floodThreshold=floodThresh) %>%
  filter(indice %in% hit_var) %>%
  mutate(Site.Agg=site,
         Sp=species)
  return(calc_allHITOut)
}

#get_dis_stats(site="test",gage="05341550",maxyear=2017,shellyear=1986,gageyear=2011)
site.dis<-NULL
for(r in c(1:3,5:14,17:30,32:40)){
  dis.t<-get_dis_stats(site=Site.data$Site.Agg[r],
                       species=Site.data$Species[r],
                       gage=Site.data$USGS.Gage[r],
                       maxyear=Site.data$maxYear[r],
                       shellyear=Site.data$minYear[r],
                       gageyear=year(Site.data$d.av.begin[r]),
                       hit_var=c("ma2","ma3","ma36","ml16", "ml19",
                                 "ml5","ml6","ml7","ml8",
                                 "mh5","mh6","mh7","mh8",
                                 "fl1","fh1","dl16","dh15",
                                 "ra8","ra9")) 
  site.dis<-bind_rows(site.dis, dis.t)
}

#warnings due to: outdated function notation and discharge = 0 (its replaced with 0.1)

#problem gages with missing data;
#will manipulate data and get HIT values
Site.data[4,] #doesn't work because gage problem in 1991; removed 1991, recalculated 
Site.data[15,] #missing 1/24/2017 womp womp; didn't include 2017
Site.data[16,] #missing one date in 2017; didn't include 2017
Site.data[31,] #needed a full year of data; included 2014 in data

hit_var1<-c("ma2","ma3","ma36","ml16", "ml19","ml5","ml6","ml7","ml8",
          "mh5","mh6","mh7","mh8","fl1","fh1","dl16","dh15","ra8","ra9")

dis.FC.post91<-get_dis_stats(site=Site.data$Site.Agg[4], species=Site.data$Species[4],
                             gage=Site.data$USGS.Gage[4],
                     maxyear=Site.data$maxYear[4], shellyear=1991,
                     gageyear=year(Site.data$d.av.begin[4]), hit_var = hit_var1) 
dis.FC.pre91<-get_dis_stats(site=Site.data$Site.Agg[4], species=Site.data$Species[4],
                            gage=Site.data$USGS.Gage[4],
                             maxyear=1990, shellyear=Site.data$minYear[4],
                             gageyear=year(Site.data$d.av.begin[4]), hit_var = hit_var1)
dis.FC<-bind_rows(dis.FC.post91, dis.FC.pre91) %>%
  group_by(indice, Site.Agg,Sp) %>% summarize(statistic=mean(statistic))
dis.Flor.no17<-get_dis_stats(site=Site.data$Site.Agg[15], species=Site.data$Species[15],
                             gage=Site.data$USGS.Gage[15],
                            maxyear=2016, shellyear=Site.data$minYear[15],
                            gageyear=year(Site.data$d.av.begin[15]), hit_var = hit_var1) 

dis.Kish<-get_dis_stats(site=Site.data$Site.Agg[16], species=Site.data$Species[16],
              gage=Site.data$USGS.Gage[16],
              maxyear=2016, shellyear=Site.data$minYear[16],
              gageyear=year(Site.data$d.av.begin[16]), hit_var = hit_var1) 

dis.W2.add2014<-get_dis_stats(site=Site.data$Site.Agg[31], species=Site.data$Species[31],
                              gage=Site.data$USGS.Gage[31],
                            maxyear=Site.data$maxYear[31], shellyear=2014,
                            gageyear=year(Site.data$d.av.begin[31]), hit_var = hit_var1) 
site.dis.com<-bind_rows(site.dis,dis.FC, dis.Flor.no17, dis.Kish, dis.W2.add2014)

site.hyd<-site.dis.com %>%
  group_by(Site.Agg,Sp) %>%
  pivot_wider(names_from=indice, values_from=statistic)

write.csv(site.hyd, "data/HITDischarge.csv")

# figure out stream characteristics -----
library(nhdplusTools); library(sf)

work_dir <- tempdir(check = TRUE)
sample_data <- system.file("extdata/sample_natseamless.gpkg", package = "nhdplusTools")

hr_gpkg <- file.path(work_dir, "hr_data.gpkg")


hr_data_dir <- download_nhdplushr(work_dir, "1114")# Download some NHDPlusHR Data

hr <- get_nhdplushr(hr_data_dir)
sapply(hr, class)
sapply(hr, nrow)
str(hr)

#MY CODE
nhdata.all<-NULL
for(s in 1:nrow(Site.data)){
  start_point <- st_sfc(st_point(c(Site.data$Long.cor[s], 
                                   Site.data$Lat.cor[s])), crs = 4269)
  start_comid <- discover_nhdplus_id(start_point)
  
  subset_file <- tempfile(fileext = ".gpkg")
  subset <- subset_nhdplus(comids = start_comid,
                           output_file = subset_file,
                           nhdplus_data = "download", 
                           return_data = TRUE)
  
  nhdata <- subset$NHDFlowline_Network %>% cbind(Site.Agg=Site.data$Site.Agg[s])
  print(Site.data$Site.Agg[s])
  nhdata.all<-rbind(nhdata.all, nhdata)
}

site.schar<- as.data.frame(nhdata.all) %>%
  mutate(ppt.cm=ppt0001/10,
         elev=minelevsmo/100,
         logDA=log(totdasqkm)) %>%
  dplyr::select(Site.Agg, gnis_name, totdasqkm,logDA, lengthkm, streamorde,
                elev, slope,q0001e,v0001e,
                temp0001, ppt.cm) 

write.csv(site.schar, 'data/SiteChar.csv')

# land use data ====
#tutorial: https://mbjoseph.github.io/posts/2018-12-27-categorical-spatial-data-extraction-around-buffered-points-in-r/
library(raster); library(FedData); library(rgdal)

landprop.sum<-NULL
for(m in unique(Site.data$Site.Agg)){
  Sp.Site.data<-Site.data %>% ungroup() %>%
    dplyr::select(Site.Agg, Lat.cor, Long.cor) %>%
    filter(Site.Agg==m) %>%
    slice(1)
  print(Sp.Site.data$Site.Agg)
  coordinates(Sp.Site.data) <- ~Long.cor + Lat.cor
  proj4string(Sp.Site.data) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  #download land cover data; takes a while (~30s per tile)
  buff_site<-buffer(Sp.Site.data, width=100)
  # older code not working
  #nlcd_raster <- get_nlcd(buff_site,
  #                        label=paste(Sp.Site.data$Site.Agg),
  #                        extraction.dir = './EXTRACTIONS/',
  #                        dataset="landcover")
  # below replaces it
  nlcd_raster<-raster(paste0("./EXTRACTIONS/", Sp.Site.data$Site.Agg,
                             "/NLCD/", Sp.Site.data$Site.Agg, 
                      "_NLCD_2011_landcover.tif"))
  
  
  # reproject your data to raster's crs
  buff_site <- spTransform(buff_site, projection(nlcd_raster))
  
  landcover <- extract(nlcd_raster, buff_site)
  landprop<-as.data.frame(prop.table(table(landcover))) %>% 
    cbind(Site.Agg=Sp.Site.data$Site.Agg)
  landprop.sum<-bind_rows(landprop.sum, landprop)
}
landprop.sum %>% pull(Site.Agg) %>% unique()
site.landuse <- landprop.sum %>%
  mutate(ID=as.numeric(paste(landcover))) %>%
  left_join(nlcd_raster@data@attributes[[1]][,c(1,8)]) 
site.landuse %>% group_by(Site.Agg) %>% #checking it sums to one
  summarize(sumF=sum(Freq)) %>%
  arrange(sumF)
site.forest<-site.landuse %>% 
  dplyr::filter(NLCD.2011.Land.Cover.Class %in%
                c("Deciduous Forest","Evergreen Forest","Mixed Forest")) %>%
  mutate(Prop=Freq*100) %>%
  group_by(Site.Agg) %>%
  summarize(Fprop100=sum(Prop))
site.urban<-site.landuse %>% 
  dplyr::filter(NLCD.2011.Land.Cover.Class %in%
                c("Developed, High Intensity","Developed, Low Intensity",
                  "Developed, Medium Intensity","Developed, Open Space")) %>%
  mutate(Prop=Freq*100) %>%
  group_by(Site.Agg) %>%
  summarize(Uprop100=sum(Prop))
site.ag<-site.landuse %>% 
  dplyr::filter(NLCD.2011.Land.Cover.Class %in%
                  c("Hay/Pasture","Cultivated Crops")) %>%
  mutate(Prop=Freq*100) %>%
  group_by(Site.Agg, NLCD.2011.Land.Cover.Class) %>%
  summarize(Prop100=sum(Prop)) %>% group_by(Site.Agg) %>%
  pivot_wider(names_from = NLCD.2011.Land.Cover.Class,
              values_from = Prop100)

site.imp.landuse<-left_join(site.forest, site.urban) %>%
  full_join(site.ag)
write.csv(site.imp.landuse, 'data/LanduseT.csv')

# combine into site x environment table -----
# apply multiple regression to estimate a & b for temp regression
# a = 0.055*logDA - 0.004BaseFlow - 0.047*ln(elev) - 0.001(prec.cm) + 0.002Forest +0.993
# b = -0.62*logDA - 0.24T + 0.15HOF - 0.06Urban + 0.04Forest +9.8

site.env<-left_join(site.schar, site.imp.landuse, by='Site.Agg') %>%
  replace(is.na(.),0) %>%
  full_join(site.hyd, by='Site.Agg') %>%
  mutate(a = 0.055*logDA - 0.004*ml19 - 0.047*log(elev) - 0.001*ppt.cm + 0.002*Fprop100 +0.993,
         b = -0.62*logDA - 0.24*temp0001 - 0.06*Uprop100 + 0.04*Fprop100 + 9.8) %>%
  filter(!duplicated(.))

# pull met data down
Site.data.met<-Site.data %>%
  mutate(nMonths=(maxYear-minYear)*12) %>%
  group_by(HUC.8) %>% 
  arrange(desc(minYear)) %>% 
  slice(1)
sum(Site.data.met$nMonths)

library(rnoaa)
options("noaakey" = Sys.getenv("noaakey"))
#ncdc_locs_cats(datasetid="GSOY")
#ncdc_datatypes(datasetid = "GSOY")
for(u in 1:14){
  y.span<-seq.Date(as.Date(paste(Site.data.met$maxYear[u], "-01-01", sep="")),
                   as.Date(paste(Site.data.met$minYear[u], "-01-01", sep="")),
                   by="-9 years")
  y.span <-c(y.span, as.Date(paste(Site.data.met$minYear[u], "-01-01", sep="")))
  ly<-length(y.span)
  out_all<-NULL
  for(l in ly:2){
      out <- ncdc(datasetid='GSOY', 
                  locationid=paste('HUC:', Site.data.met$HUC.8[u], sep=""), 
                  datatypeid=c('TAVG','TMAX','TMIN','HDSD','CDSD'), 
                  startdate =paste(y.span[l]),
                  enddate=paste(y.span[l-1]),
                  limit=1000)
      out_all<-rbind(out_all, out$data) 
    }
    
    print(paste(Site.data.met$HUC.8[u],paste(u, l)))
    write.csv(out_all, file=paste("data/weather/",
                                paste(Site.data.met$HUC.8[u],
                                      ".csv", sep = ""), sep=""))
  
}

wfiles<-paste("data/weather/",
               list.files(path="data/weather/"), sep="")
weather<-NULL
for (i in 2:length(wfiles)){
  ws <- read.csv(wfiles[i]) %>% 
    mutate(HUC.8=substr(wfiles[i],14,21),
           Date.r=as.Date(date))
  weather<-rbind(weather, ws)
} 

#problem childs: 03160107;12100202
Site.data.met[which(is.na(match(unique(Site.data.met$HUC.8),
                                unique(weather$HUC.8)))),1:7]
weath.W2<-ncdc(datasetid='GSOY', locationid='HUC:031601', 
               datatypeid=c('TAVG','TMAX','TMIN','HDSD','CDSD'), 
               startdate='2014-01-01', enddate='2016-01-01',
               limit=1000)
weath.W2$data<-weath.W2$data %>% mutate(HUC.8="03160107", Date.r=as.Date(date))
weath.DW1<-ncdc(datasetid='GSOY', locationid='HUC:12100202', 
               datatypeid=c('TAVG','TMAX','TMIN','HDSD','CDSD'), 
               startdate='2006-01-01', enddate='2016-01-01',
               limit=1000)
weath.DW2<-ncdc(datasetid='GSOY', locationid='HUC:12100202', 
                datatypeid=c('TAVG','TMAX','TMIN','HDSD','CDSD'), 
                startdate='2004-01-01', enddate='2005-01-01',
                limit=1000)
weath.DW<-bind_rows(weath.DW1, weath.DW2) %>%
  mutate(HUC.8="12100202", Date.r=as.Date(date))
site.airT<-bind_rows(weather, weath.DW,weath.W2) %>%
  dplyr::select(date, datatype, station, value, HUC.8, Date.r) %>%
  filter(!is.na(datatype))

write.csv(site.airT, 'data/weather/weatherAll.csv')

pdf('data/weather/weathercheck.pdf')
ggplot(site.airT[site.airT$datatype=="TAVG",], aes(x=Date.r,y=value))+
    geom_point()+
    scale_x_date(date_labels="%y")+
    theme(axis.text.x = element_text(angle=30))+
  facet_wrap(~HUC.8, scales="free")
ggplot(weather[weather$datatype=="TMIN",], aes(x=Date.r,y=value))+
  geom_point()+
  scale_x_date(date_labels="%y")+
  theme(axis.text.x = element_text(angle=30))+
  facet_wrap(~HUC.8, scales="free")
ggplot(weather[weather$datatype=="TMAX",], aes(x=Date.r,y=value))+
  geom_point()+
  scale_x_date(date_labels="%y")+
  theme(axis.text.x = element_text(angle=30))+
  facet_wrap(~HUC.8, scales="free")
dev.off()

# apply regression to calculate water temperature
site.airT %>% 
  group_by(datatype, Date.r) %>%
  filter(datatype=="TMIN") %>%
  summarize(mean(value))

w.sum<-site.airT %>% 
  group_by(HUC.8, datatype, Date.r) %>%
  filter(value>0) %>%
  summarize(meanval=mean(value)) %>%
  spread(datatype,meanval) %>%
  right_join(Site.data) %>%
  left_join(site.env) %>%
  mutate(year=year(Date.r),
         wtavg=a*TAVG+b,
         wtmax=a*TMAX+b,
         wtmin=a*TMIN+b,
         Sp=case_when(Species %in% c("LCAR","LORN")~"LAMP",
                             Species =="APLI"~"APLI"),
         SpF=factor(Sp)) %>%
  dplyr::select(Site.Agg,SpF,HUC.8,year,
                TAVG,wtavg, TMAX,wtmax,TMIN,wtmin, minYear, maxYear) %>%
  group_by(Site.Agg, SpF) %>%
  summarize(Xwtavg=mean(wtavg, na.rm=T),
            Xwtmax=mean(wtmax, na.rm=T),
            Xwtmin=mean(wtmin, na.rm=T),
            maxy=max(year),
            miny=min(year),
            nobs=n(),
            minYear=mean(minYear),
            maxYear=mean(maxYear))


write.csv(w.sum, 'data/weatherSum.csv')

# see how well it predicted based on usgs gages 
#pull stream temperature from usgs gages
  
dailyT.all<-NULL
for(k in 1:40){
  dailyT <- readNWISdata(siteNumber=Site.data$USGS.Gage[k], service="dv",
                         parameterCd="00010",
                         startDate=paste(Site.data$minYear-1,"10-01", sep="-")[k],
                         endDate=paste(Site.data$maxYear,"09-30", sep="-")[k]) %>%
    mutate(Site.Agg=Site.data$Site.Agg[k],
           Species=Site.data$Species[k])
  dailyT.all<-bind_rows(dailyT.all, dailyT)
}
  
dailyT.all %>% group_by(site_no, Site.Agg) %>%
  summarize(meanT=mean(X_00010_00003, na.rm=T),
            nobs=n())

airTemp<-site.airT %>% filter(datatype=="TAVG") %>% 
  mutate(year=year(Date.r)) %>%
  group_by(HUC.8, year) %>%
  summarize(mAT=mean(value, na.rm=T)) %>%
  ungroup()%>%
  left_join(Site.data[,-1], by="HUC.8") %>%  
  left_join(site.env[!duplicated(site.env[,c(1)]),c(1,37,38)]) %>% 
  dplyr::select(Site.Agg,HUC.8,year,mAT,a,b) %>%
  filter(!duplicated(.))
checkTReg<-dailyT.all %>% left_join(Site.data[,-1]) %>%
  dplyr::select(Site.Agg, HUC.8,dateTime, X_00010_00003) %>%
  ungroup() %>%
  rename(ActWT=X_00010_00003) %>% 
  mutate(year=year(dateTime)) %>%
  group_by(HUC.8, year) %>%
  summarize(myWT=mean(ActWT, na.rm=T)) %>%
  left_join(airTemp) %>% 
  mutate(wtavg=round((a*mAT+b),2))%>%
  filter(!is.na(wtavg))

  
which(duplicated(checkTReg))
summary(lm(wtavg~myWT-1, data=checkTReg))
cor.test(checkTReg$wtavg,checkTReg$myWT) #0.794
length(unique(checkTReg$HUC.8));length(unique(checkTReg$year))

# look at the nutrient data available in USGS ---
head(Site.data)
nuttest<-NULL
for(n in 1:40){
test<-readNWISqw(Site.data$USGS.Gage[n], 
          "NUT", 
          startDate = paste(Site.data$minYear[n], "-01-01", sep=""), 
          endDate = paste(Site.data$maxYear[n], "-01-01", sep=""))
if(nrow(test)==0){
  test<-data.frame(site_no=Site.data$USGS.Gage[n], 
          sample_dt=NA, parm_cd=NA, result_va=NA)
  nuttest<-bind_rows(nuttest,test)
}else{
  test<-test %>% dplyr::select(site_no, sample_dt, parm_cd, result_va)
  nuttest<-bind_rows(nuttest,test)
}
}

View(nuttest %>% filter(is.na(result_va)) %>%
       group_by(site_no) %>% tally() %>%
  left_join(Site.data[,2:5], by=c("site_no"="USGS.Gage")))
