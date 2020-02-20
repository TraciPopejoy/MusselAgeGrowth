library(maps);library(maptools);library(sp)
library(rgdal);library(hydroMap)

#Data: National Hydrological Database 
#Metadata: Downloaded 9November2016
KiamichiHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/KiamichiShape", layer="NHDFlowline")
LittleHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/LittleShape", layer="NHDFlowline")
SipseyHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/SipseyShape", layer="NHDFlowline")
GrandHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/GrandShape", layer="NHDFlowline")
StCroixHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/StCroixShape", layer="NHDFlowline")
MissMNHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/MississippiShape", layer="NHDFlowline")
TexasHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/TexasShape", layer="NHDFlowline")
KishwakHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07090006_HU8_Shape/Shape", layer = "NHDFlowline")
IowaHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07080101_HU8_Shape/Shape", layer = "NHDFlowline")
BigRiverHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07140104_HU8_Shape/Shape", layer="NHDFlowline")
IllinoisHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07130011_HU8_Shape/Shape", layer = "NHDFlowline")
ElkRiverHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07110006_HU8_Shape/Shape", layer="NHDFlowline")
FrenchHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_05010004_HU8_Shape/Shape", layer="NHDFlowline")
KankakeeHUC8<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07120002_HU8_Shape/Shape", layer="NHDFlowline")

WS1<-rbind(KiamichiHUC8[grep("Kiamichi",KiamichiHUC8$GNIS_NAME),],
           LittleHUC8[grep("Little River",LittleHUC8$GNIS_NAME),])
WS2<-rbind(WS1, SipseyHUC8[grep("Sipsey", SipseyHUC8$GNIS_NAME),-c(13:15)])
WS3<-rbind(WS2, GrandHUC8[grep("Grand",GrandHUC8$GNIS_NAME),-c(13:15)])
WS4<-rbind(WS3, StCroixHUC8[grep("Croix", StCroixHUC8$GNIS_NAME),-c(13:15)])
WS5<-rbind(WS4, MissMNHUC8[grep("Mississippi", MissMNHUC8$GNIS_NAME),-c(13:15)])
WS6<-rbind(WS5, TexasHUC8[grep("Guadalupe River", TexasHUC8$GNIS_NAME),-c(13:15)])
WS7<-rbind(WS6, KishwakHUC8[KishwakHUC8$GNIS_NAME=="Kishwaukee River",-c(13:15)])
WS8<-rbind(WS7, IowaHUC8[grep("Mississippi",IowaHUC8$GNIS_NAME),-c(13:15)])
WS9<-rbind(WS8, BigRiverHUC8[grep("Big River", BigRiverHUC8$GNIS_NAME),-c(13:15)])
WS10<-rbind(WS9,IllinoisHUC8[grep("Illinois",IllinoisHUC8$GNIS_NAME),-c(13:15)])
WS11<-rbind(WS10, ElkRiverHUC8[grep("Elk Fork Salt River", ElkRiverHUC8$GNIS_NAME), -c(13:15)])
WS12<-rbind(WS11, LittleHUC8[grep("Glover",LittleHUC8$GNIS_NAME),])
WS13<-rbind(WS12, TexasHUC8[grep("Colorado River", TexasHUC8$GNIS_NAME), -c(13:15)])
WS14<-rbind(WS13, TexasHUC8[grep("Neches",TexasHUC8$GNIS_NAME), -c(13:15)])
WS15<-rbind(WS14, FrenchHUC8[FrenchHUC8$GNIS_NAME=="French Creek", -c(13:15)])


bigR<-getFlowLines(c(-102.5,-71.8, 27,50), 7, filePath="7thOrderR")
bigR<-readOGR(dsn="data/7thOrderR")

TSsiteLocations<-SiteID 
coordinates(TSsiteLocations)<- ~Long.cor+Lat.cor

library(wesanderson)

musbioreg<-data.frame(regions=c(rep("Interior Highlands",2),
                     "Mobile Basin",
                     rep("St Lawrence-Great Lakes",2),
                     rep("Upper Mississippi",7),
                     rep("Western Gulf",2)),
           River=c('Kiamichi River', 'Little River', 'Sipsey River',
                  'French Creek', 'Grand River', 'Big River',
                   'Elk Fork Salt River', 'Illinois River',
                   'Kankakee River','Kishwaukee River',
                   'Mississippi River', 'Saint Croix River',
                   'Colorado River', 'Guadalupe River'))


makemap<-function(species, lab=F, point=T){
  SPP<-NULL
  for(s in 1:length(species)){
    SPP1<-ShellDim[ShellDim$Species==species[s],]
    SPP<-rbind(SPP,SPP1)
  }
  sites<-subset(TSsiteLocations, TSsiteLocations$SiteID %in% unique(SPP$Site))
    st<-map('state', fill=T, col="transparent", plot=F)
  IDs<-sapply(strsplit(st$names, ":"), function(x) x[1])
  states_sp<-map2SpatialPolygons(st, IDs=IDs,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))  
  pointsSP<-SpatialPoints(sites[,c(3,2)], 
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
  indices<-over(pointsSP, states_sp)
  stateNames<-sapply(states_sp@polygons, function(x) x@ID)
  states<-stateNames[indices]
  par(mar=c(.5,.5,.5,.5))
  map('state', region=states, col="darkgrey", lwd=2)
  map('state', region = c('texas','oklahoma','kansas','south dakota','north dakota',
                          'arkansas','louisiana','indiana','kentucky','nebraska',
                          'minnesota','iowa','michigan','missouri', 'illinois', 'wisconsin',
                          'ohio', 'mississippi','alabama','tennessee', 'georgia','new york',
                          'west virginia','maryland','virginia','north carolina',
                          'south carolina','florida','new jersey','delaware'), 
      lwd=2.5, col="darkgray", add=T) 
  plot(bigR,col="lightgrey", add=T)
  for(j in 1:length(unique(sites$River))){
    ind<-unique(sites$River)[j]
    plot(WS15[grep(ind,WS15$GNIS_NAME),], col=sites$color[sites$River==ind], lwd=6, add=T)
  }
  if(point==T) {
    plot(sites, pch=19, add=T)
  }
  if(lab==T) {
    pointLabel(as.vector(coordinates(sites)[,1]), 
             as.vector(coordinates(sites)[,2]), labels=sites$SiteID)
  }
}

makemap("APLI")
makemap(c("LCAR","LORN"))
makemap("QVER")
makemap(c("FFLA","FCER"))
makemap(c("QPUS","QASP", "QMOR"))
makemap("OREF")
map.scale(-95,28, relwidth = .23,ratio=F)
legend("topleft", 
       legend=unique(SPP$River),
        col=c("black"), lty=c(0,0), 
        pch=c(19),lwd=3)
legend("topleft",
       legend=c("Sites","Kiamichi, OK","Sipsey, AL","St. Croix, MN","Mississippi, MN","Guadalupe, TX"),
       col=c("black",col[c(1:5)]),lty=c(0,1,1,1,1,1), pch=c(19,NA,NA,NA,NA,NA,NA,NA,NA), lwd=3)
points(TSsiteLocations$Longitude[c(1,2,7,11:14,17,5)], TSsiteLocations$Latitude[c(1,2,7,11:14,17,5)], 
     pch=19,cex=.8)
text(TSsiteLocations$Longitude[c(7,20,17,5)], TSsiteLocations$Latitude[c(7,20,17,5)], 
     label=TSsiteLocations$SiteID[c(7,20,17,5)],cex=.8, pos=2)
text(TSsiteLocations$Longitude[1:2], TSsiteLocations$Latitude[1:2], 
     label=TSsiteLocations$SiteID[1:2],cex=.8, pos=4)

# publication ready plot -----
usa <- map_data("usa")
states <- map_data("state",c('texas','oklahoma','kansas','south dakota','north dakota',
                             'arkansas','louisiana','indiana','kentucky','nebraska',
                             'minnesota','iowa','michigan','missouri', 'illinois', 'wisconsin',
                             'ohio', 'mississippi','alabama','tennessee', 'georgia','new york',
                             'west virginia','maryland','virginia','north carolina',
                             'south carolina','florida','new jersey','delaware', 
                             'pennsylvania'))
interior<-rbind(KiamichiHUC8[grep("Kiamichi",KiamichiHUC8$GNIS_NAME),],
           LittleHUC8[grep("Little River",LittleHUC8$GNIS_NAME),]) %>%
  broom::tidy() %>% mutate(region="Interior Highlands")
mobile<-SipseyHUC8[grep("Sipsey", SipseyHUC8$GNIS_NAME),-c(13:15)] %>%
  broom::tidy() %>% mutate(region="Mobile Basin")
stlaw<- rbind(GrandHUC8[grep("Grand",GrandHUC8$GNIS_NAME),-c(13:15)],
              FrenchHUC8[FrenchHUC8$GNIS_NAME=="French Creek", -c(13:15)]) %>%
  broom::tidy() %>% mutate(region="St Lawrence-Great Lakes")
upmis<- rbind(StCroixHUC8[grep("Croix", StCroixHUC8$GNIS_NAME),-c(13:15)],
              MissMNHUC8[grep("Mississippi", MissMNHUC8$GNIS_NAME),-c(13:15)]) %>%
  rbind(KishwakHUC8[KishwakHUC8$GNIS_NAME=="Kishwaukee River",-c(13:15)]) %>%
  rbind(IowaHUC8[grep("Mississippi",IowaHUC8$GNIS_NAME),-c(13:15)]) %>%
  rbind(BigRiverHUC8[grep("Big River", BigRiverHUC8$GNIS_NAME),-c(13:15)]) %>%
  rbind(ElkRiverHUC8[grep("Salt River", ElkRiverHUC8$GNIS_NAME), -c(13:15)]) %>%
  rbind(IllinoisHUC8[grep("Illinois",IllinoisHUC8$GNIS_NAME),-c(13:15)])%>%
  rbind(KankakeeHUC8[grep("Iroquois",KankakeeHUC8$GNIS_NAME), -c(13:15)]) %>%
  broom::tidy() %>% mutate(region="Upper Mississippi")
txts<-rbind(TexasHUC8[grep("Guadalupe River", TexasHUC8$GNIS_NAME),-c(13:15)],
           TexasHUC8[grep("Colorado River", TexasHUC8$GNIS_NAME), -c(13:15)])%>%
  broom::tidy() %>% mutate(region="Western Gulf")

tssites<-TSsiteLocations[TSsiteLocations$Site.Agg %in% c(Asites, Lsites),-c(1,3,4,11)]
tssites<-tssites[which(!duplicated(tssites@data)),]
tidysites<-data.frame(tssites) 

#usstreams<-broom::tidy(bigR)
library(ggplot2); library(ggrepel) 
ggplot()+
  geom_polygon(data = states, aes(x = long, y = lat, group=group), 
               color="black",fill="white")+
  #geom_path(data = usstreams, aes(x=long, y = lat, group = group), 
  #          color="darkcyan", size=0.4) +
  geom_path(data=txts, aes(x=long, y=lat, group=group, color=region),
            size=1.5)+
  geom_path(data=mobile, aes(x=long, y=lat, group=group, color=region),
            size=1.5)+
  geom_path(data=interior, aes(x=long, y=lat, group=group, color=region),
            size=1.5)+
  geom_path(data=stlaw, aes(x=long, y=lat, group=group, color=region),
            size=1.5)+
  geom_path(data=upmis, aes(x=long, y=lat, group=group, color=region),
            size=1.5)+
  geom_point(data=tidysites, aes(x=Long.cor, y=Lat.cor))+
  #geom_text_repel(data=tidysites, aes(x=Long.cor, y=Lat.cor, 
  #                                    label=Site.Agg), size=3)+
  scale_colour_manual(name="Biogeographic Province",
                      breaks=c("Upper Mississippi","St Lawrence-Great Lakes",
                               "Interior Highlands","Mobile Basin","Western Gulf"),
                      values=wes_palette("Zissou1")[c(3,4,1,2,5)])+
  labs(x='Longitude', y='Latitude')+theme_minimal()+
  theme(legend.position ="bottom", 
        legend.box.background = element_rect(fill="white"))+
  guides(col = guide_legend(nrow = 2, title.position="top"))+
  coord_cartesian()
ggsave('figures/map.tiff', width=5, height=5.8)

# pulling out the HUC12 to access weather data ----
KiamichiHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/KiamichiShape", layer="WBDHU12")
LittleHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/LittleShape", layer="WBDHU12")
SipseyHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/SipseyShape", layer="WBDHU12")
GrandHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/GrandShape", layer="WBDHU12")
StCroixHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/StCroixShape", layer="WBDHU12")
MissMNHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/MississippiShape", layer="WBDHU12")
TexasHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/TexasShape", layer="WBDHU12")
KishwakHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07090006_HU8_Shape/Shape", layer = "WBDHU12")
IowaHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07080101_HU8_Shape/Shape", layer = "WBDHU12")
BigRiverHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07140104_HU8_Shape/Shape", layer="WBDHU12")
IllinoisHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07130011_HU8_Shape/Shape", layer = "WBDHU12")
ElkRiverHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07110006_HU8_Shape/Shape", layer="WBDHU12")
FrenchHUC12<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_05010004_HU8_Shape/Shape", layer="WBDHU12")

huc12<-rbind(KiamichiHUC12,LittleHUC12,SipseyHUC12,GrandHUC12,StCroixHUC12,
          MissMNHUC12,TexasHUC12,KishwakHUC12,IowaHUC12,BigRiverHUC12,
          IllinoisHUC12,ElkRiverHUC12,FrenchHUC12)
SpSite<-Site.data
coords(SpSite) <- 