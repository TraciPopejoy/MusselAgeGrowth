library(maps);library(maptools);library(sf)
library(rgdal);library(hydroMap)


#bigR<-getFlowLines(c(-102.5,-71.8, 27,50), 7, filePath="7thOrderR")
bigR<-read_sf(dsn="data/7thOrderR")

source('1-BCShellLengths.R')
Site.data<-AxL %>% group_by(Species, Site) %>%
  rename(SiteID=Site)%>%
  summarize(maxYear=max(Year),
            minYear=min(Year)) %>% 
  left_join(SiteID, by="SiteID") %>%
  dplyr::select(Species,Site.Agg, HUC.8,HUC.12,USGS.Gage, maxYear, minYear, Lat.cor, Long.cor) %>%
  filter(!duplicated(Lat.cor))%>%
  rowwise()
SiteID<-read_excel('data/!GrowthRawData.xlsx', sheet="Location") %>%
  group_by(Site.Agg) %>%
  mutate(Lat.cor=mean(Latitude),
         Long.cor=mean(Longitude))


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

# publication ready plot -----
library(tidyverse)
usa <- st_as_sf(map("usa", plot=F, fill=T))
states <- st_as_sf(map("state", c('texas','oklahoma','kansas','south dakota','north dakota',
                             'arkansas','louisiana','indiana','kentucky','nebraska',
                             'minnesota','iowa','michigan','missouri', 'illinois', 'wisconsin',
                             'ohio', 'mississippi','alabama','tennessee', 'georgia','new york',
                             'west virginia','maryland','virginia','north carolina',
                             'south carolina','florida','new jersey','delaware', 
                             'pennsylvania'),
                       plot = FALSE, fill = TRUE))

#Data: National Hydrological Database 
#Metadata: Downloaded 9November2016
interior<-bind_rows(read_sf(dsn="C:/Users/Owner/Documents/GISfile/KiamichiShape", layer="NHDFlowline") %>%
                      filter(grepl('Kiamichi', GNIS_NAME)),
                    read_sf(dsn="C:/Users/Owner/Documents/GISfile/LittleShape", layer="NHDFlowline") %>%
                      filter(grepl('Little River', GNIS_NAME))) %>%
  mutate(region="Interior Highlands")
intWS<-rbind(read_sf(dsn="C:/Users/Owner/Documents/GISfile/KiamichiShape", layer="WBDHU8") %>% 
               filter(grepl("Kiamichi",NAME)|
                      grepl("Upper Little",NAME))) %>%
  mutate(region="Interior Highlands")
mobile<-read_sf(dsn="C:/Users/Owner/Documents/GISfile/SipseyShape", layer="NHDFlowline") %>%
                  filter(grepl("Sipsey", GNIS_NAME)) %>%
                  mutate(region="Mobile Basin")
mobWS<-read_sf(dsn="C:/Users/Owner/Documents/GISfile/SipseyShape", layer="WBDHU8") %>%
  filter(grepl("Sipsey", NAME)) %>%
  mutate(region="Mobile Basin")
txts<-read_sf(dsn="C:/Users/Owner/Documents/GISfile/TexasShape", layer="NHDFlowline") %>%
  filter(GNIS_NAME %in% c("Guadalupe River","Colorado River"))%>%
  mutate(region="Western Gulf")
txWS<-read_sf(dsn="C:/Users/Owner/Documents/GISfile/TexasShape", layer="WBDHU8") %>%
  filter(grepl("Guadalupe", NAME) |
           grepl("Colorado", NAME) |
           HUC8 %in% c("12090201", "12090205")) %>%
  mutate(region="Western Gulf")
stlaw<-bind_rows(read_sf(dsn="C:/Users/Owner/Documents/GISfile/GrandShape", layer="NHDFlowline") %>%
                   filter(grepl("Grand",GNIS_NAME)),
                 read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_05010004_HU8_Shape/Shape", layer="NHDFlowline") %>%
                           filter(GNIS_NAME=="French Creek")) %>%
  mutate(region="St Lawrence-Great Lakes")
stlawWS<-bind_rows(read_sf(dsn="C:/Users/Owner/Documents/GISfile/GrandShape", layer="WBDHU8"), 
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_05010004_HU8_Shape/Shape", layer="WBDHU8")) %>%
  mutate(region="St Lawrence-Great Lakes")
upmis<- bind_rows(read_sf(dsn="C:/Users/Owner/Documents/GISfile/StCroixShape", layer="NHDFlowline") %>%
                    filter(grepl('Croix', GNIS_NAME)),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/MississippiShape", layer="NHDFlowline") %>%
                    filter(grepl("Mississippi", GNIS_NAME)),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07090006_HU8_Shape/Shape", layer = "NHDFlowline") %>%
                    filter(GNIS_NAME=="Kishwaukee River"),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07080101_HU8_Shape/Shape", layer = "NHDFlowline") %>%
                    filter(grepl("Mississippi",GNIS_NAME)),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07140104_HU8_Shape/Shape", layer="NHDFlowline") %>%
                    filter(grepl("Big River", GNIS_NAME)),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07130011_HU8_Shape/Shape", layer = "NHDFlowline") %>%
                    filter(grepl('Illinois', GNIS_NAME)) ,
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07110006_HU8_Shape/Shape", layer="NHDFlowline") %>%
                    filter(grepl("Salt River", GNIS_NAME)),
                  read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07120002_HU8_Shape/Shape", layer="NHDFlowline") %>%
                    rename(GNIS_NAME=GNIS_Name) %>% st_zm() %>%
                    filter(grepl("Iroquois", GNIS_NAME))) %>%
  mutate(region="Upper Mississippi")
upmisWS<-bind_rows(read_sf(dsn="C:/Users/Owner/Documents/GISfile/StCroixShape", layer="WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/MississippiShape", layer="WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07090006_HU8_Shape/Shape", layer = "WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07080101_HU8_Shape/Shape", layer = "WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07140104_HU8_Shape/Shape", layer="WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07130011_HU8_Shape/Shape", layer = "WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07110006_HU8_Shape/Shape", layer="WBDHU8"),
                   read_sf(dsn="C:/Users/Owner/Documents/GISfile/NHD_H_07120002_HU8_Shape/Shape", layer="WBDHU8")) %>%
  mutate(region="Upper Mississippi")

Lsites<- AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species %in% c("LCAR")) %>%
  pull(Site.Agg) %>% unique()
Asites<- AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species=="APLI") %>%
  pull(Site.Agg) %>% unique()
TSsiteLocations<-SiteID %>%
  right_join(Site.data) %>%
  ungroup() %>%
  filter(!duplicated(.)) 
jittPTS<-TSsiteLocations %>%
  filter(!(Site.Agg %in% Asites[Asites %in% Lsites])) %>%
  mutate(Long.cor=ifelse(Species=="LCAR", Long.cor-.13, Long.cor+.13)) %>%
  st_as_sf(coords=c('Long.cor','Lat.cor'), crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))

normPTS<-TSsiteLocations %>%
  filter(Site.Agg %in% Asites[Asites %in% Lsites]) %>%
  st_as_sf(coords=c('Long.cor','Lat.cor'), crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))

ggplot(jittPTS) +geom_sf(aes(color=Site.Agg))

colorvector<-c('#b3bd28','#ca6833','#00822d','#0069d6','#c64579')

library(ggplot2); library(ggrepel) 
ggplot()+
  geom_sf(data = states, color="black",fill="white")+
  geom_sf(data=intWS, color=NA,aes( fill=region), alpha=.5)+
  geom_sf(data=mobWS,color=NA, aes( fill=region),alpha=.5)+
  geom_sf(data=stlawWS,color=NA, aes( fill=region),alpha=.5)+
  geom_sf(data=upmisWS, color=NA,aes( fill=region),alpha=.5)+
  geom_sf(data=txWS, color=NA,aes( fill=region),alpha=.5)+
  geom_sf(data = bigR, color="dark grey") +
  geom_sf(data=txts, aes(color=region))+
  geom_sf(data=mobile, aes(color=region))+
  geom_sf(data=interior, aes(color=region))+
  geom_sf(data=stlaw, aes(color=region))+
  geom_sf(data=upmis, aes(color=region))+
  geom_sf(data=jittPTS, aes(shape=Species))+
  geom_sf(data=normPTS, aes(shape='Both'))+
  scale_shape_manual(name='Species', values=c(24,23,25))+
  #geom_text_repel(data=tidysites, aes(x=Long.cor, y=Lat.cor, 
  #                                    label=Site.Agg), size=3)+
  scale_colour_manual(name="Biogeographic Province",
                      breaks=c("Upper Mississippi","St Lawrence-Great Lakes",
                               "Interior Highlands","Mobile Basin","Western Gulf"),
                      values=colorvector[c(4,3,1,2,5)],
                      aesthetics=c('fill','color'))+
  labs(x='Longitude', y='Latitude')+theme_minimal()+
  theme(legend.position ="bottom", 
        legend.box.background = element_rect(fill="white"))+
  guides(fill = guide_legend(nrow = 2, title.position="top"),
         col = guide_legend(nrow = 2, title.position="top"),
         shape = guide_legend(nrow = 2, title.position="top"))
ggsave('figures/mapFull2_2111.svg', width=6, height=6)

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
