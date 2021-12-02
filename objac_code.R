library(tidyverse); library(readxl); library(sf)

crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

ObJac_pts<-read_xlsx('C:/Users/Owner/Downloads/2021_SiteCoordinates.xlsx') %>%
  st_as_sf(coords=c('Long_dec', 'Lat_dec'), 
           crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))
st_layers('C:/Users/Owner/Documents/GISfile/SEOkNHD')
huc12<-read_sf(dsn='C:/Users/Owner/Documents/GISfile/SEOkNHD/Shape',
        layer='WBDHU12')
#seOK<-read_sf(dsn='C:/Users/Owner/Documents/GISfile/SEOkNHD/Shape',
#              layer="NHDFlowline")

objp<-ObJac_pts %>% 
  st_transform(st_crs(huc12))%>%
  st_join(huc12)
unique(substr(objp$HUC12, 1,8))
ws_bound<-huc12 %>% mutate(huc8=substr(HUC12, 1, 8)) %>%
  filter(huc8 %in%c("11140107", "11140108")) %>%
  st_union()

#tutorial: https://mbjoseph.github.io/posts/2018-12-27-categorical-spatial-data-extraction-around-buffered-points-in-r/
library(raster); library(FedData); library(stars)
library(nhdplusTools)

nlcd2001<-get_nlcd(ws_bound, label="SEOK2001", year=2001)
nlcd2001stars<-st_as_stars(nlcd2001)

rm(nlcd2001, huc12, seOK)
buff_out<-NULL
start.time<-Sys.time()
for(m in objp$Paper_ID){
  pt<-objp %>% filter(Paper_ID==m)
  
  start_coid<-discover_nhdplus_id(pt)
  nldi_feature <- list(featureSource = "comid", 
                       featureID = start_coid)
  pt_nhd<-get_nldi_feature(nldi_feature)
  
  for(scale in c(1000, 1)){
    flowline_nldi <- navigate_nldi(nldi_feature, 
                                   mode = "upstreamTributaries", 
                                   distance_km = scale)
    fl_buf_small<-st_geometry(flowline_nldi$UT) %>% 
      st_transform(crs.albers) %>%
      st_buffer(.1) %>% st_union() %>% 
      st_transform(st_crs(nlcd2001stars)) 
    rm(flowline_nldi)
    fl_raster_small<-nlcd2001stars[fl_buf_small] 
    
    ex_res<-fl_raster_small$SEOK2011_NLCD_Land_Cover_2011 %>%
      summary()
    
    buff_out<-bind_rows(buff_out,
                           data.frame(n_cell=ex_res,
                                      type=names(ex_res)) %>%
                             filter(n_cell!=0, type != "NA's") %>%
                             mutate(land_cover_percent=n_cell/sum(n_cell)*100,
                                    Paper_ID=m,
                                    scale=scale))
  }
}
end.time<-Sys.time() 
end.time-start.time  
write.csv(buff_out, 'SE_OK_landuse_survey2011.csv')

buff_out %>%
  group_by(Paper_ID, scale) %>%
  select(-n_cell) %>%
  pivot_wider(names_from=type, values_from=land_cover_percent) %>%
  mutate(spatial_scale=recode(scale, '1000'='100m buffer across upstream tributaries',
                              '1'='100m buffer 1km upstream')) %>%
  select(Paper_ID, spatial_scale,scale, everything()) %>%
  write.csv('SE_OK_Survey_landuse_percents2011.csv')


huc_out<-NULL
for(m in objp$Paper_ID){
  pt<-objp %>% filter(Paper_ID==m)
  start_coid<-discover_nhdplus_id(pt)
  nldi_feature <- list(featureSource = "comid", 
                       featureID = start_coid)
  pt_nhd<-get_nldi_feature(nldi_feature)
  basin<-get_nldi_basin(nldi_feature) %>% 
      st_transform(st_crs(nlcd2001stars)) 
  basin_rast<-nlcd2001stars[basin] 
  ex_res<-basin_rast$SEOK2001_NLCD_Land_Cover_2001 %>%
      summary()
    
    huc_out<-bind_rows(huc_out,
                        data.frame(n_cell=ex_res,
                                   type=names(ex_res)) %>%
                          filter(n_cell!=0, type != "NA's") %>%
                          mutate(land_cover_percent=n_cell/sum(n_cell)*100,
                                 Paper_ID=m,
                                 scale='watershed'))
    print(m)
}
huc_out %>%
  mutate(year=2001) %>% 
  write.csv('SE_OK_landuse_watershed_2001.csv')

## bring all years back together ---
landuse_long<-NULL
for(i in list.files(pattern='SE_OK_land')){
  landuse_long<-bind_rows(landuse_long, 
                          read.csv(i) %>%
                            mutate(year=gsub('.csv','',
                                             gsub('SE_OK_landuse_watershed_','',
                                             gsub('SE_OK_landuse_survey','',i))),
                                   spatial_scale=recode(scale, '1000'='100m buffer across upstream tributaries',
                                                               '1'='100m buffer 1km upstream')) %>%
                                     select(-scale))
  
}
View(landuse_long)

landuse_long %>% select(-X, -n_cell) %>%
  group_by(Paper_ID, year, spatial_scale) %>%
  pivot_wider(names_from=type, values_from=land_cover_percent) %>%
  arrange(Paper_ID, spatial_scale) %>%
  write.csv('SE OK landuse for long term data.csv', row.names=F)

nrow(objp)*3*3
# Example plot


ex_pt<-objp %>% sample_n(1)
start_coid<-discover_nhdplus_id(ex_pt)
nldi_feature <- list(featureSource = "comid", 
                     featureID = start_coid)
basin<-get_nldi_basin(nldi_feature) 
flowline_nldi_1000 <- navigate_nldi(nldi_feature, 
                               mode = "upstreamTributaries", 
                               distance_km = 1000)
flow_1000<-st_geometry(flowline_nldi_1000$UT) %>% 
  st_transform(crs.albers) %>%
  st_buffer(.1) %>% st_union() 
flowline_nldi_1 <- navigate_nldi(nldi_feature, 
                                    mode = "upstreamTributaries", 
                                    distance_km = 1)
flow_1<-st_geometry(flowline_nldi_1$UT) %>% 
  st_transform(crs.albers) %>%
  st_buffer(.1) %>% st_union() 
ggplot()+
  geom_sf(data=basin)+
  geom_sf(data=flow_1000, color='darkgrey', fill='darkgrey')+
  geom_sf(data=flow_1, color='red', fill='red')+
  geom_sf(data=ex_pt, size=2)+
  geom_sf_text(data=ex_pt, aes(label=Paper_ID), 
               nudge_y=0.021, size=3)+
  labs(y='Longitude',x='Latitude')+
  theme_classic()+
  theme(axis.text = element_text(size=8),
        axis.title=element_text(size=8))
ggsave('SE_OK_landuse_example.jpg', width=3.5, height = 3.5)

