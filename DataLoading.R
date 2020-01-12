library(tidyverse);library(readxl)
library(dplR)

#load in data
rw<-read_excel('data/growth.xlsx', sheet="rwl")
SlideData<-read_excel('data/growth.xlsx', sheet="SlideData")
ShellDim<-read_excel('data/!GrowthRawData.xlsx', sheet="SampleLengths")

sxs<-ShellDim %>% 
  left_join(SlideData, by=c('Shell.ID'='FileName'))%>%
  filter(is.na(PS)) %>%
  select(Shell.ID, Site.Agg, Species, RawName,PS)
rw.long<-rw %>% gather(Shell.ID, annuli.growth, -Year) %>%
  left_join(sxs, by=c("Shell.ID"='RawName')) %>%
  filter(is.na(PS))

View(rw.long %>% group_by(Site.Agg, Species) %>% 
  summarize(length(unique(Shell.ID))))

rw.long.apli<-rw.long %>% filter(Species=="APLI") %>%
  select(-Shell.ID.y, -PS)

apli.sites<-NULL
for(u in 1:length(unique(rw.long.apli$Site.Agg))) {
  apli.sites[[u]]<-rw.long.apli %>% 
    filter(Site.Agg==unique(Site.Agg)[u]) %>%
    spread(Shell.ID, annuli.growth)%>% 
    column_to_rownames("Year") %>% 
    select(-Species,-Site.Agg) 
}
names(apli.sites)<-unique(rw.long.apli$Site.Agg)

# exporting excel to do physical crossdating
for(i in 1:length(apli.sites)){
  write.csv(as.data.frame(apli.sites[[i]]),
            file=paste("apli_raw/excel/",
                       paste(names(apli.sites)[i],
                       ".csv", sep = ""), sep=""))
}

temp<-paste("data/apli_xd/",list.files(path="data/apli_xd"), sep="")
for (i in 1:length(temp)){
  apli_xd = lapply(temp, read.csv,row.names = "X")
} 
names(apli_xd)<-gsub(".csv","",gsub("data/apli_xd/","",temp))

#exporting each species at each site to raw
for(i in 1:length(apli_xd)){
  write.rwl(as.data.frame(apli_xd[[i]]),
            fname=paste("data/apli_xd/",paste(names(apli_xd)[i],
                                          ".raw", sep = ""), sep=""))
}

#### Lampsilis ------
rw.long.lamp<-rw.long %>% 
  filter(Species %in% c("LCAR","LORN")) %>%
  select(-Shell.ID.y, -PS)

lamp.sites<-NULL
for(u in 1:length(unique(rw.long.lamp$Site.Agg))) {
  lamp.sites[[u]]<-rw.long.lamp %>% 
    filter(Site.Agg==unique(Site.Agg)[u]) %>%
    spread(Shell.ID, annuli.growth)%>% 
    column_to_rownames("Year") %>% 
    select(-Species,-Site.Agg) 
}
names(lamp.sites)<-unique(rw.long.lamp$Site.Agg)

# exporting excel to do physical crossdating
for(i in 1:length(lamp.sites)){
  write.csv(as.data.frame(lamp.sites[[i]]),
            file=paste("lamp_raw/excel/",paste(names(lamp.sites)[i],
                                               ".csv", sep = ""), sep=""))
}
