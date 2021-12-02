library(tidyverse);library(readxl);library(dplR)

#load in data
rw<-read_excel('data/growth.xlsx', sheet="rwl")
SlideData<-read_excel('data/growth.xlsx', sheet="SlideData")
ShellDim<-read_excel('data/!GrowthRawData.xlsx', sheet="SampleLengths")

sxs<-ShellDim %>% 
  left_join(SlideData, by=c('Shell.ID'='FileName'))%>%
  filter(is.na(PS)) %>%
  select(Shell.ID, Site.Agg, Species, RawName,PS)%>%
  mutate(Site.Agg=str_to_upper(Site.Agg))
rw.long<-rw %>% gather(Shell.ID, annuli.growth, -Year) %>%
  left_join(sxs, by=c("Shell.ID"='RawName')) %>%
  filter(is.na(PS))

View(rw.long %>% group_by(Site.Agg, Species) %>% 
  summarize(length(unique(Shell.ID))))
# NA because problem shells or duplicate pictures

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
            file=paste("data/apli_raw/excel/",
                       paste(names(apli.sites)[i],
                       ".csv", sep = ""), sep=""))
}

apfiles<-paste("data/apli_raw/physxd_fin/",
            list.files(path="data/apli_raw/physxd_fin/"), sep="")
apli_xd = lapply(apfiles, read.csv,row.names = "X")
names(apli_xd)<-gsub(".csv","",gsub("data/apli_raw/physxd_fin/","",apfiles))
#exporting each species at each site to raw
for(i in 1:length(apli_xd)){
  write.rwl(as.data.frame(apli_xd[[i]]),
            fname=paste("data/apli_raw/",paste(names(apli_xd)[i],
                                          ".raw", sep = ""), sep=""),
            format="compact")
}

#### Lampsilis ------
rw.long.lamp<-rw.long %>% 
  filter(Species == "LCAR") %>%
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
            file=paste("data/lamp_raw/excel/",paste(names(lamp.sites)[i],
                                               ".csv", sep = ""), sep=""))
}

#identyifying just those that are lcar
rw.long %>% filter(Species=='LCAR') %>% pull(Site.Agg) %>% unique()
rw.long %>% filter(Species=='LORN') %>% pull(Site.Agg) %>% unique()

lcfiles<-paste("data/lamp_raw/physxd_fin/",
               list.files(path="data/lamp_raw/physxd_fin/"), sep="") 
lcfiles<-lcfiles[-c(1,12,15,16)]
for (i in 1:length(lcfiles)){
  lamp_xd = lapply(lcfiles, read.csv,row.names = "X")
} 
names(lamp_xd)<-gsub(".csv","",gsub("data/lamp_raw/physxd_fin/","",lcfiles))

#exporting each species at each site to raw
for(i in 1:length(lamp_xd)){
  write.rwl(as.data.frame(lamp_xd[[i]]),
            fname=paste("data/lamp_raw/",paste(names(lamp_xd)[i],
                                               ".raw", sep = ""), sep=""),
            format="compact")
}

# Tritigonia verrucosa ------
rw.long.qver<-rw.long %>% 
  filter(Species=="QVER") %>%
  select(-Shell.ID.y, -PS)

qver.sites<-NULL
for(u in 1:length(unique(rw.long.qver$Site.Agg))) {
  qver.sites[[u]]<-rw.long.qver %>% 
    filter(Site.Agg==unique(Site.Agg)[u]) %>%
    spread(Shell.ID, annuli.growth)%>% 
    column_to_rownames("Year") %>% 
    select(-Species,-Site.Agg) 
}
names(qver.sites)<-unique(rw.long.qver$Site.Agg)

# exporting excel to do physical crossdating
for(i in 1:length(qver.sites)){
  write.csv(as.data.frame(qver.sites[[i]]),
            file=paste("data/qver_raw/excel/",paste(names(qver.sites)[i],
                                                    ".csv", sep = ""), sep=""))
}
