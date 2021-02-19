library(readxl);library(ggplot2); library(tidyverse)
#library(dplR)
#load in data
rw<-read_excel('data/growth.xlsx', sheet="rwl")
SlideData<-read_excel('data/growth.xlsx', sheet="SlideData")
ShellDim<-read_excel('data/!GrowthRawData.xlsx', sheet="SampleLengths")
SiteID<-read_excel('data/!GrowthRawData.xlsx', sheet="Location") %>%
  group_by(Site.Agg) %>%
  mutate(Lat.cor=mean(Latitude),
         Long.cor=mean(Longitude))

# pull in adjusted, final crossdated shell growth annuli
apfiles1<-paste("data/apli_raw/physxd_fin/",
               list.files(path="data/apli_raw/physxd_fin/"), sep="")
lcfiles1<-paste("data/lamp_raw/physxd_fin/",
                list.files(path="data/lamp_raw/physxd_fin/"), sep="")
AxW1<-NULL
for (i in 1:length(apfiles1)){
  axws <- read.csv(apfiles1[i]) %>% gather(ShellIDy,AnGrowth,-X)
  AxW1<-rbind(AxW1, axws)
} 

for (i in 1:length(lcfiles1)){
  axws <- read.csv(lcfiles1[i]) %>% gather(ShellIDy,AnGrowth,-X)
  AxW1<-rbind(AxW1, axws)
} 

# use annual growth annuli to build age x width table
AxW<-AxW1 %>% filter(!is.na(AnGrowth)) %>% rename(Year=X) %>%
  group_by(ShellIDy) %>%
  mutate(RawName=case_when(substr(ShellIDy,1,1)=="X"~substr(ShellIDy,2,10),
                            T~ShellIDy),
         cumGrowth=cumsum(AnGrowth), #cumulative growth along the margin
         # divide individual growth annuli by sum of all growth annuli
         propGrowth=cumGrowth/sum(AnGrowth)*100,
         Age=Year-min(Year)+1) %>%
  left_join(SlideData[,-c(1:6)], by=c("RawName")) %>%
  left_join(ShellDim, by=c("FileName"="Shell.ID")) %>%
  # multiply that proportion by the 'shell extent' to get ~width
  mutate(est.Width=propGrowth*Width.y/100) %>%
  rename(Shell.ID=FileName) %>%
  group_by(Shell.ID) %>%
  dplyr::select(-ShellIDy, -RawName,-Width.x,-curv.width,-ruler.width)

# use shell dim to get regression coefficient to translate width to length
#### anchoring regressions 
anchor<-ShellDim %>% 
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  slice(1) %>%
  dplyr::select(Site.Agg, Species, Length, Width) %>%
  mutate(Length=0,Width=0)

reg<-ShellDim %>% bind_rows(anchor) %>%
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  summarize(Ha=lm(Length~Width)$coefficients[1],
         Hb=lm(Length~Width)$coefficients[2],
         R2=summary(lm(Length~Width))$adj.r.squared*100)

# apply regression to build age x length table
AxL<-AxW %>% 
  left_join(reg, by=c("Site.Agg","Species")) %>%
  mutate(est.Length=Ha+est.Width*Hb,
         Age.adj=Age+as.numeric(start.annuli)) %>%
  dplyr::select(Species,Site,Year, Shell.ID,  AnGrowth,
         Age.adj,Age, age, est.Length,Length, est.Width) %>%
  rename(id=Shell.ID, L=est.Length, Age.unadj=Age,Age=Age.adj)

#check the shells that had eroded umbos (age underestimate)
AxL %>% filter(Age.unadj-Age < 0) %>% slice(1)

# assess difference between final est.Length and real Length
AxL %>% dplyr::select(Species, Site, id,Year, L, Length, age) %>%
  mutate(dif.Length=Length-L) %>%
  slice(n()) %>% arrange(desc(dif.Length))


## Table 1 ------
is_spline<-read_xlsx("data/SplineResults.xlsx", sheet="FINAL") %>%
  rename(Site.Agg=Population)
head(is_spline)
tb1<-AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  dplyr::select(-Site) %>%
  filter(!duplicated(.))%>%
  group_by(Site.Agg, Species, Lat.cor, Long.cor) %>% 
  summarize(shells=length(unique(id)),
            rings=n(),
            maxAge=max(age),
            minAge=min(age)) %>%
  left_join(reg) %>%
  left_join(is_spline)
write.csv(tb1, "figures/table1shell.csv")
min(tb1$R2, na.rm=T) #adj r
min(tb1$Cor, na.rm=T); max(tb1$Cor, na.rm=T) #interseries correlation range

ss_mat<-AxL %>% left_join(SiteID, by=c('Site'='SiteID')) %>%
  group_by(River) %>% arrange(desc(Lat.cor)) %>%
  summarize(Sitesss=paste(unique(Site.Agg), collapse = ', '),
            Spss=paste(unique(Species), collapse = ', '))
write.csv(ss_mat, "figures/table1sites.csv")


#plot for powerpoint
head(AxL)
AxL %>%  select(Site,id) %>%
  slice(1) %>% group_by(Site) %>% tally() %>%
  arrange(desc(n))

AxLplotdf<-AxL %>% filter(Site=="Hudson")

library(cowplot)
alapli<-ggplot()+
  stat_function(data=AxLplotdf[AxLplotdf$Species=="APLI",],
                aes(x=Age, y=L),
              fun=function(x) 148.6*(1-exp(-0.052*(x - 1.12))), 
              size=2)+
  geom_point(data=AxLplotdf[AxLplotdf$Species=="APLI",],
             aes(x=Age, y=L), shape=0, alpha=.5)+
  scale_y_continuous("")+
  scale_x_continuous("", breaks=seq(0,30, by=5))+
  ggtitle(expression(italic("Amblema plicata")))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust=.5),
        axis.title.y = element_text(size=0))
allcar<-ggplot()+
  stat_function(data=AxLplotdf[AxLplotdf$Species=="LCAR",],
                aes(x=Age, y=L),
                fun=function(x) 110.9*(1-exp(-.388*(x - .12))), 
                size=2)+
  geom_point(data=AxLplotdf[AxLplotdf$Species=="LCAR",],
             aes(x=Age, y=L), shape=1, alpha=.5)+
  scale_y_continuous("")+
  scale_x_continuous("", breaks=seq(0,30, by=3))+
  ggtitle(expression(italic("Lampsilis cardium")))+
  theme_cowplot()+
  theme(plot.title = element_text(hjust=.5),
        axis.title.y = element_text(size=0))
plot_grid(alapli, allcar)
