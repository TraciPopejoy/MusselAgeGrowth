library(readxl);library(dplR);library(ggplot2); library(tidyverse)

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
  # multiply that proportion by the 'shell extent' to get ~width
  mutate(straitGrowth=propGrowth*Width/100) %>%
  rename(Shell.ID=FileName)

# use shell dim to get regression coefficient to translate width to length
#### anchoring regressions 
anchor<-ShellDim %>% 
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  slice(1) %>%
  select(Site.Agg, Species, Length, Width) %>%
  mutate(Length=0,Width=0)

reg<-ShellDim %>% bind_rows(anchor) %>%
  group_by(Site.Agg, Species) %>%
  filter(!is.na(Length),!is.na(Width),
         Species %in% c("LORN","LCAR","APLI","QVER"))%>%
  summarize(Ha=lm(Length~Width)$coefficients[1],
         Hb=lm(Length~Width)$coefficients[2])

# apply regression to build age x length table
AxL<-AxW %>% left_join(reg, by=c("Site.Agg","Species")) %>%
  mutate(est.Length=Ha+est.Width*Hb)

# assess difference between final est.Length and real Length
AxL %>% left_join(ShellDim) %>%
  select(Shell.ID, Site.Agg, Species, Length, est.Length) %>%
  mutate(dif.Length=Length-est.Length)

# vonBertanlanffys
VB

