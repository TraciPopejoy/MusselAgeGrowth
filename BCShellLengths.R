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
  mutate(est.Width=propGrowth*Width/100) %>%
  rename(Shell.ID=FileName) %>%
  group_by(Shell.ID) %>%
  select(-ShellIDy, -RawName,-Width,-curv.width,-ruler.width)

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
AxL<-AxW %>% left_join(ShellDim) %>%
  left_join(reg, by=c("Site.Agg","Species")) %>%
  mutate(est.Length=Ha+est.Width*Hb,
         Age.adj=Age+as.numeric(start.annuli)) %>%
  select(Species,Site,Year, Shell.ID,  AnGrowth,
         Age, Age.adj, age, est.Length,Length, est.Width) %>%
  rename(id=Shell.ID, L=est.Length)

#check the shells that had eroded umbos (age underestimate)
AxL %>% filter(Age.adj-Age > 0) %>% slice(1)

# assess difference between final est.Length and real Length
AxL %>% select(Species, Site, id,Year, L, Length, age) %>%
  mutate(dif.Length=Length-L) %>%
  slice(n()) %>% arrange(desc(dif.Length))

APLIaxl<-AxL %>% filter(Species=="APLI")
LAMPaxl<-AxL %>% filter(Species %in% c("LCAR","LORN"))

# recommended MF back calculation model
#Li = x.75Lop + exp(log(Lop-x.75Lop)+
#                (log(Length-x.75Lop)-log(Lop-x.75Lop)*(log(Radi)-log(Radop)))/
#                (log(Radcap)-log(Radop)))

# vonBertanlanffys - Vigliola & Meekan 2009
library(grid);library(lattice);library(stats);library(nlme)

#my data
APdatgr<-groupedData(L~Age|id, data=APLIaxl,
                   labels=list(x="Age",y="Size"),
                   units=list(x="Years",y="mm"))
LCdatgr<-groupedData(L~Age|id, data=LAMPaxl,
                     labels=list(x="Age",y="Size"),
                     units=list(x="Years",y="mm"))

LVB=function(x,t0,Lmax,K){
  y=Lmax*(1-exp(-K*(x-t0)))
  y
}
EXP=function(x,L0,K){
  y=L0*exp(K*x)
  y
}

for(q in 1:4){#length(unique(LCdatgr$Site))){
  subLCdata<- LAMPaxl %>% filter(Site==unique(LCdatgr$Site)[q])
  LCdatgr<-groupedData(L~Age|id, data=subLCdata,
              labels=list(x="Age",y="Size"),
              units=list(x="Years",y="mm"))
  LVB.nlme.LC<-nlme(L~LVB(Age,t0,Lmax,K), data=LCdatgr,
                    fixed=list(t0~1, Lmax~1, K~1),
                    random=t0+Lmax+K~1,
                    start=list(fixed=c(to=-0.02, Lmax=120,K=0.4)))
  
}
EXP.nlme.LC=nlme(L~EXP(age,L0,K), data=LCdatgr,
              fixed=list(L0~1, K~1),
              random=L0+K~1,
              start=list(fixed=c(L0=2,K=0.06)))

summary(LVB.nlme)
intervals(LVB.nlme)
anova(LVB.nlme)
plot(LVB.nlme)
plot(augPred(LVB.nlme, level=0:1))
coef(LVB.nlme)

#### Predicting shell growth with non-linear models ####
library(nlme)
library(FSA)
library(nlstools)

VBcoef<-NULL
agesites<-NULL
agesites<-unique(LAMPaxl$Site)
vbTypical <- L~Linf*(1-exp(-K*(Age-L0)))
for(k in 1:length(agesites)){
  drac<-subset(LAMPaxl, site==as.character(agesites[k]))
  yay<-vbStarts(drac$length~drac$age, 
                param = "vonBertalanffy", 
                methLinf="Walford")
  blo<-data.frame(site=agesites[k], t(unlist(yay)))
  VBcoef<-rbind(VBcoef, blo)
  yay<-NULL
}  

fit87<-nls(vbTypical,data=AgeLengths[AgeLengths$site==agesites[1],], start=VBcoef[1,2:4])
plot(AgeLengths[AgeLengths$site==agesites[1],4], AgeLengths[AgeLengths$site==agesites[1],5],
     xlab="Age", ylab="Shell Length (mm)", main="Lampsilis cardium from Kishwaukee River", pch=19)
lines(seq(1:18), predict(fit87)[200:217], lwd=2)
sum87<-summary(fit87)
sum87