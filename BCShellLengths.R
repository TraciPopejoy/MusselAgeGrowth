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
  left_join(ShellDim, by=c("ShellIDy"="Shell.ID")) %>%
  # multiply that proportion by the 'shell extent' to get ~width
  mutate(est.Width=propGrowth*Width.y/100) %>%
  rename(Shell.ID=FileName) %>%
  group_by(Shell.ID) %>%
  select(-ShellIDy, -RawName,-Width.x,-curv.width,-ruler.width)

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
AxL<-AxW %>% 
  left_join(reg, by=c("Site.Agg","Species")) %>%
  mutate(est.Length=Ha+est.Width*Hb,
         Age.adj=Age+as.numeric(start.annuli)) %>%
  select(Species,Site,Year, Shell.ID,  AnGrowth,
         Age.adj,Age, age, est.Length,Length, est.Width) %>%
  rename(id=Shell.ID, L=est.Length, Age.unadj=Age,Age=Age.adj)

#check the shells that had eroded umbos (age underestimate)
AxL %>% filter(Age.unadj-Age < 0) %>% slice(1)

# assess difference between final est.Length and real Length
AxL %>% select(Species, Site, id,Year, L, Length, age) %>%
  mutate(dif.Length=Length-L) %>%
  slice(n()) %>% arrange(desc(dif.Length))

APLIaxl<-AxL %>% filter(Species=="APLI")
LAMPaxl<-AxL %>% filter(Species %in% c("LCAR","LORN")) %>%
   #left_join(SiteID, by=c("Site"="SiteID")) %>% 
  select(Site, Species, id, Age, L)

# Bayesian von Bertanlanffy's
library(brms)
set.seed(6363) #set random number seed for replicable analysis

### testing priors 
hist(rstudent_t(1000, 3,0,10)) #default prior on sd(TankIntercept)
hist(rgamma(1000,0.01,0.01), xlim=c(0,70)) # prior on shape of Treat:Day relationship
hist(rnorm(1000, 0,1)) # prior on beta (Treat:Day interaction)
hist(rnorm(1000, 100,50)) # prior on beta (Treat:Day interaction)
hist(rnorm(1000, .5,.001)) # prior on beta (Treat:Day interaction)


# Von Bertanlanffy's model
# L = Lmax*(1-exp(-K*(Age-t0)))
LampD<-LAMPaxl %>% left_join(SiteID, by=c('Site'='SiteID')) %>%
  dplyr::select(Site, id, Age, L, Latitude)

##### Model Code #####
model_string<- "model {
  # likelihood
  for (j in 1:length(id)) {  
    # varying Lmax/K/t0 for every site
      L[j] ~ dnorm(mu[j],tau) 
      mu[j]<-Lmax*(1-exp(-K*(Age[j]-T0)))
      #mu[j]<-Lmax[Site[m]] * (1-exp(-K[Site[m]]*(Age[j]-T0[Site[m]])))
      #Lmax[Site[j]]<-Latitude-1
  }
  # priors
  tau ~ dgamma(0.001, 0.001) #variance on L
  Lmax ~ dnorm(100, 50)
  K ~ dnorm(0.4, 0.5)
  T0 ~ dnorm(0, 0.01)
  }"


inits<-list(Lmax=100,
            K=.5,
            T0=.1)

#run the model
library(rjags);library(MCMCpack)
NHmodel<-jags.model(textConnection(model_string), 
                    data=LampD, 
                    inits=inits,
                    n.chains=3, n.adapt=30000)
update(NHmodel, 10000) # burn in for 2000 samples
NHmcmc<-coda.samples(NHmodel,
                     variable.names=c("Lmax","K", "T0"), 
                     n.iter=5000, thin=40)
pdf('mcmcDiagnostic/nh4mcmcdiag.pdf')
plot(NHmcmc)
gelman.plot(NHmcmc)
dev.off()

summary(NHmcmc)
NH4mcmc.data<-as.matrix(NHmcmc); colnames(NH4mcmc.data)




# quantifying the % probability response increased after impact
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1 ) 
nhBACIgraph %>% group_by(ratio) %>% dplyr::summarize(meanBACI=mean(value))
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1.5 ) 
marginal_effects(nhBmodel)$`TreF:DayF` %>% 
  select(TreF, DayF, estimate__, lower__,upper__) %>%
  filter(DayF==4)




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



LCAR.nlme.res<-list()
LCAR.nlme.sum<-NULL
for(q in 2:length(unique(LAMPaxl$Site))){
  subLCdata<- LAMPaxl %>% filter(Site==unique(LAMPaxl$Site)[q])
  LCdatgr<-groupedData(L~Age|id, data=subLCdata,
              labels=list(x="Age",y="Size"),
              units=list(x="Years",y="mm"))
  
  LVB.nlme.LC<-nlme(L~LVB(Age,t0,Lmax,K), data=LCdatgr,
                    fixed=list(t0~1, Lmax~1, K~1),
                    random=t0+Lmax+K~1,
                    start=list(fixed=c(t0=.2, Lmax=100,K=0.45)),
                    control=nlmeControl(maxIter = 200, msMaxIter = 200))
  LCAR.nlme.res[[unique(LAMPaxl$Site)[q]]] <- LVB.nlme.LC
  lcv<-data.frame(Site=unique(LAMPaxl$Site)[q],
                  t0=LVB.nlme.LC$coefficients$fixed[1], 
                  Lmax=LVB.nlme.LC$coefficients$fixed[2],
                  K=LVB.nlme.LC$coefficients$fixed[3])
  LCAR.nlme.sum<-rbind(LCAR.nlme.sum, lcv)
}


EXP.nlme.LC=nlme(L~EXP(age,L0,K), data=LCdatgr,
              fixed=list(L0~1, K~1),
              random=L0+K~1,
              start=list(fixed=c(L0=2,K=0.06)))

summary(LVB.nlme.LC)
intervals(LVB.nlme)
anova(LVB.nlme.LC)
plot(LVB.nlme.LC)
plot(augPred(LVB.nlme.LC, level=0:1))
coef(LVB.nlme.LC)