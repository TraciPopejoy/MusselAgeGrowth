library(BayesFactor); library(tidyverse)
library(dataRetrieval)
#load data
source('1-BCShellLengths.R')
site.hyd<-read.csv('data/HITDischarge.csv') %>% dplyr::select(-X) %>%
  mutate(Sp=case_when(Sp %in% c("LCAR")~"LAMP",
               Sp =="APLI"~"APLI")) %>%
  filter(!is.na(Sp))
site.schar<-read.csv('data/SiteChar.csv')%>% dplyr::select(-X)%>%
  filter(!duplicated(.))
site.imp.landuse<-read.csv('data/LanduseT.csv')%>% dplyr::select(-X)
w.sum<-read.csv('data/weatherSum.csv')%>% dplyr::select(-X) 
Site.data<-AxL %>% group_by(Species, Site) %>%
  rename(SiteID=Site)%>%
  summarize(maxYear=max(Year),
            minYear=min(Year)) %>% 
  left_join(SiteID, by="SiteID") %>%
  dplyr::select(Species,Site.Agg, HUC.8,HUC.12,USGS.Gage, maxYear, minYear, Lat.cor, Long.cor) %>%
  filter(!duplicated(Lat.cor))%>%
  rowwise() #%>%
  #mutate(d.av.begin=min(whatNWISdata(siteNumber=USGS.Gage, parameterCd="00060", service="dv")$begin_date),
  #       prob.year=minYear-year(d.av.begin)) %>%
  #arrange(prob.year)

env.model.data<-Site.data %>% ungroup() %>%
  mutate(Sp=case_when(Species %in% c("LCAR")~"LAMP",
                      Species =="APLI"~"APLI"),
         SpF=factor(Sp))%>%
  dplyr::select(Site.Agg, SpF, HUC.8, Lat.cor, Long.cor) %>%
  left_join(lmax.long) %>%
  left_join(site.schar, by=c("Site.Agg")) %>% 
  filter(!duplicated(.),
         Site.Agg != 'Wendell2') %>%
  left_join(site.imp.landuse)%>% 
  left_join(w.sum %>%
                                             mutate(Site.Agg=recode(Site.Agg,'Wendel3'='Wendell3'))) %>%
  left_join(site.hyd%>%
              mutate(Site.Agg=recode(Site.Agg,'Wendel3'='Wendell3')),
            by=c("Site.Agg","SpF"="Sp")) %>%
  replace(is.na(.),0) %>%
  group_by(Site.Agg, SpF) %>%
  mutate(summerLflow=mean(c(ml5,ml6,ml7,ml8)),
         summerHflow=mean(c(mh5,mh6,mh7,mh8)),
         Appt=ppt.cm*totdasqkm)
aptable<-env.model.data %>% ungroup() %>% group_by(Site.Agg)%>%
  select(Site.Agg,  Long.cor, Lat.cor, logDA, HUC.8,
         ppt.cm,Xwtavg,Cultivated.Crops,ml16) %>%
   summarize_if(is.numeric, mean) %>%
  left_join(Site.data[,2:5]) %>% filter(!duplicated(.)) %>%
  arrange(Lat.cor)
                                     
write.csv(aptable, "figures/AppendixTable_102721.csv")

em.data<-env.model.data %>% ungroup() %>%
  mutate(SpN=factor(SpF),
         Lat.round=round(Lat.cor, 2))%>%
  dplyr::select(SpN, x50,Lat.round, 
                ppt.cm, streamorde,logDA,elev, slope,
                Uprop100,Cultivated.Crops,Hay.Pasture,Fprop100,
                Xwtavg,Xwtmax,
                q0001e, v0001e, ma2, ma3, ma36,
                ml5, ml6, ml7, ml8,
                mh5, mh6, mh7, mh8, 
                ml16,ml19, fl1, fh1, dl16, dh15,
                ra8, summerLflow, summerHflow)

bf.allsing<-generalTestBF(x50~., data=em.data, whichRandom = "SpN",
                          whichModels = 'bottom')

head(bf.allsing, n=10)
tail(bf.allsing)

ds<-bf.allsing/bf.allsing["SpN"]
plot(ds)
head(ds, n=10)
tail(bf.allsing, n=10)/bf.allsing["SpN"]
View(em.data)

BF.plotdata<-extractBF(ds) %>% 
  rownames_to_column() %>% 
  dplyr::select(rowname, bf, error) %>%
  mutate(Variable=fct_reorder(rowname, bf),
         VType=case_when(rowname %in% c("Xwtavg","Xwtmax")~'Temperature',
                         rowname %in% c("Lat.round", "ppt.cm",
                                        "streamorde","logDA","elev", "slope")~'Characteristics',
                         rowname %in% c("ra8","ml16","ml19","mh5","mh6","summerHflow","q0001e",
                                        "ml6","mh7","fl1","mh8","ma2", "dh15","ma3","ma36",
                                        'ml5', 'fh1', 'summerLflow', 'ml7', 'ml8','dl16', 'v0001e')~'Flow',
                         rowname %in% c("Uprop100","Cultivated.Crops","Hay.Pasture","Fprop100")~'Land Cover'),
         BFH0SpN=bf-1)
line_table<-data.frame(BF=c('1-3','3-20','20-150','150'),
                       evidence=c('Negligible',
                                  'Positive','Strong',
                                  'Very strong'))
n_cat_breaks<-hist(BF.plotdata$bf,
     breaks = c(0,3,20,60), plot=F)
n_cat_breaks$breaks
n_cat_breaks$counts
ggplot()+
  geom_col(data=BF.plotdata[BF.plotdata$rowname!="SpN",], 
           aes(y=Variable, x=bf, fill=VType))+
  geom_hline(yintercept=(c(26, 26+7, 26+7+2)-.5), size=1.2)+
  annotate(geom="text",x=20,
            y=(c(20, 26+3.5, 26+7+1)-1),
            label=rev(c("Strong Evidence",
                    "Positive Evidence",
                    "Negligible Evidence")),
           color=c("black","black","white"),
           hjust=1, size=3)+
  scale_x_continuous(trans="log1p", expand = c(0,0),
                     breaks=c(-.5,0,1,3,5,10,25, 50, 75, 100),
                     name=expression(paste('Bayes Factor'~~H[0]:L[max]%~%Species)))+
  scale_fill_grey(name="Variable Type")+
  theme_classic()+
  theme(legend.position = c(.8,.175))
ggsave('figures/BayesModelSelect_11.svg', width=3.5, height=6)

chains.temp = posterior(bf.allsing["Xwtavg"], iterations = 10000)
chains.ppt = posterior(bf.allsing["ppt.cm"], iterations = 10000)
chains.crop = posterior(bf.allsing["Fprop100"], iterations = 10000)
chains.dis = posterior(bf.allsing["ml16"], iterations = 10000)
library(bayesplot)
all.c<-bind_rows(mcmc_intervals_data(chains.temp, regex_pars = "Xwtavg"),
                 mcmc_intervals_data(chains.ppt, regex_pars = "ppt.cm"),
                 mcmc_intervals_data(chains.crop, regex_pars = "Fprop100"), 
                 mcmc_intervals_data(chains.dis, regex_pars = "ml16")) %>%
mutate(parF=factor(parameter, levels=c("Fprop100","ml16","ppt.cm","Xwtavg"),
                   labels=c('Streamside Forest',
                            'Minimal Flow\ncoef/10', 'Precipitation',
                            'Avg. Water Temperature')))
all.c[4,5:9]<-all.c[4,5:9]/10

bayesres<-all.c %>%
  ggplot() +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_linerange(aes(ymin = ll, ymax = hh, x = parF))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = parF), size = 2)+ #inner line
  geom_point(aes(y = m, x = parF), size = 3, shape= 21,
             fill="darkgrey")+
  scale_x_discrete(name="")+
  scale_y_continuous("Coefficient")+
  coord_flip()+
  theme_cowplot()+
  theme(axis.title.x = element_text(size=0))
bayesres  
pptgraph<-ggplot(data=env.model.data)+
  geom_smooth(method="lm",aes(x=ppt.cm, y=x50, color=SpF),
              fill=NA)+
  geom_point(aes(x=ppt.cm, y=x50, fill=SpF, shape=SpF),
              color="black", size=2.5)+
  scale_fill_manual(name="Species",values=c("black","white"),
                    labels=c(expression(italic("Amblema plicata")),
                             expression(italic("Lampsilis cardium"))))+
  scale_shape_manual(name="Species",values=c(21,24),
                     labels=c(expression(italic("Amblema plicata")),
                              expression(italic("Lampsilis cardium"))))+
  scale_color_manual(name="Species",values=c("black","grey"),
                     labels=c(expression(italic("Amblema plicata")),
                              expression(italic("Lampsilis cardium"))))+
  labs(y="Potential Max Length (mm)", x="Watershed Precipitation (cm)")+
  theme_classic()+theme(legend.direction="horizontal")
pptgraph

library(ggrepel)
pptwflowgraph<-ggplot(data=env.model.data)+
  geom_smooth(method="lm",aes(x=ppt.cm, y=ma2, color=SpF))+
  geom_point(aes(x=ppt.cm, y=ma2, fill=SpF, shape=SpF),
             color="black", size=2.5)+
  scale_fill_manual(name="Species",values=c("black","white"))+
  scale_shape_manual(name="Species",values=c(21,24))+
  labs(y="Median Annual Flow", x="Watershed Precipitation")+
  scale_y_log10()+
  scale_color_manual(name="Species",values=c("black","white"))+
  geom_text_repel(aes(x=ppt.cm, y=ma2, label=Site.Agg),
                  size=2)+
  theme_classic()+
  theme(legend.position = "none")
pptwflowgraph

Xwtavggraph<-ggplot(data=env.model.data)+
  geom_smooth(method="lm",aes(x=Xwtavg, y=x50, color=SpF),
              fill=NA)+
  geom_point(aes(x=Xwtavg, y=x50, fill=SpF, shape=SpF),
             color="black", size=2.5)+
  scale_color_manual(name="Species",values=c("black","grey"))+
  scale_fill_manual(name="Species",values=c("black","white"))+
  scale_shape_manual(name="Species",values=c(21,24))+
  labs(y="Potential Max Length (mm)", x="Avg. Water Temperature (°C)")+
  theme_classic()+
  theme(legend.position = "none")
Xwtavggraph

minQgraph<-ggplot(data=env.model.data)+
  geom_smooth(method="lm",aes(x=ml16, y=x50, color=SpF),
              fill=NA)+
  geom_point(aes(x=ml16, y=x50, fill=SpF, shape=SpF),
             color="black", size=2.5)+ 
  scale_color_manual(name="Species",values=c("black","grey"))+
  scale_fill_manual(name="Species",values=c("black","white"))+
  scale_shape_manual(name="Species",values=c(21,24))+
  labs(y="Potential Max Length(mm)", 
       x="Minimum Flows")+
  theme_classic()+
  theme(legend.position="none")
minQgraph

ccgraph<-ggplot(data=env.model.data)+
  geom_smooth(method="lm",aes(x=Fprop100, y=x50, color=SpF),
              fill=NA)+
  geom_point(aes(x=Fprop100, y=x50, fill=SpF, shape=SpF),
             color="black", size=2.5)+ 
  scale_color_manual(name="Species",values=c("black","grey"))+
  scale_fill_manual(name="Species",values=c("black","white"))+
  scale_shape_manual(name="Species",values=c(21,24))+
  labs(y="Potential Max Length (mm)", x="Percent Streamside Forest Cover")+
  theme_classic()+
  theme(legend.position="none")
ccgraph

mex<-plot_grid( Xwtavggraph,pptgraph+theme(legend.position = "none"),
               minQgraph, ccgraph, 
          labels=toupper(letters[2:5]))
mexLeg<-get_legend(pptgraph)
plot_grid(bayesres,mexLeg, mex, ncol = 1,
          labels=c("A",NA), rel_heights = c(.375,.06,1))
ggsave("figures/modelexp_111121.tiff", width=6, height=8)

1-ecdf(chains.crop[,2])(0)


lat<-ggplot(data=env.model.data)+
  geom_point(aes(x=ppt.cm, y=Lat.cor), color="red",size=2)+
  theme_classic()
long<-ggplot(data=env.model.data)+
  geom_point(aes(x=ppt.cm, y=Long.cor), color="blue",size=2)+
  theme_classic()
library(cowplot)
plot_grid(lat,long)

em.data.red<-em.data %>%
  dplyr::select(ppt.cm, Xwtavg, elev, streamorde, 
                `Cultivated Crops`, ml16, SpN, logDA, Uprop100, x50)

bf.all<-generalTestBF(x50~., data=em.data.red, noSample=T)
length(bf.all)
head(bf.all)/bf.all["SpN"]
head(bf.all)/bf.all["ppt.cm"]
bf.all["ppt.cm"]/bf.all["logDA"]
bf.all["ppt.cm"]/bf.all["elev + streamorde + logDA"]

bf.latsp<-generalTestBF(x50~Lat.round*SpN,
                        data=em.data)
head(bf.latsp)

chains = posterior(bf.latsp["Lat.round + SpN"], iterations = 10000)
summary(chains)

latsp.c<-mcmc_intervals_data(chains) %>%
  filter(parameter %in% c("Lat.round-Lat.round", "SpN-APLI","SpN-LAMP"))
ggplot(latsp.c) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_linerange(aes(ymin = ll, ymax = hh, x = parameter))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = parameter), size = 2)+ #inner line
  geom_point(aes(y = m, x = parameter), size = 3, shape= 21,
             fill="darkgrey")+
  scale_y_continuous("Coefficient")+
  coord_flip()+
  theme_cowplot()

#appendix graph ----
em.data.reduced<-env.model.data %>% ungroup() %>%
  select(Site.Agg,logDA, streamorde, ppt.cm,Appt, Lat.cor) %>%
  filter(!duplicated(.))
pairs(em.data.reduced[,-c(1,2)])

library(ggrepel)
cor.test(em.data.reduced$logDA, em.data.reduced$streamorde)
daso<-ggplot(data=em.data.reduced, aes(x=logDA, y=streamorde))+
  #geom_smooth(method="lm", color="black")+
  geom_point(size=2.5, alpha=0.8)+
  labs(x="Drainage Area", y="Stream Order")+
  geom_text(x=7.5, y=8, label="r(26) = 0.82, p < 0.001")+
  theme_cowplot()
cor.test(em.data.reduced$logDA, em.data.reduced$ppt.cm)
dapt<-ggplot(data=em.data.reduced,aes(x=logDA, y=ppt.cm))+
  #geom_smooth(method="lm", color="black")+
  geom_point(size=2.5, alpha=0.8)+
  labs(x="Drainage Area", y="Precipitation (cm)")+
  geom_text(x=11, y=137, label="r(26) = -0.58, p = 0.001")+
  theme_cowplot()
cor.test(em.data.reduced$Lat.cor, em.data.reduced$ppt.cm)
latpt<-ggplot(data=em.data.reduced, aes(x=Lat.cor, y=ppt.cm))+
  #geom_smooth(method="lm", color="black")+
  geom_point(size=2.5, alpha=0.8)+
  labs(x="Latitude", y="Precipitation (cm)")+
  geom_text(x=33, y=80, label="r(26) = -0.79, p < 0.001")+
  theme_cowplot()
plot_grid(daso,dapt,latpt, ncol=1, labels="AUTO")
ggsave("figures/Appendixfigure.tiff", width=5, height=8)
