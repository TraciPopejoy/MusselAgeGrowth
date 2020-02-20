library(BayesFactor)
#load data
site.hyd<-read.csv('data/HITDischarge.csv') %>% dplyr::select(-X) %>%
  mutate(Sp=case_when(Sp %in% c("LCAR","LORN")~"LAMP",
               Sp =="APLI"~"APLI"))
site.schar<-read.csv('data/SiteChar.csv')%>% dplyr::select(-X)
site.imp.landuse<-read.csv('data/LanduseT.csv')%>% dplyr::select(-X)
w.sum<-read.csv('data/weatherSum.csv')%>% dplyr::select(-X)

env.model.data<-Site.data %>% ungroup() %>%
  mutate(Sp=case_when(Species %in% c("LCAR","LORN")~"LAMP",
                      Species =="APLI"~"APLI"),
         SpF=factor(Sp))%>%
  dplyr::select(Site.Agg, SpF, HUC.8, Lat.cor, Long.cor) %>%
  full_join(lmax.long) %>%
  left_join(site.schar, by=c("Site.Agg")) %>%
  filter(!duplicated(.)) %>%
  left_join(site.imp.landuse)%>% left_join(w.sum[,-2]) %>%
  filter(!duplicated(.)) %>% 
  left_join(site.hyd, by=c("Site.Agg","SpF"="Sp")) %>%
  replace(is.na(.),0) %>%
  group_by(Site.Agg, SpF) %>%
  mutate(summerLflow=mean(c(ml5,ml6,ml7,ml8)),
         summerHflow=mean(c(mh5,mh6,mh7,mh8)))

em.data<-env.model.data %>% ungroup() %>%
  mutate(SpN=factor(SpF),
         Lat.round=round(Lat.cor, 2))%>%
  dplyr::select(-SpF, -Lat.cor, -mu_l_z, -HUC.8,-gnis_name,-miny,-maxy,
                -Site.Agg, -Long.cor,-totdasqkm, -lengthkm,-nobs, -Latscale,
                -temp0001, -Xwtmin, -ra9)

bf.allsing<-generalTestBF(x50~., data=em.data, whichRandom = "SpN",
                          whichModels = 'bottom')

head(bf.allsing)
tail(bf.allsing)

ds<-head(bf.allsing, n=15)/bf.allsing["SpN"]
plot(ds)
head(ds)
tail(bf.allsing)/bf.allsing["SpN"]


BF.plotdata<-extractBF(ds) %>% 
  rownames_to_column() %>% 
  dplyr::select(rowname, bf, error) %>%
  mutate(Variable=fct_reorder(rowname, bf),
         VType=case_when(rowname %in% c("Xwtavg","Xwtmax")~'Temperature',
                         rowname %in% c("Lat.round", "ppt.cm",
                                        "streamorde","logDA","elev")~'Characteristics',
                         rowname %in% c("ra8","ml16","ml19","mh5","mh6","summerHflow","q0001e",
                                        "ml6","mh7","fl1","mh8","ma2")~'Flow',
                         rowname %in% c("Uprop100")~'Land Cover'),
         BFH0SpN=bf-1)

ggplot()+
  geom_col(data=BF.plotdata[BF.plotdata$rowname!="SpN",], 
           aes(x=Variable, y=BFH0SpN, fill=VType))+
  geom_vline(xintercept=c(13.5,12.5,9.5,6.5), size=1.4)+
  annotate(geom="text",y=c(4,4,27,27),
            x=c(14,13,10,7),
            label=c("Very Strong Evidence",
                    "Strong Evidence",
                    "Moderate Evidence",
                    "Anecdotal Evidence"),
           color=c("white","white","black","black"))+
  scale_y_continuous(trans="log1p", 
                     breaks=c(-.5,0,1,3,10,30,60),
                     name=expression(paste('Bayes Factor'~~H[0]:L[max]%~%Species)))+
  scale_fill_grey(name="Variable Type")+
  theme_classic()+
  theme(legend.position = c(.8,.175))+
  coord_flip()
  
ggplot(data=env.model.data)+
  geom_point( aes(x=ppt.cm, y=x50, fill=SpF, shape=SpF),
              color="black", size=2)+
  scale_fill_manual(name="Species",values=c("black","white"))+
  scale_shape_manual(name="Species",values=c(21,24))+
  theme_classic()+
  theme(legend.position = c(.9,.9))


lat<-ggplot(data=env.model.data)+
  geom_point(aes(x=ppt.cm, y=Lat.cor), color="red",size=2)+
  theme_classic()
long<-ggplot(data=env.model.data)+
  geom_point(aes(x=ppt.cm, y=Long.cor), color="blue",size=2)+
  theme_classic()
library(cowplot)
plot_grid(lat,long)


# Model comparison ------------
bf<-generalTestBF(x50~., data=em.data,
                  progress=T, neverExclude = "SpN") #341
length(bf)
head(bf, n=10)

head(bf)/bf["SpN"]


me.hyp<-head(bf)/bf["Lat.cor"]
plot(me.hyp)

bf.sp<- regressionBF(x50~ . , data=em.data[,1:4])

bf.Lat <- regressionBF(x50 ~ Lat.round,  data=env.model.data)
bf.Sp <- regressionBF(x50 ~ SpF, data = env.model.data)
bf.streamorde <- regressionBF(x50 ~ streamorde, data = env.model.data)
bf.Xwtavg <- regressionBF(x50 ~ Xwtavg, data = env.model.data)
bf.Xwtmax <- regressionBF(x50 ~ Xwtmax, data = env.model.data)
bf.summerLflow <- regressionBF(x50 ~ summerLflow, data = env.model.data)
bf.summerHflow <- regressionBF(x50 ~ summerHflow, data = env.model.data)

bf.all<-c(bf.Lat,bf.Sp,bf.streamorde,bf.Xwtavg,bf.Xwtmax,
          bf.summerLflow, bf.summerHflow)
length(bf.all)
bf.all['Lat.round']
head(bf.all)
## Compare the 5 best models to the best
bf2 = head(bf.all) / max(bf.all)
bf2
plot(bf2)

chains = posterior(complaintsLearningBf, iterations = 10000)
summary(chains)


