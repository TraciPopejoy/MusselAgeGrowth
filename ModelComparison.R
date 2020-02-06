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
                -temp0001, -Xwtmax, -Xwtmin,-elev)

bf.allsing<-generalTestBF(x50~.*SpN, data=em.data, 
                          whichModels = 'bottom')
    
head(bf.allsing)
tail(bf.allsing)

ds<-head(bf.allsing, n=10)/bf.allsing["SpN"]
plot(ds)

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


