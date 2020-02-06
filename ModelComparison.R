library(BayesFactor)
#load data
site.hyd<-read.csv('data/HITDischarge.csv') %>% dplyr::select(-X) %>%
  mutate(Sp=case_when(Sp %in% c("LCAR","LORN")~"LAMP",
               Sp =="APLI"~"APLI"))
site.schar<-read.csv('data/SiteChar.csv')%>% dplyr::select(-X)
site.imp.landuse<-read.csv('data/LanduseT.csv')
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
  dplyr::select(x50,SpF,Lat.cor,streamorde, Xwtavg, Xwtmax,
                ma2, fl1,fh1,dl16,dl17,ra8,summerLflow,summerHflow)

em.data<-as.data.frame(em.data)
# Model comparison ------------
cor.test(y=em.data$x50, x=em.data$Lat.cor)
correlationBF(y = em.data$x50, x = em.data$Lat.cor)

bf.Lat <- regressionBF(x50 ~ Lat.cor, data = em.data)
bf.Sp <- regressionBF(x50 ~ SpF, data = em.data)
bf.streamorde <- regressionBF(x50 ~ streamorde, data = em.data)
bf.Xwtavg <- regressionBF(x50 ~ Xwtavg, data = em.data)
bf.Xwtmax <- regressionBF(x50 ~ Xwtmax, data = em.data)
bf.summerLflow <- regressionBF(x50 ~ summerLflow, data = em.data)
bf.summerHflow <- regressionBF(x50 ~ summerHflow, data = em.data)

bf.all<-c(bf.Lat,bf.Sp,bf.streamorde,bf.Xwtavg,bf.Xwtmax,
          bf.summerLflow, bf.summerHflow)
length(bf.all)
bf.all['Lat.cor']
head(bf.all)
## Compare the 5 best models to the best
bf2 = head(bf.all) / max(bf.all)
bf2
plot(bf2)


chains = posterior(complaintsLearningBf, iterations = 10000)
summary(chains)


