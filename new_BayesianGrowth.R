# Bayesian von Bertanlanffy's
set.seed(6363) #set random number seed for replicable analysis

### testing priors 
#hist(runif(1000, 10,200)) # prior on Lmax
#hist(runif(1000,0,1)) # prior on K
#hist(rnorm(1000, 0,.1)) # prior on T0
#hist(rgamma(1000,.1,.1)) #variance on L
#hist(rgamma(1000,1,1)) #variance on mu_l
#hist(rgamma(1500,.001,.001),breaks=1300,xlim=c(0,1)) #variance on tauk
#hist(rgamma(1000,.01,.01),breaks=800,xlim=c(0,1)) #variance on tau_t0

# Von Bertanlanffy's model
# L = Lmax*(1-exp(-K*(Age-t0)))

library(rjags);library(MCMCpack)

nadapt=5000
burnin=10000
nit=70000
thin=4

##### Model Code #####
model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
L[i] ~ dnorm(mu[i],tau) 
mu[i] <- Lmax[idF[i]] * (1-exp(-K[idF[i]] * (Age[i] - T0[idF[i]])))
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L

for(u in 1:nid){
# so Lmax, K and T0 will be unique for each individual
Lmax[u]~ dnorm(mu_l, tau_l)
K[u] ~ dnorm(mu_k, tau_k)
T0[u] ~ dnorm(mu_t0, tau_t0)
}

# hyper parameter priors
mu_l ~ dunif(60, 160)
tau_l ~ dgamma(.01,.01)

mu_k ~ dunif(0.01,1)
tau_k ~ dgamma(.0001,.0001)

mu_t0 ~ dnorm (0, .1)
tau_t0 ~ dgamma(.0001,.0001)
}"
v<-c("mu_l","mu_k","mu_t0")

#getting a vector of sites so my factor levels line up
Lsites<- AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species %in% c("LCAR","LORN")) %>%
  pull(Site.Agg) %>% unique()

lamp.mcmc.data<-rep(NA,nit/thin*3) #start the final mcmc data set
lamp.mcmc.sum<-NULL #start the final mcmc summary stats
lamp.gr<-NULL

#loop through each site, pull data, run model, save mcmc
for(s in 1:length(Lsites)) {
  
  data<-AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species %in% c("LCAR","LORN")) %>%
  filter(Site.Agg==Lsites[s])%>%
  mutate(idF=as.factor(id)) %>%
  dplyr::select(idF,Age,L, Latitude)
  
  nid<-length(unique(data$idF))
  nobs<-nrow(data)
  
  model<-jags.model(textConnection(model_string), 
                      data=list(L = data$L,
                                Age = data$Age,
                                idF = data$idF,
                                nid=nid,
                                nobs=nobs),
                      inits=list(Lmax=rep(100,nid),
                                 K=rep(.5,nid),
                                 T0=rep(0,nid)),
                      n.chains=3, n.adapt=nadapt)
  update(model, burnin) # burn in for 200000 samples
  mcmc<-coda.samples(model,
                      variable.names=c("mu_l","mu_k","mu_t0"), 
                      n.iter=nit, thin=thin)
  #save the summary statistics
  mcmc.sum<-data.frame(variable=rownames(summary(mcmc)$statistics),
                       mean=summary(mcmc)$statistics[,1],
                       SD=summary(mcmc)$statistics[,2],
                       SiteID=Lsites[s],
                       x2.5=summary(mcmc)$quantiles[,1],
                       x25=summary(mcmc)$quantiles[,2],
                       x50=summary(mcmc)$quantiles[,3], 
                       x75=summary(mcmc)$quantiles[,4],
                       x97.5=summary(mcmc)$quantiles[,5])
  lamp.mcmc.sum<-rbind(lamp.mcmc.sum, mcmc.sum)
  
  #save my mcmc chains for later use
  mcmc.data<-rbind(as.matrix(mcmc[1]),as.matrix(mcmc[2]),
                   as.matrix(mcmc[3]))
  names(v)<-c(paste("mu_l",s, sep="_"),
                 paste("mu_k",s, sep="_"),
                 paste("mu_t0",s, sep="_"))
  mcmc.data <- as.data.frame(mcmc.data) %>% 
    dplyr::rename(all_of(v))
  lamp.mcmc.data<-cbind(lamp.mcmc.data, mcmc.data)
  #print the mcmc diagnostics to check
  gr<-cbind(gelman.diag(mcmc)[[1]],gelman.diag(mcmc)[[2]], Lsites[s])
  lamp.gr<-rbind(lamp.gr, gr)
  
  pdf(paste(paste('mcmc_output/LampMcmcDiag_Lmax',s, sep="_"),'.pdf',sep=""))
  plot(mcmc)
  gelman.plot(mcmc)
  dev.off()
}

names(lamp.mcmc.data)
library(bayesplot)
mcmc_intervals(lamp.mcmc.data[,-1], regex_pars = "mu_l") +
  coord_flip()
head(lamp.mcmc.sum)
lamp.lmax.graph<-lamp.mcmc.sum %>% left_join(SiteID)
ggplot(data=lamp.lmax.graph)+
  geom_crossbar(aes(x=Latitude,y=x50,ymin=x2.5, ymax=x97.5), linetype="dashed")+
  geom_crossbar(aes(x=Latitude,y=x50,ymin=x25, ymax=x75))+
  geom_point(aes(x=Latitude, y=mean),size=2)+
  facet_wrap(~variable, scales="free_y")+
  geom_smooth(aes(x=Latitude, y=mean),method="lm", se=F)+
  theme_bw()

write.csv(lamp.mcmc.data, "Lamp_Lmax_mcmcres.csv")
write.csv(lamp.mcmc.sum, "Lamp_Lmax_mcmc_sum.csv")
View(lamp.gr)

# Amblema -------
#getting a vector of sites so my factor levels line up
Asites<- AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species=="APLI") %>%
  pull(Site.Agg) %>% unique()

amb.mcmc.data<-rep(NA,nit/thin*3) #start the final mcmc data set
amb.mcmc.sum<-NULL #start the final mcmc summary stats
amb.gr<-NULL


#loop through each site, pull data, run model, save mcmc
for(s in 1:length(Asites)) {
  
  data<-AxL %>% ungroup() %>%
    left_join(SiteID, by=c('Site'='SiteID')) %>%
    filter(Species=="APLI") %>%
    filter(Site.Agg==Asites[s])%>%
    mutate(idF=as.factor(id)) %>%
    dplyr::select(idF,Age,L, Latitude)
  
  nid<-length(unique(data$idF))
  nobs<-nrow(data)
  
  model<-jags.model(textConnection(model_string), 
                    data=list(L = data$L,
                              Age = data$Age,
                              idF = data$idF,
                              nid=nid,
                              nobs=nobs),
                    inits=list(Lmax=rep(100,nid),
                               K=rep(.5,nid),
                               T0=rep(0,nid)),
                    n.chains=3, n.adapt=nadapt)
  update(model, burnin) # burn in for 200000 samples
  mcmc<-coda.samples(model,
                     variable.names=c("mu_l","mu_k","mu_t0"), 
                     n.iter=nit, thin=thin)
  #save the summary statistics
  mcmc.sum<-data.frame(variable=rownames(summary(mcmc)$statistics),
                       mean=summary(mcmc)$statistics[,1],
                       SD=summary(mcmc)$statistics[,2],
                       SiteID=Asites[s],
                       x2.5=summary(mcmc)$quantiles[,1],
                       x97.5=summary(mcmc)$quantiles[,5],
                       x50=summary(mcmc)$quantiles[,3]) 
  amb.mcmc.sum<-rbind(amb.mcmc.sum, mcmc.sum)
  
  #save my mcmc chains for later use
  mcmc.data<-rbind(as.matrix(mcmc[1]),as.matrix(mcmc[2]),
                   as.matrix(mcmc[3]))
  names(v)<-c(paste("mu_l",s, sep="_"),
              paste("mu_k",s, sep="_"),
              paste("mu_t0",s, sep="_"))
  mcmc.data <- as.data.frame(mcmc.data) %>% 
    dplyr::rename(all_of(v))
  amb.mcmc.data<-cbind(amb.mcmc.data, mcmc.data)
  #print the mcmc diagnostics to check
  #print the mcmc diagnostics to check
  gr<-cbind(gelman.diag(mcmc)[[1]],gelman.diag(mcmc)[[2]], Asites[s])
  amb.gr<-rbind(amb.gr, gr)
  
  pdf(paste(paste('mcmc_output/AMBMcmcDiag_Lmax',s, sep="_"),'.pdf',sep=""))
  plot(mcmc)
  gelman.plot(mcmc)
  dev.off()
}
names(amb.mcmc.data)
mcmc_intervals(amb.mcmc.data[,-1], regex_pars = "mu_l")+
  scale_y_discrete(labels=Asites, name="Sites")+
  coord_flip()
write.csv(amb.mcmc.data, "Amb_Lmax_mcmcres.csv")
write.csv(amb.mcmc.sum, "Amb_Lmax_mcmc_sum.csv")
View(amb.gr)


amb.mcmc.data<-read.csv("Amb_Lmax_mcmcres.csv") %>% 
  dplyr::select(-X)
lamp.mcmc.data<-read.csv("Lamp_Lmax_mcmcres.csv") %>% 
  dplyr::select(-X)
# plots ------
color_scheme<-color_scheme_get("gray")
#lampsilis
lamp.lmax.graph1<-mcmc_intervals_data(lamp.mcmc.data[,-1], regex_pars = "mu_l") %>%
  mutate(Site.Agg=Lsites[as.numeric(substr(parameter,6,7))]) %>%
  left_join(SiteID[,c(2,13)]) %>% filter(!duplicated(.)) %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
llmax.con<-ggplot(lamp.lmax.graph1) +
  geom_ribbon(aes(x=Lat.cor,ymin=2.6003*Lat.cor, ymax=3.0162*Lat.cor),
              fill="lightgrey", alpha=.5)+
  geom_linerange(aes(ymin = ll, ymax = hh, x = Lat.cor))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = Lat.cor), size = 2)+ #inner line
  geom_point(aes(y = m, x = Lat.cor), size = 3, shape= 21,
             fill="darkgrey")+
  scale_x_continuous(name="Latitude", limits=c(29.33, 45.4))+
  scale_y_continuous(name="Maximum Length (mm)", limits = c(60,160))+
  geom_abline(intercept=0,slope=2.7917)+
  theme_classic()+
  theme(axis.title.x = element_text(size=0))
llmax.dis<-ggplot(lamp.lmax.graph1) +
  geom_segment(aes(x = ll, xend = hh, y = SiteLat, yend=SiteLat),
               color=color_scheme[[3]])+  #outer line
  geom_segment(aes(x = l, xend = h, y = SiteLat, yend = SiteLat),
               size = 2, color=color_scheme[[5]])+ #inner line
  geom_point(aes(x = m, y = SiteLat), size = 4, shape= 21,
             color=color_scheme[[6]],fill=color_scheme[[1]])+
  scale_x_continuous(name="Maximum Length (mm)")+
  scale_y_discrete(name="Sites")+
  theme_classic()+coord_flip()+
  theme(axis.text.x = element_text(angle = 40, hjust=.9),
        axis.title.y = element_text(size=0))
#amblema
amb.lmax.graph1<-mcmc_intervals_data(amb.mcmc.data[,-1], regex_pars = "mu_l") %>%
  mutate(Site.Agg=Asites[as.numeric(substr(parameter,6,7))]) %>%
  left_join(SiteID[,c(2,13)]) %>%
  filter(!duplicated(.)) %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
amb.lmax.graph1$image <- "figures/amblemaP.png"
library(ggimage)
#abiorect<-data.frame(xmin=c(29.3,31.5,33.5),
#                     xmax=c(31.5,33.5,37), 
#                     ymin=rep(59,3), ymax=rep(62,3))
almax.con<-ggplot(amb.lmax.graph1) +
  #geom_rect(data=abiorect, aes(xmin=xmin,xmax=xmax,
  #                               ymin=ymin,ymax=ymax), 
  #          fill=wes_palette("Darjeeling1")[c(1,4,3)])+
  #geom_image(aes(image = image), x=30.5, y=153, size=.15)+
  #geom_text(y=145, x=31, label=expression(italic("A. plicata")))+
  geom_ribbon(aes(x=Lat.cor,ymin=2.95*Lat.cor, ymax=3.34*Lat.cor),
              fill="lightgrey", alpha=.5)+
  geom_linerange(aes(ymin = ll, ymax = hh, x = Lat.cor))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = Lat.cor),
               size = 2)+ #inner line
  geom_point(aes(y = m, x = Lat.cor), size = 3, shape= 21, fill="white")+
  scale_x_continuous(name="Latitude")+
  scale_y_continuous(name="Maximum Length (mm)")+
  geom_abline(intercept=0,slope=3.135)+
  theme_classic()
almax.con

almax.dis<-ggplot(amb.lmax.graph1) +
  geom_segment(aes(x = ll, xend = hh, y = SiteLat, yend=SiteLat),
               color=color_scheme[[3]])+  #outer line
  geom_segment(aes(x = l, xend = h, y = SiteLat, yend = SiteLat),
               size = 2, color=color_scheme[[5]])+ #inner line
  geom_point(aes(x = m, y = SiteLat), size = 4, shape= 21,
             color=color_scheme[[6]],fill=color_scheme[[1]])+
  scale_x_continuous(name="Maximum Length (mm)")+
  scale_y_discrete(name="Sites")+
  theme_classic()+coord_flip()+
  theme(axis.text.x = element_text(angle = 40, hjust=.9))

library(cowplot)
plot_grid(llmax.con, almax.con, ncol=1, labels="AUTO")
ggsave("figures/LmaxContinuous.tiff", width=6.5, height=7)

plot_grid(llmax.dis, almax.dis, ncol=1, labels="AUTO")
ggsave("figures/LmaxDiscrete.tiff", width=6.5, height=7)


lmax.graph<-bind_rows(amb.lmax.graph1 %>% mutate(Species="A. plicata"),
          lamp.lmax.graph1 %>% mutate(Species="Lampsilis spp.")) %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
bioregrect<-data.frame(xmin=c(0,2.5,6.5,12.5,24.5,22.5),
                       xmax=c(2.5,6.5,12.5,28.5,25.5,23.5), 
                       ymin=rep(59,6), ymax=rep(62,6))
lmax.dis<-ggplot(lmax.graph) +
  geom_rect(data=bioregrect, aes(xmin=xmin,xmax=xmax,
                                 ymin=ymin,ymax=ymax), 
            fill=wes_palette("Darjeeling1")[c(1,4,3,5,2,2)])+
  geom_linerange(aes(ymin = ll, ymax = hh, x = SiteLat,
                     group=Species), position=position_dodge(.75))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = SiteLat, group=Species),
               size = 2,position=position_dodge(.75))+ #inner line
  geom_point(aes(x = SiteLat, y = m, fill=Species, group=Species), 
             size = 3, shape= 21, 
             position=position_dodge(.75))+
  scale_y_continuous(name="Maximum Length (mm)")+
  scale_x_discrete(name="Sites")+
  scale_fill_manual(values=c("white","darkgrey"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, hjust=1),
        legend.position = 'top')
lmax.dis
plot_grid(plot_grid(llmax.con, almax.con, ncol=1, 
                    labels=c("A","B")),
          lmax.dis, nrow=1, labels=c("","C"),
          rel_widths = c(.6,1.1))
ggsave('figures/lmaxbio3.tiff', width=9, height=5)

# results numbers ------
amb.mcmc.sum %>% filter(variable=="mu_l") %>% summarise(mean(x50), min(x50), max(x50))
amb.mcmc.sum %>% filter(variable=="mu_k") %>% summarise(mean(x50), min(x50), max(x50))
amb.mcmc.sum %>% filter(variable=="mu_t0") %>% summarise(mean(x50), min(x50), max(x50))

lamp.mcmc.sum %>% filter(variable=="mu_l") %>% summarise(mean(x50), min(x50), max(x50))
lamp.mcmc.sum %>% filter(variable=="mu_k") %>% summarise(mean(x50), min(x50), max(x50))
lamp.mcmc.sum %>% filter(variable=="mu_t0") %>% summarise(mean(x50), min(x50), max(x50))

# compare between sites with both species
ggplot(lmax.graph[lmax.graph$Site.Agg %in% c("K2","KT","K4","KS","K7",
                                             "MS1","Upstream","Sylvan",
                                             "STF","Hudson",
                                             "Barry","Wendel3","MusselMania","Wendel2"),]) +
  geom_linerange(aes(ymin = ll, ymax = hh, x = SiteLat,
                     group=Species), position=position_dodge(.75))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = SiteLat, group=Species),
                 size = 2,position=position_dodge(.75))+ #inner line
  geom_point(aes(x = SiteLat, y = m, fill=Species, group=Species), 
             size = 3, shape= 21, 
             position=position_dodge(.75))+
  #geom_text(aes(x = SiteLat, y = 60, label=seq(1:40), group=Species),
  #          position=position_dodge(.75), size=3)+
  scale_y_continuous(name="Maximum Length (mm)")+
  scale_x_discrete(name="Sites")+
  scale_fill_manual(values=c("white","darkgrey"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, hjust=.9),
        legend.position = 'top')
