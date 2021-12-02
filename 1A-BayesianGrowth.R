# Bayesian von Bertanlanffy's
source('1-BCShellLengths.R')
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

library(tidyverse);library(rjags)

nadapt=6000
burnin=11000
nit=100000
thin=8

##### Model Code #####
#sex differences
model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
L[i] ~ dnorm(mu[i],tau) 
mu[i] <- Lmax[idF[i],sex[i]] * (1-exp(-K[idF[i],sex[i]] * (Age[i] - T0[idF[i]])))
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L


for(u in 1:nid){
for(s in 1:ns){
# so Lmax, K and T0 will be unique for each individual
Lmax[u,s]~ dnorm(mu_l[s], tau_l)T(0,1e7)
K[u,s] ~ dnorm(mu_k[s], tau_k)T(0,1e7)
}
T0[u] ~ dnorm(mu_t0, tau_t0)T(0,1e7)
}

# hyper parameter priors
# mean lmax and k are allowed to vary by sex
for(n in 1:ns){
mu_l[n] ~ dnorm(h_mu_l, h_tau_l)T(0,1e6)
mu_k[n] ~ dnorm(h_mu_k, h_tau_k)T(0,1e6)
}

h_mu_l ~ dunif(60,200)
h_tau_l ~dgamma(.01,.01)

h_mu_k ~ dunif(0.001,1)
h_tau_k ~ dgamma(.0001,.0001)

tau_l ~ dgamma(.01,.01)
tau_k ~ dgamma(.0001,.0001)

mu_t0 ~ dnorm (0, .1)
tau_t0 ~ dgamma(.0001,.0001)
}"
v<-c("h_mu_l", "h_mu_k",
     "mu_l[1]","mu_k[1]","mu_l[2]","mu_k[2]","mu_l[3]","mu_k[3]",
     'mu_l','mu_k',"mu_t0")

#getting a vector of sites so my factor levels line up
Lsites<- AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species %in% c("LCAR")) %>%
  pull(Site.Agg) %>% unique()

lamp.mcmc.data<-NULL #start the final mcmc data set
lamp.mcmc.sum<-NULL #start the final mcmc summary stats
lamp.gr<-NULL

# have to set up a dataframe to identify pars not to keep
all_sex_vars<-data.frame(name=c(rep('par_lmax',4), rep('par_k',4)),
           value=c('mu_l','mu_l[1]','mu_l[2]','mu_l[3]',
                   'mu_k','mu_k[1]','mu_k[2]','mu_k[3]'))

#loop through each site, pull data, run model, save mcmc
for(s in 1:length(Lsites)) {
  
  data<-AxL %>% ungroup() %>%
  left_join(SiteID, by=c('Site'='SiteID')) %>%
  filter(Species %in% c("LCAR")) %>%
  filter(Site.Agg==Lsites[s])%>%
  mutate(idF=as.factor(id), sexF=as.factor(Sex)) %>%
  dplyr::select(idF,Age,L, Latitude, sexF)
  
  nid<-length(unique(data$idF))
  ns<-length(unique(data$sexF))
  nobs<-nrow(data)
  
  model<-jags.model(textConnection(model_string), 
                      data=list(L = data$L,
                                Age = data$Age,
                                idF = data$idF,
                                sex = data$sexF,
                                ns=ns,
                                nid=nid,
                                nobs=nobs),
                     # inits=list(Lmax=rep(100,nid),
                     #            K=rep(.5,nid),
                     #            T0=rep(0,nid)),
                      n.chains=3, n.adapt=nadapt)
  update(model, burnin) # burn in for 200000 samples
  mcmc<-coda.samples(model,
                      variable.names=c("h_mu_l", "h_mu_k",
                                       "mu_l","mu_k","mu_t0"),
                                       #"Lmax","K"), 
                      n.iter=nit, thin=thin)
  #gelman.diag(mcmc)
  if(length(unique(data$sexF))==1){
    sex_vars_real<-data.frame(fact_level=as.numeric(unique(data$sexF)),
                              sex_level=as.character(unique(data$sexF))) %>%
      mutate(par_lmax='mu_l',
             par_k='mu_k') %>%
      pivot_longer(cols=starts_with('par'))
  }else{
  sex_vars_real<-data.frame(fact_level=as.numeric(unique(data$sexF)),
            sex_level=as.character(unique(data$sexF))) %>%
    mutate(par_lmax=paste0('mu_l[',fact_level,']'),
           par_k=paste0('mu_k[', fact_level,']')) %>%
    pivot_longer(cols=starts_with('par'))
  }
  dest_sex<-all_sex_vars %>% anti_join(sex_vars_real, by=c('name','value'))
  #save the summary statistics
  mcs<-summary(mcmc)
  mcmc.sum<-data.frame(variable=rownames(mcs$statistics),
                       mean=mcs$statistics[,1],
                       SD=mcs$statistics[,2],
                       SiteID=Lsites[s],
                       x2.5=mcs$quantiles[,1],
                       x25=mcs$quantiles[,2],
                       x50=mcs$quantiles[,3], 
                       x75=mcs$quantiles[,4],
                       x97.5=mcs$quantiles[,5],
                       s=s) %>%
    #remove parameters that could not be informed by the data
    filter(!(variable %in% c(dest_sex$value))) %>%
    left_join(sex_vars_real, by=c('variable'='value'))
  lamp.mcmc.sum<-rbind(lamp.mcmc.sum, mcmc.sum)
  
  #save my mcmc chains for later use
  mcmc.data<-rbind(as.matrix(mcmc[1][,which(colnames(mcmc[[1]]) %in% mcmc.sum$variable)]),
                   as.matrix(mcmc[2][,which(colnames(mcmc[[1]]) %in% mcmc.sum$variable)]),
                   as.matrix(mcmc[3][,which(colnames(mcmc[[1]]) %in% mcmc.sum$variable)]))
  names(v)<-c(paste0("h_mu_l_", s), paste0("h_mu_k_",s),
              paste0("mu_l[1]_",s), paste0("mu_k[1]_",s),
              paste0("mu_l[2]_",s), paste0("mu_k[2]_",s),
              paste0("mu_l[3]_",s), paste0("mu_k[3]_",s),
              paste0("mu_l_",s), paste0("mu_k_",s),
                 paste0("mu_t0_",s))
  v_x<-v[v %in% mcmc.sum$variable]
  mcmc.data <- as.data.frame(mcmc.data) %>% 
    dplyr::rename(all_of(v_x))
  lamp.mcmc.data<-bind_cols(lamp.mcmc.data, mcmc.data)
  #print the mcmc diagnostics to check
  gr_d<-gelman.diag(mcmc)
  gr<-cbind(gr_d[[1]],gr_d[[2]], Lsites[s])
  lamp.gr<-rbind(lamp.gr, gr)
  
  pdf(paste(paste('mcmc_output/LampMcmcDiag_Lmax',s, sep="_"),'.pdf',sep=""))
  plot(mcmc)
  gelman.plot(mcmc)
  dev.off()
  print(paste('completed', s, 'of', length(Lsites)))
}

names(lamp.mcmc.data)
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
Lmax[u]~ dnorm(mu_l, tau_l)T(0,1e7)
K[u] ~ dnorm(mu_k, tau_k)T(0,1e7)
T0[u] ~ dnorm(mu_t0, tau_t0)T(0,1e7)
}
# hyper parameter priors
mu_l ~ dunif(60, 200)
tau_l ~ dgamma(.01,.01)
mu_k ~ dunif(0.001,1)
tau_k ~ dgamma(.0001,.0001)
mu_t0 ~ dnorm (0, .1)
tau_t0 ~ dgamma(.0001,.0001)
}"
v<-c("mu_l","mu_k","mu_t0")

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
lamp.lmax.graph1<-mcmc_intervals_data(lamp.mcmc.data[,-1], regex_pars = "h_mu_l") %>%
  mutate(Site.Agg=Lsites[as.numeric(substr(parameter,8,9))]) %>%
  left_join(SiteID[,c(2,13)]) %>% filter(!duplicated(.)) %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
lamp.lmax.graph1
llmax.con<-ggplot(lamp.lmax.graph1) +
  #geom_ribbon(aes(x=Lat.cor,ymin=2.6003*Lat.cor, ymax=3.0162*Lat.cor),
  #            fill="lightgrey", alpha=.5)+
  geom_linerange(aes(ymin = ll, ymax = hh, x = Lat.cor))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = Lat.cor), size = 2)+ #inner line
  geom_point(aes(y = m, x = Lat.cor), size = 3, shape= 21,
             fill="darkgrey")+
  scale_x_continuous(name="Latitude", limits=c(29.33, 45.4))+
  scale_y_continuous(name="Potential\nMaximum Length (mm)", limits = c(60,160))+
  #geom_abline(intercept=0,slope=2.7917)+
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
  filter(!duplicated(.),
         Site.Agg != 'Wendell2') %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
almax.con<-ggplot(amb.lmax.graph1) +
  #geom_rect(data=abiorect, aes(xmin=xmin,xmax=xmax,
  #                               ymin=ymin,ymax=ymax), 
  #          fill=wes_palette("Darjeeling1")[c(1,4,3)])+
  #geom_image(aes(image = image), x=30.5, y=153, size=.15)+
  #geom_text(y=145, x=31, label=expression(italic("A. plicata")))+
  #geom_ribbon(aes(x=Lat.cor,ymin=2.95*Lat.cor, ymax=3.34*Lat.cor),
              #fill="lightgrey", alpha=.5)+
  geom_linerange(aes(ymin = ll, ymax = hh, x = Lat.cor))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = Lat.cor),
               size = 2)+ #inner line
  geom_point(aes(y = m, x = Lat.cor), size = 3, shape= 21, fill="white")+
  scale_x_continuous(name="Latitude")+
  scale_y_continuous(name="Potential\nMaximum Length (mm)")+
  #geom_abline(intercept=0,slope=3.135)+
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

colorvector<-c('#b3bd28','#ca6833','#00822d','#0069d6','#c64579')
lmax.graph<-bind_rows(amb.lmax.graph1 %>% mutate(Species="A. plicata"),
          lamp.lmax.graph1 %>% mutate(Species="Lampsilis spp.")) %>%
  mutate(SiteLat=fct_reorder(Site.Agg,Lat.cor))
bioregrect<-data.frame(xmin=c(0,  2.5,4.5,12.5,22.5,20.5,23.5),
                       xmax=c(2.5,4.5,12.5,22.5,23.5,21.5, 26.5), 
                       ymin=rep(55,7), ymax=rep(60,7))
lmax.dis<-ggplot(lmax.graph) +
  geom_rect(data=bioregrect, aes(xmin=xmin,xmax=xmax,
                                 ymin=ymin,ymax=ymax), 
            fill=colorvector[c(5,2,1,4,3,3,4)])+
  geom_linerange(aes(ymin = ll, ymax = hh, x = SiteLat,
                     group=Species), position=position_dodge(.75))+  #outer line
  geom_linerange(aes(ymin = l, ymax = h, x = SiteLat, group=Species),
               size = 2,position=position_dodge(.75))+ #inner line
  geom_point(aes(x = SiteLat, y = m, fill=Species, group=Species), 
             size = 3, shape= 21, 
             position=position_dodge(.75))+
  scale_y_continuous(name="Potential\nMaximum Length (mm)",
                     breaks=c(60,100,150,200),expand=c(0,0), limits=c(55,200))+
  scale_x_discrete(name="Sites")+
  scale_fill_manual(values=c("white","darkgrey"),
                    labels=c(expression(italic("A. plicata")),
                             expression(italic("L. cardium"))))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, size=8.5, hjust=1),
        legend.position = 'top')
lmax.dis
plot_grid(plot_grid(llmax.con, almax.con, ncol=1, 
                    labels=c("A","B")),
          lmax.dis, nrow=1, labels=c("","C"),
          rel_widths = c(.57,1.1))
ggsave('figures/fig2_211111.tiff', width=7, height=4)

# results numbers ------
amb.mcmc.sum %>% filter(variable=="mu_l", 
                        SiteID != 'Wendel2') %>% summarise(mean(x50), min(x50), max(x50))
amb.mcmc.sum %>% filter(variable=="mu_k", 
                        SiteID != 'Wendel2') %>% summarise(mean(x50), min(x50), max(x50))
amb.mcmc.sum %>% filter(variable=="mu_t0", 
                        SiteID != 'Wendel2') %>% summarise(mean(x50), min(x50), max(x50))

lamp.mcmc.sum %>% filter(variable=="h_mu_l") %>% summarise(mean(x50), min(x50), max(x50))
lamp.mcmc.sum %>% filter(variable=="h_mu_k") %>% summarise(mean(x50), min(x50), max(x50))
lamp.mcmc.sum %>% filter(variable=="mu_t0") %>% summarise(mean(x50), min(x50), max(x50))

# Table 2
lkt.table<-lamp.mcmc.sum %>%
  mutate(Sp="LAMP") %>%
  bind_rows(amb.mcmc.sum %>% mutate(Sp="APLI"))%>%
  mutate(Site.Agg=recode(SiteID, 
                         Wendel2="Wendell2", Wendel3="Wendell3")) %>%
  select(-SiteID) %>%
  left_join(SiteID) %>%
  select(Site.Agg, Lat.cor, variable, Sp, x50) %>%
  pivot_wider(names_from=variable, values_from=x50, values_fn=list(x50=mean)) %>%
  arrange(Sp,Lat.cor)
write.csv(lkt.table, "figures/Table2KLT.csv")


# Rensch's figure ###
#raw plots
lamp.mcmc.sum %>% filter(!is.na(sex_level), name=='par_k') %>%
  ggplot()+
  geom_linerange(aes(x=SiteID, ymin=x25, y=x50, ymax=x75, color=sex_level), 
                 position=position_dodge(width=.4), alpha=.3)+
  geom_point(aes(x=SiteID, y=x50, color=sex_level),
             position=position_dodge(width=.4))


ggsave('figures/sex_latitude.tiff', width=4.5, height=3)
# ratio at sites
#female / male ratio
lamp.mcmc.data<-read.csv("Lamp_Lmax_mcmcres.csv")
lamp.mcmc.sum<-read.csv("Lamp_Lmax_mcmc_sum.csv")
sex_rat<-lamp.mcmc.data %>%
  rowid_to_column() %>%
  select(rowid,starts_with('mu_l')) %>%
  pivot_longer(-rowid, names_to = 'parameter') %>%
  #mutate(nnn=nchar(parameter)) %>% arrange(desc(nnn))
  mutate(s=gsub(".*_\\s*|_.*","", parameter),
         variable=case_when(nchar(parameter)==10~
                              substr(parameter, 1,(nchar(parameter)-3)),
                            nchar(parameter)==7~
                              substr(parameter, 1,(nchar(parameter)-3)),
                            T~
                              substr(parameter, 1,(nchar(parameter)-2)))) %>%
  
  left_join(lamp.mcmc.sum %>% 
              filter(grepl('mu_l', variable)) %>%
              mutate(s=as.character(s))) %>%
  #filter(is.na(sex_level))
  group_by(rowid, SiteID) %>%
  select(rowid, SiteID,sex_level, value) %>%
  pivot_wider(names_from=sex_level, values_from=value)%>%
  mutate(FtM=`F`/M) 
plot_grid(
  lamp.mcmc.sum %>%  filter(grepl('mu_l', variable),
                            variable != 'h_mu_l',
                            sex_level !='U') %>%
    left_join(SiteID) %>%
    mutate(Sex=recode(sex_level, `F`='Female','M'='Male'))%>%
    ggplot()+
    geom_linerange(aes(x=Lat.cor, ymin=x25, y=x50, ymax=x75, color=Sex), 
                   position=position_dodge(width=.2))+
    geom_point(aes(x=Lat.cor, y=x50, color=Sex),
               position=position_dodge(width=.2))+
    facet_wrap(~Sex)+
    scale_color_grey('Sex')+
    scale_y_continuous('Maximum Length (mm)')+
    scale_x_continuous('Latitude')+
    theme_classic()+theme(legend.position = 'none'),
  sex_rat %>%
  filter(!is.na(FtM)) %>%
  group_by(SiteID) %>%
  select(SiteID, FtM) %>%
  summarize_all(quantile) %>%
  mutate(var=c('minimum','x25','x50','x75','maximum')) %>%
  pivot_wider(names_from=var, values_from=FtM) %>%
  left_join(SiteID) %>%
  ggplot()+
  geom_hline(yintercept=1, linetype=2, color='grey')+
    geom_linerange(aes(x=Lat.cor,y=x50, ymin=x25, ymax=x75))+
  geom_point(aes(x=Lat.cor, y=x50))+
  scale_x_continuous('Latitude')+
  scale_y_continuous('Female : Male Linf ratio')+
  theme_classic(),
  labels='AUTO', rel_widths = c(.65,.35))
ggsave('figures/lamp_sex_latitude.tiff', width=6.5, height=3)

lamp.mcmc.sum %>%
  filter(!is.na(sex_level)) %>%
  mutate(cell=paste0(round(x50,1), ' (',round(x2.5,1),', ', round(x97.5,1), ')')) %>%
  select(SiteID, sex_level, name, cell) %>%
  pivot_wider(names_from=name, values_from = cell) %>%
  left_join(SiteID) %>%
  select(SiteID, Latitude, Longitude, sex_level, par_k, par_lmax) %>%
  arrange(Latitude) %>% 
  write.csv('Lamp_sex_dif_table.csv')
AxL %>% group_by(id) %>% slice(1) %>% ungroup() %>% 
    filter(Species=='LCAR') %>% 
    group_by(Site) %>% count(Sex) %>%View()