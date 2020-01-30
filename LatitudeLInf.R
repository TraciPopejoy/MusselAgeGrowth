vl<-c(paste("mu_l",1:16, sep="_"))
names(vl)<-Lsites
va<-c(paste("mu_l", 1:24, sep="_"))
names(va)<-Asites
amb.mcmc.chains<-read.csv("Amb_Lmax_mcmcres.csv") %>%
  dplyr::sample_frac(.10) %>%
  dplyr::select(starts_with("mu_l")) %>%
  dplyr::rename(all_of(va)) %>%
  gather(Site.Agg, x50) %>%
  mutate(Sp="APLI", SpF=factor(Sp)) %>%
  left_join(SiteID) %>%
  dplyr::select(SpF,Site.Agg, Lat.cor, x50)
lam.mcmc.chains<-read.csv("Lamp_Lmax_mcmcres.csv") %>%
  dplyr::sample_frac(.10) %>%
  dplyr::select(starts_with("mu_l")) %>%
  dplyr::rename(all_of(vl)) %>%
  gather(Site.Agg, x50) %>%
  left_join(SiteID) %>%
  mutate(Sp="LAMP", SpF=factor(Sp)) %>%
  dplyr::select(SpF, Site.Agg, Lat.cor, x50)
mcmc.chains<-rbind(amb.mcmc.chains, lam.mcmc.chains) %>%
  dplyr::sample_frac(.33)
mcmc.chains %>% group_by(SpF) %>%
  tally()
head(mcmc.chains)

amb.mcmc.sum<-read.csv("Amb_Lmax_mcmc_sum.csv")
lamp.mcmc.sum<-read.csv("Lamp_Lmax_mcmc_sum.csv")
lmax.long<-lamp.mcmc.sum %>%
  mutate(Sp="LAMP") %>%
  bind_rows(amb.mcmc.sum %>% mutate(Sp="APLI"))%>%
  rename(Site.Agg=SiteID) %>%
  left_join(SiteID) %>%
  filter(variable=="mu_l")%>%
  mutate(SpF=factor(Sp)) %>%
  dplyr::select(SpF,Site.Agg, x50, Lat.cor) %>%
  filter(!duplicated(.))
#fire emblem harry potter vs. got (politics)
#rimworld - like stardew valley

### testing priors 
hist(rnorm(1000,0,5)) # prior on beta
hist(rgamma(1000,0.01, 0.01), breaks=800, xlim=c(0,2)) #precision on Lat


beta_model_string<- "model {
# likelihood
for (i in 1:nobs){ 
Latitude[i] ~ dnorm(lat[i],tau) 
lat[i] <- beta[SpF[i]]*mu_l[i]
}
# priors
tau ~ dgamma(0.01, 0.01) #variance on Lat
for(s in 1:2){
beta[s] ~ dnorm(mu_b, tau_b)
}
#hyper priors
mu_b ~ dnorm(0,5)
tau_b ~ dgamma(.01,.01)
#transformation
invbeta<-1/beta
difbeta<-invbeta[1]-invbeta[2]
}"

# Amblema x Latitude model ----------

beta.model<-jags.model(textConnection(beta_model_string), 
                    data=list(mu_l = lmax.long$x50,
                              Latitude = lmax.long$Lat.cor,
                              SpF=lmax.long$SpF,
                              nobs=40),
                    inits=list(beta=c(.5,.5)),
                    n.chains=3, n.adapt=1000)
#testing using mcmc chains
#beta.model<-jags.model(textConnection(beta_model_string), 
#                       data=list(mu_l = mcmc.chains$x50,
#                                 Latitude = mcmc.chains$Lat.cor,
#                                 SpF=mcmc.chains$SpF,
#                                 nobs=54285),
#                       inits=list(beta=c(.5,.5)),
#                       n.chains=3, n.adapt=1000)
update(beta.model, 2000) # burn in for 2000 samples
beta.mcmc<-coda.samples(beta.model,
                    variable.names=c("beta","invbeta","difbeta"), 
                    n.iter=10000, thin=3)
pdf('betaMCMCdiag.pdf')
plot(beta.mcmc)
gelman.plot(beta.mcmc)
dev.off()

gelman.diag(beta.mcmc[,-3])[[1]]

summary(beta.mcmc)
library(bayesplot)
color_scheme_set("blue")
mcmc_areas(beta.mcmc, regex_pars = c("invbeta"))+
  scale_y_discrete(labels=c("A. plicata","Lampsilis spp."))+
  theme_classic()

beta.plot<-mcmc_intervals(beta.mcmc, regex_pars = c("invbeta"))+
  scale_y_discrete(labels=c("A. plicata","Lampsilis spp.","Difference"))+
  theme_classic()+
  xlab("Lmax = Beta*Latitude")

difbeta.plot<-mcmc_areas(beta.mcmc, pars="difbeta")+
  theme_classic()+
  xlab("Difference in Beta")+
  scale_y_discrete(labels="Lampsilis spp.")+
  theme(axis.text.y = element_text(color="white"))

library(cowplot)
plot_grid(beta.plot,difbeta.plot, ncol=1, labels="AUTO")
ggsave("Figure3.tiff", width=3.3, height=4)


beta.mcmc.data<-as.matrix(beta.mcmc)
beta.cum<-ecdf(beta.mcmc.data[,3])
summary(beta.cum)
1-beta.cum(0)


