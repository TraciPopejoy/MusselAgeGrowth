
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
}"

# Amblema x Latitude model ----------

beta.model<-jags.model(textConnection(beta_model_string), 
                    data=list(mu_l = lmax.long$x50,
                              Latitude = lmax.long$Lat.cor,
                              SpF=lmax.long$SpF,
                              nobs=40),
                    inits=list(beta=c(.5,.5)),
                    n.chains=3, n.adapt=1000)
update(beta.model, 2000) # burn in for 2000 samples
beta.mcmc<-coda.samples(beta.model,
                    variable.names=c("beta","invbeta"), 
                    n.iter=10000, thin=3)
pdf('betaMCMCdiag.pdf')
plot(beta.mcmc)
gelman.plot(beta.mcmc)
dev.off()

gelman.diag(beta.mcmc)[[1]]

summary(beta.mcmc)
library(bayesplot)
color_scheme_set("red")
pdf('betares.pdf')
mcmc_areas(beta.mcmc, regex_pars = "invbeta")+
  scale_y_discrete(labels=levels(lmax.long$SpF))
dev.off()

beta.mcmc.data<-as.matrix(beta.mcmc)
1-ecdf(beta.mcmc.data[,2])( 1 ) #100%

lmax.long %>% 
  group_by(SpF) %>%
  summarize(meanLmax=mean(x50),
            maxLmax=max(x50))

lm(x50~Lat.cor-1, data=lmax.long)
pdf('testing.pdf')
ggplot(lmax.long,aes(x=Lat.cor, y=x50))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~SpF)
dev.off()
