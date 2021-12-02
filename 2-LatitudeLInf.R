vl<-c(paste("h_mu_l",1:12, sep="_"))
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
  dplyr::select(starts_with("h_mu_l")) %>%
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
  mutate(Sp="LAMP") %>% select(-X) %>%
  filter(variable=="h_mu_l")%>%
  bind_rows(amb.mcmc.sum %>% mutate(Sp="APLI") %>% 
              select(-X) %>%
              filter(variable=="mu_l"))%>%
  mutate(Site.Agg=recode(SiteID, 
                       Wendel2="Wendell2", Wendel3="Wendell3")) %>%
  filter(Site.Agg != 'Wendell2') %>%
  select(-SiteID) %>%
  left_join(SiteID) %>%
  group_by(Sp) %>%
  mutate(mu_l_z=scale(x50, center=F)[,1],
         Latscale=scale(Lat.cor, center=F)[,1],
         precision=1/SD^2) %>%
  ungroup() %>%
  mutate(SpF=factor(Sp), 
         deltax50=case_when(SpF=="APLI"~x50/173*100,
                            SpF=="LAMP"~x50/131*100))%>%
  dplyr::select(SpF,Site.Agg, x50, Lat.cor, mu_l_z, 
                Latscale, deltax50, precision)%>%
  filter(!duplicated(.))

lhist<-lmax.long %>% 
  dplyr::select(SpF, Lat.cor, x50, mu_l_z, deltax50) %>%
  gather(variable, value, -SpF, -Lat.cor) 
ggplot(lhist)+geom_histogram(aes(x=value))+
  facet_wrap(~SpF+variable, scales="free")

### testing priors 
hist(rnorm(1000,0,10)) # prior on beta 
hist(rnorm(1000,0,50)) # prior on alpha 

old_model_string<- "model {
# likelihood
for (i in 1:nobs){ 
Latitude[i] ~ dnorm(lat[i],tau) 
lat[i] <- beta[SpF[i]]*mu_l[i]
}
# priors
tau <- dgamma(0.01, 0.01) #variance on Lat
for(s in 1:2){
beta[s] ~ dnorm(mu_b, tau_b)
}

#hyper priors
mu_b ~ dnorm(0,10)
tau_b ~ dgamma(.01,.01)

#transformation
invbeta<-1/(beta)
difbeta<-invbeta[1]-invbeta[2]
}"

test_model_string<- "model {
# likelihood
for (i in 1:nobs){ 
mu_l[i] ~ dnorm(mu_e[i],tau[i]) 
mu_e[i] <- beta[SpF[i]]*Latitude[i] + alpha
}
# priors
alpha ~ dnorm(0,50)T(-100,100)
for(s in 1:2){
beta[s] ~ dnorm(mu_b, tau_b)
}

#hyper priors
mu_b ~ dnorm(0,10)
tau_b ~ dgamma(.01,.01)

#transformation
difbeta<-beta[1]-beta[2]
}"

# Amblema x Latitude model ----------
library(rjags)
beta.model<-jags.model(textConnection(test_model_string), 
                    data=list(mu_l = lmax.long$deltax50,
                              Latitude = lmax.long$Lat.cor,
                              SpF=lmax.long$SpF,
                              tau=lmax.long$precision,
                              nobs=35),
                    n.chains=3, n.adapt=10000)
update(beta.model, 30000) # burn in for 2000 samples
beta.mcmc<-coda.samples(beta.model,
                    variable.names=c("beta","difbeta","alpha"), 
                    n.iter=500000, thin=3)
pdf('betaMCMCdiag.pdf')
plot(beta.mcmc)
gelman.plot(beta.mcmc)
dev.off()

gelman.diag(beta.mcmc, multivariate = F)

summary(beta.mcmc)
library(bayesplot)
color_scheme_set("gray")

mcmc_intervals(beta.mcmc, regex_pars=c('alpha'))+ 
  theme_classic()
mcmc_intervals(beta.mcmc, pars=c('beta[1]', 'beta[2]'))+
  scale_y_discrete(labels=c("A. plicata","Lampsilis spp."))+
  theme_classic()

beta.plot<-mcmc_intervals(beta.mcmc, pars = c("beta[1]", "beta[2]"))+
  scale_y_discrete(labels=c(expression(italic("A. plicata")),
                            expression(italic("L. cardium"))))+
  theme_classic()+
  ggtitle("% Max. Length = slope * Latitude + intercept")+
  xlab(expression("% Max. Length"%.%"Latitude"^-1))+
  theme(title=element_text(size=8))

difbeta.plot<-mcmc_areas(beta.mcmc, pars="difbeta")+
  theme_classic()+
  ggtitle(expression(italic("A. plicata")*" slope - "*italic("L. cardium")*" slope"))+
  scale_y_discrete(labels="Lampsilis spp.")+
  theme(axis.text.y = element_text(color="white"),
        title = element_text(size=9))

plot_grid(beta.plot,difbeta.plot, ncol=1, labels="AUTO",
          rel_heights= c(1,.75), align=c('v'))
ggsave("figures/Figure4beta.tiff", width=3.5, height=3.5)

summary(lmax.long[lmax.long$SpF=="APLI",]$x50);summary(lmax.long[lmax.long$SpF=="APLI",]$mu_l_z)
summary(lmax.long[lmax.long$SpF=="LAMP",]$x50);summary(lmax.long[lmax.long$SpF=="LAMP",]$mu_l_z)

beta.mcmc.data<-as.matrix(beta.mcmc)
beta.cum<-ecdf(beta.mcmc.data[,"difbeta"])
summary(beta.cum)
beta.cum(0) #area under the curve <0


#Appendix - bootstrapped linear model
#from tutorial: https://www.statmethods.net/advstats/bootstrapping.html
library(boot);library(broom)
# function to obtain R-Squared from the data 
coefff <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(c(coef(fit),summary(fit)$r.square,
                     glance(fit)$p.value))
} 
# bootstrapping with 1000 replications 
results <- boot(data=lmax.long, statistic=coefff, 
                R=1000, formula=deltax50~Lat.cor:SpF)

# view results
results 
plot(results, index=1) #intercept
plot(results, index=2) #amb slope
plot(results, index=3) #lam slope
plot(results, index=4) #rsquared
plot(results, index=5) #pvalue

# get 95% confidence interval 
boot.ci(results, type="bca", index=1)
boot.ci(results, type="bca", index=2)
boot.ci(results, type="bca", index=3)
boot.ci(results, type="bca", index=4)
boot.ci(results, type="bca", index=5)
