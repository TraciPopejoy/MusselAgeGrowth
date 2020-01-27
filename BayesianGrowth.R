# Bayesian von Bertanlanffy's
set.seed(6363) #set random number seed for replicable analysis

### testing priors 
hist(rnorm(1000, 1,20)) # prior on beta
hist(runif(1000, 10,200)) # prior on Lmax
hist(rnorm(1000, .5,.1)) # prior on K
hist(rnorm(1000, 0,.1)) # prior on T0
hist(rgamma(1000,1,10)) #variance on Latitude
hist(rgamma(1000,.1,.1)) #variance on L
hist(rgamma(1000,1,1)) #variance on mu_l
hist(rgamma(1000,.1,.1),breaks=800,xlim=c(0,1))
hist(rgamma(1000,.01,.01),breaks=800,xlim=c(0,1)) #variance on tauk


# Von Bertanlanffy's model
# L = Lmax*(1-exp(-K*(Age-t0)))
LampD<-LAMPaxl %>% left_join(SiteID, by=c('Site'='SiteID')) %>%
  dplyr::select(Site, id, Age, L, Latitude) %>%
  ungroup()%>%
  mutate(SiteF=as.factor(Site),
         idF=as.factor(id)) %>%
  dplyr::select(-Site, -id)

library(rjags);library(MCMCpack)
##### Model Code #####
model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
L[i] ~ dnorm(mu[i],tau) 
mu[i] <- Lmax[Species[i],SiteF[i],idF[i]] * (1-exp(-K[Species[i],SiteF[i],idF[i]] * (Age[i] - T0[Species[i],SiteF[i],idF[i]])))
Latitude[i] ~dnorm(Lat_mu[i], tau_lat)
Lat_mu[i]<-beta[Species[i]]*Lmax_bar[Species[i],SiteF[i]]
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L
tau_lat~dgamma(1, 10) #variance on Latitude

for(y in 1:2){
for (s in 1:nsite) { # loop over sites
for(u in 1:nid){
# so Lmax, K and T0 will be unique for each site
Lmax[y,s,u] ~ dnorm(mu_l[y,s], tau_l[y,s])
K[y,s,u] ~ dnorm(muk[y,s], tauk[y,s])
T0[y,s,u] ~ dnorm(mu_t0[y,s], tau_t0[y,s])
}

#get average Lmax, K, T0 for each site
Lmax_bar[y,s]<- mean(Lmax[y,s, 1:210])
K_bar[y,s]<- mean(K[y,s, 1:210])
T0_bar[y,s]<- mean(T0[y,s, 1:210])

# hyper parameter priors
mu_l[y,s] ~ dunif(10, 200)
tau_l[y,s] ~ dgamma(1,1)

muk[y,s] ~ dnorm (0.5, 1)
tauk[y,s] ~ dgamma(.1,.1)

mu_t0[y,s] ~ dnorm (0, 0.1)
tau_t0[y,s] ~ dgamma(.01,.01)
}
beta[y] ~ dnorm(0,20)
}
}"
#loop lengths
nid<-length(unique(data$idF))
nobs<-nrow(data)
nsite<-length(unique(data$SiteF))

#can't figure out the dimensions for these so removing
inits<-list(Lmax=matrix(rep(100,2839),nrow=28),
            K=rep(.5,2839),
            T0=rep(.001,2839),nrow=28))

#run the model
NHmodel<-jags.model(textConnection(model_string), 
                    data=list(Species=data$SpF,
                              L=data$L,
                              Age=data$Age,
                              Latitude=data$Latitude,
                              SiteF=data$SiteF,
                              idF=data$idF,
                              nid=nid,
                              nobs=nobs,
                              nsite=nsite),
                    #  inits=inits,
                    n.chains=3, n.adapt=1000)
update(NHmodel, 1000) # burn in for 2000 samples
NHmcmc<-coda.samples(NHmodel,
                     variable.names=c("Lmax_bar", "beta",
                                      "K_bar","T0_bar"), 
                     n.iter=5000, thin=40)
pdf('mcmcdiag.pdf')
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

# PLOTS ######
#pull results from summary(mcmc) to get into long dataframe
# want to make a point with mean of each sitexspecies
# want to plot 95% credible intervals with those points
ggplot(data=long_mcmc)+
  geom_point()+
  geom_hline()