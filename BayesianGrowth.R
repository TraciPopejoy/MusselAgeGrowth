# Bayesian von Bertanlanffy's
set.seed(6363) #set random number seed for replicable analysis

### testing priors 
hist(rnorm(1000, 1,20)) # prior on beta
hist(runif(1000, 10,200)) # prior on Lmax
hist(rnorm(1000, .5,.5)) # prior on K
hist(rnorm(1000, 0,.1)) # prior on T0
hist(rgamma(1000,1,10)) #variance on Latitude
hist(rgamma(1000,.1,.1)) #variance on L
hist(rgamma(1000,1,1)) #variance on mu_l
hist(rgamma(1000,.1,.1),breaks=800,xlim=c(0,1))
hist(rgamma(1000,.01,.01),breaks=800,xlim=c(0,1)) #variance on tauk


# Von Bertanlanffy's model
# L = Lmax*(1-exp(-K*(Age-t0)))

library(rjags);library(MCMCpack)
##### Model Code #####
model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
L[i] ~ dnorm(mu[i],tau) 
mu[i] <- Lmax[idF[i],SiteF[i]] * (1-exp(-K[idF[i],SiteF[i]] * (Age[i] - T0[idF[i],SiteF[i]])))
Latitude[i] ~dnorm(Lat_mu[i], tau_lat)
Lat_mu[i]<-beta*Lmax_bar[SiteF[i]]
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L
tau_lat~dgamma(1,10) #variance on Latitude
beta ~ dnorm(0,20)

for (s in 1:nsite) { # loop over sites
for(u in 1:nid){
# so Lmax, K and T0 will be unique for each site
Lmax[u,s] ~ dnorm(mu_l[s], tau_l[s])
K[u,s] ~ dnorm(muk[s], tauk[s])
T0[u,s] ~ dnorm(mu_t0[s], tau_t0[s])
}

#get average Lmax, K, T0 for each site
Lmax_bar[s]<- mean(Lmax[1:nid,s])
K_bar[s]<- mean(K[1:nid,s])
T0_bar[s]<- mean(T0[1:nid,s])

# hyper parameter priors
mu_l[s] ~ dunif(40, 200)
tau_l[s] ~ dgamma(.1,.1)

muk[s] ~ dnorm (.5, .5)
tauk[s] ~ dgamma(.1,.1)

mu_t0[s] ~ dnorm (0, .1)
tau_t0[s] ~ dgamma(.01,.01)
}
}"

# Lampsilis model -----------
#loop lengths
nidL<-length(unique(dataL$idF))
nobsL<-nrow(dataL)
nsiteL<-length(unique(dataL$SiteF))

#run the model
L.model<-jags.model(textConnection(model_string), 
                    data=list(L = dataL$L,
                              Age = dataL$Age,
                              Latitude = dataL$Latitude,
                              SiteF = dataL$SiteF,
                              idF = dataL$idF,
                              nid=nidL,
                              nobs=nobsL,
                              nsite=nsiteL),
                     inits=list(Lmax=matrix(rep(100,1536),nrow=96),
                                K=matrix(rep(.5,1536),nrow=96),
                                T0=matrix(rep(.001,1536),nrow=96)),
                    n.chains=3, n.adapt=2000)
update(L.model, 4000) # burn in for 2000 samples
Lmcmc1<-coda.samples(L.model,
                     variable.names=c("Lmax_bar", "beta", "K_bar","T0_bar"), 
                     n.iter=10000, thin=6)
pdf('LampMcmcDiag_sim.pdf')
plot(Lmcmc1)
gelman.plot(Lmcmc1)
dev.off()

str(summary(Lmcmc1))
postLamp<-data.frame(variable_long=rownames(summary(Lmcmc)$statistics),
                     variable=substr(rownames(summary(Lmcmc)$statistics),1,4),
                     mean=summary(Lmcmc)$statistics[,1],
                     SD=summary(Lmcmc)$statistics[,2],
                     SiteID=c(rep(levels(unique(dataL$SiteF)),3),"beta"),
                     x2.5=summary(Lmcmc)$quantiles[,1],
                     x97.5=summary(Lmcmc)$quantiles[,5],
                     x50=summary(Lmcmc)$quantiles[,3]) %>%
  left_join(SiteID)
Lmcmc.data<-as.matrix(Lmcmc); colnames(Lmcmc.data)

ggplot(data=postLamp[postLamp$variable!="beta",])+
  geom_crossbar(aes(x=Latitude,y=x50,ymin=x2.5, ymax=x97.5))+
  geom_point(aes(x=Latitude, y=mean))+
  facet_wrap(~variable, scales="free_y")

# Amblema model -----------
#loop lengths
nidA<-length(unique(dataA$idF))
nobsA<-nrow(dataA)
nsiteA<-length(unique(dataA$SiteF))

#run the model
A.model<-jags.model(textConnection(model_string), 
                    data=list(L = dataA$L,
                              Age = dataA$Age,
                              Latitude = dataA$Latitude,
                              SiteF = dataA$SiteF,
                              idF = dataA$idF,
                              nid=nidA,
                              nobs=nobsA,
                              nsite=nsiteA),
                    inits=list(Lmax=matrix(rep(100,2736),nrow=24),
                               K=matrix(rep(.5,2736),nrow=24),
                               T0=matrix(rep(.001,2736),nrow=24)),
                    n.chains=3, n.adapt=100000)
update(A.model, 10000) # burn in for 2000 samples
Amcmc<-coda.samples(A.model,
                     variable.names=c("Lmax_bar", "beta", "K_bar","T0_bar"), 
                     n.iter=5000, thin=40)
pdf('AMcmcDiag.pdf')
plot(Amcmc)
gelman.plot(Amcmc)
dev.off()

summary(Amcmc)
Amcmc.data<-as.matrix(Amcmc); colnames(Amcmc.data)


# PLOTS ######
#pull results from summary(mcmc) to get into long dataframe
# want to make a point with mean of each sitexspecies
# want to plot 95% credible intervals with those points
ggplot(data=long_mcmc)+
  geom_point()+
  geom_hline()