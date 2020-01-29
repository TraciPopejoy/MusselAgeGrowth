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
hist(rgamma(1500,.001,.001),breaks=1300,xlim=c(0,1)) #variance on tauk
hist(rgamma(1000,.01,.01),breaks=800,xlim=c(0,1)) #variance on tau_t0


# Von Bertanlanffy's model
# L = Lmax*(1-exp(-K*(Age-t0)))

library(rjags);library(MCMCpack)
##### Model Code #####
model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
L[i] ~ dnorm(mu[i],tau) 
mu[i] <- Lmax[idF[i],SiteF[i]] * (1-exp(-K[idF[i],SiteF[i]] * (Age[i] - T0[idF[i],SiteF[i]])))
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L

for (s in 1:nsite) { # loop over sites
for(u in 1:nid){
# so Lmax, K and T0 will be unique for each site and individual
Lmax[u,s]~ dnorm(mu_l[s], tau_l[s])
K[u,s] ~ dnorm(muk[s], tauk[s])
T0[u,s] ~ dnorm(mu_t0[s], tau_t0[s])
}

# hyper parameter priors
mu_l[s] ~ dunif(50, 180)
tau_l[s] ~ dgamma(.01,.01)

muk[s] ~ dnorm (.5, .5)
tauk[s] ~ dgamma(.0001,.0001)

mu_t0[s] ~ dnorm (0, .1)
tau_t0[s] ~ dgamma(.0001,.0001)
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
                              SiteF = dataL$SiteF,
                              idF = dataL$idF,
                              nid=nidL,
                              nobs=nobsL,
                              nsite=nsiteL),
                     inits=list(Lmax=matrix(rep(100,1536),nrow=96),
                                K=matrix(rep(.5,1536),nrow=96),
                                T0=matrix(rep(.001,1536),nrow=96)),
                    n.chains=3, n.adapt=2000)
update(L.model, 20000) # burn in for 200000 samples
Lmcmc<-coda.samples(L.model,
                     variable.names=c("mu_l","muk","mu_t0"), 
                     n.iter=100000, thin=3)
pdf('LampMcmcDiag.pdf')
plot(Lmcmc)
gelman.plot(Lmcmc)
dev.off()

str(summary(Lmcmc))
postLamp<-data.frame(variable_long=rownames(summary(Lmcmc)$statistics),
                     variable=substr(rownames(summary(Lmcmc)$statistics),1,4),
                     mean=summary(Lmcmc)$statistics[,1],
                     SD=summary(Lmcmc)$statistics[,2],
                     SiteID=c(rep(levels(unique(dataL$SiteF)),3)),
                     x2.5=summary(Lmcmc)$quantiles[,1],
                     x97.5=summary(Lmcmc)$quantiles[,5],
                     x50=summary(Lmcmc)$quantiles[,3]) %>%
  left_join(SiteID)
Lmcmc.data<-as.matrix(Lmcmc)

ggplot(data=postLamp)+
  geom_crossbar(aes(x=Latitude,y=x50,ymin=x2.5, ymax=x97.5))+
  geom_point(aes(x=Latitude, y=mean),size=2)+
  facet_wrap(~variable, scales="free_y")+
  geom_smooth(aes(x=Latitude, y=mean),method="lm", se=F)+
  theme_bw()

# Amblema model -----------
#loop lengths
nidA<-length(unique(dataA$idF))
nobsA<-nrow(dataA)
nsiteA<-length(unique(dataA$SiteF))

#run the model
A.model<-jags.model(textConnection(model_string), 
                    data=list(L = dataA$L,
                              Age = dataA$Age,
                              SiteF = dataA$SiteF,
                              idF = dataA$idF,
                              nid=nidA,
                              nobs=nobsA,
                              nsite=nsiteA),
                    inits=list(Lmax=matrix(rep(100,2736),nrow=114),
                               K=matrix(rep(.5,2736),nrow=114),
                               T0=matrix(rep(.001,2736),nrow=114)),
                    n.chains=3, n.adapt=2000)
update(A.model, 20000) # burn in for 2000 samples
Amcmc<-coda.samples(A.model,
                     variable.names=c("mu_l","muk","mu_t0"), 
                     n.iter=100000, thin=3)
pdf('AMcmcDiag.pdf')
plot(Amcmc)
gelman.plot(Amcmc)
dev.off()

str(summary(Amcmc))
postAmb<-data.frame(variable_long=rownames(summary(Amcmc)$statistics),
                     variable=substr(rownames(summary(Amcmc)$statistics),1,4),
                     mean=summary(Amcmc)$statistics[,1],
                     SD=summary(Amcmc)$statistics[,2],
                     SiteID=c(rep(levels(unique(dataL$SiteF)),3)),
                     x2.5=summary(Amcmc)$quantiles[,1],
                     x97.5=summary(Amcmc)$quantiles[,5],
                     x50=summary(Amcmc)$quantiles[,3]) %>%
  left_join(SiteID)
Amcmc.data<-as.matrix(Amcmc)

ggplot(data=postAmb)+
  geom_crossbar(aes(x=Latitude,y=x50,ymin=x2.5, ymax=x97.5))+
  geom_point(aes(x=Latitude, y=mean),size=2)+
  facet_wrap(~variable, scales="free_y")+
  geom_smooth(aes(x=Latitude, y=mean),method="lm", se=F)+
  theme_bw()