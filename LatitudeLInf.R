head(Amcmc.data); head(Lmcmc.data)

### testing priors 
hist(rnorm(1000,0,10)) # prior on beta

beta_model_string<- "model {
# likelihood
for (i in 1:nobs){ # each row
Latitude[i] ~ dnorm(lat[i],tau) 
lat[i]<- beta*mu_l[i]
}
# priors
tau ~ dgamma(0.1, 0.1) #variance on L
beta ~ dnorm(0,20)
}"

# Lampsilis x Latitude model -----------
#loop lengths
nobsL<-nrow(Lmcmc.data)

#run the model
L.beta.model<-jags.model(textConnection(beta_model_string), 
                    data=list(mu_l=lamp.mcmc.data$mu_l,
                              Latitude = lamp.mcmc.data$Latitude,
                              nobs=nobsL),
                    n.chains=3, n.adapt=2000)
update(L.beta.model, 20000) # burn in for 2000 samples
Lbeta.mcmc<-coda.samples(L.beta.model,
                    variable.names=c("beta"), 
                    n.iter=100000, thin=3)
pdf('LbetaMCMCdiag.pdf')
plot(Lbeta.mcmc)
gelman.plot(Lbeta.mcmc)
dev.off()