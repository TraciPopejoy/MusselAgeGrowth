library(readxl);library(dplR);library(ggplot2)

SiteID<-read_excel('data/!GrowthRawData.xlsx', sheet="Location")

# take Linf for each population and match it with Latitude
latxVB <- VB %>% SiteID

library(brms)
set.seed(6363) #set random number seed for replicable analysis

### testing priors 
hist(rstudent_t(1000, 3,0,10)) #default prior on sd(TankIntercept)
hist(rgamma(1000,0.01,0.01), xlim=c(0,70)) # prior on shape of Treat:Day relationship
hist(rnorm(1000, 0,5)) # prior on beta (Treat:Day interaction)

#### Ammonium analysis ####
model<-brm(Linf~Latitude-1, data=latxVB,
              family=Gamma(link='log'), 
              prior=c(prior(normal(0,5), class="b"),
                      prior(gamma(0.01,0.01), class="shape")),
              chains=4, iter=2000)

print(model, prior=T)
plot(marginal_effects(model)) #check it is reproducing data well
#pull out posterior samples for each parameter
post_model<-posterior_samples(model) 

# quantifying the % probability response increased after impact
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1 ) 
nhBACIgraph %>% group_by(ratio) %>% dplyr::summarize(meanBACI=mean(value))
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1.5 ) 
marginal_effects(nhBmodel)$`TreF:DayF` %>% 
  select(TreF, DayF, estimate__, lower__,upper__) %>%
  filter(DayF==4)