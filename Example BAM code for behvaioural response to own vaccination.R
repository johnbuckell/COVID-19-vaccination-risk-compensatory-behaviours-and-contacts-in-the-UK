
#####################################################################################################################
# Example of code written by John Buckell and Koen Pouwels (Health Economics Research Centre, University of Oxford) #
# for analyses of behavioural outcomes and vaccination                                                              #
# Accompanying paper/preprint: example R scripts used for paper COVID-19 vaccination, risk-compensatory behaviours, #
# and contacts in the UK                                                                                            #
#####################################################################################################################

#####################
###### load packages
##################### 

library(mgcv)
library(dplyr)

#####################
###### estimate model
##################### 

prepost_bam <- bam(
  contact_physical_18y_scalar~
                    s(time_from_vaccination,k=20)
                   +s(age,k=30)  
                   +ever_lthc
                   +sex+non_white+country_imd_rank+ur2+ur3+ur4+hhsize2+hhsize3+hhsize4+hhsize5+multigen_collapsed
                   +gor9d
                   +vaccine_type
                   +lockdown_stage
                   +s(time_cal,by=gor9d,k=50),
                   data=prepost_data,
                   family=ocat(R=5),
                   method="fREML",
                   discrete=TRUE,
                   nthreads=8
)
summary(prepost_bam)
AIC(prepost_bam)
plot(prepost_bam,pages=1,shade=2,scale=0)


#################################################################
###### predicting behavior over time with posterior probabilities
#################################################################

### first vaccination
newdat<-with(prepost_data,expand.grid(
  ever_lthc=0,
  vaccine_type="oxford_az",
  sex="Female",
  white=1,
  lockdown_stage="easing_2",
  ever_hsc=0,
  ur2=0,
  ur3=0,
  ur4=0,
  hhsize2=1,
  hhsize3=0,
  hhsize4=0,
  hhsize5=0,
  multigen_collapsed=0,
  age=mean(age,na.rm=TRUE),
  gor9d="E12000007",
  country_imd_rank=mean(country_imd_rank,na.rm=TRUE),
  time_cal=mean(time_cal,na.rm=TRUE),
  time_from_vaccination=seq(min(time_from_vaccination),max(time_from_vaccination))  
  
))

fit <- predict(prepost_bam,newdata=newdat,type="response",se=TRUE,discrete=FALSE)

newdat <-cbind(newdat,fit)
newdat <- newdat %>%
  dplyr::mutate(some_contacts=1-fit.1) %>%
  dplyr::mutate(lcb_some_contacts=some_contacts-(1.96*se.fit.1)) %>%
  dplyr::mutate(ucb_some_contacts=some_contacts+(1.96*se.fit.1)) 
names(newdat)
summary(newdat$some_contacts)
temp_prob<-summary(newdat$some_contacts)
assign(paste0("prob_",dep),temp_prob)

ggplot(newdat,aes(x=time_from_vaccination,y=some_contacts))+
  geom_ribbon(aes(ymin=lcb_some_contacts,ymax=ucb_some_contacts),alpha=0.2) +
  coord_cartesian(ylim=c(0,1)) +
  geom_vline(xintercept=0) +
  geom_vline(xintercept=84,linetype="dotted") +
  labs(x="Time (days) from first vaccination", y='Probabiity of some contacts') +
  theme_bw() +
  geom_line(color="black")


############################################################
#### slope testing across range(s) of time since vaccination
############################################################
# code adapted from: https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# function generating random values from a multivariate normal
rmvn<-function(n,mu,sig){
  L<-mroot(sig) 
  m<-ncol(L) 
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

# extract elements from fitted BAM
Vb<-vcov(prepost_bam)
newd<-with(prepost_data,data.frame(
  ever_lthc=0,
  vaccine_type="oxford_az",
  sex="Female",
  white=1,
  lockdown_stage="easing_2",
  ever_hsc=0,
  ur2=0,
  ur3=0,
  ur4=0,
  hhsize2=1,
  hhsize3=0,
  hhsize4=0,
  hhsize5=0,
  multigen_collapsed=0,
  age=mean(age,na.rm=TRUE),
  gor9d="E12000007",
  country_imd_rank=mean(country_imd_rank,na.rm=TRUE),
  time_cal=mean(time_cal,na.rm=TRUE),
  time_from_vaccination=seq(min(time_from_vaccination),max(time_from_vaccination))
))
pred<-predict(prepost_bam,newd,se.fit=TRUE)
se.fit<-pred$se.fit

# generate simulations of the max abs standardised deviation of the fitted from the true model
set.seed(1354)
N<-10000
BUdiff<-rmvn(N,mu=rep(0,nrow(Vb)),sig=Vb)
Cg<- predict(prepost_bam, newd,type="lpmatrix")

# use to compute simulated confidence interval and compare to the model's predicted
simDev<-Cg %*% t(BUdiff)
absDev<- abs(sweep(simDev, 1, se.fit, FUN="/"))
masd<-apply(absDev, 2L, max)
crit<- quantile(masd, prob=0.95, type=8)
pred<- transform(cbind(data.frame(pred),newd),
                 uprP=fit+(1.96*se.fit),
                 lwrP=fit-(1.96*se.fit),
                 uprS=fit+(crit*se.fit),
                 lwrS=fit-(crit*se.fit)
)

# simulations based on draws from the posterior distribution of the fitted model
set.seed(1354)
N<-10000
Cg<- predict(prepost_bam, newd,type="lpmatrix")
pp<- Cg %*%coef(prepost_bam)
sims<- rmvn(N, mu = coef(prepost_bam), sig=Vb)
fits<- Cg %*% t(sims)

nrnd<-10000
rnd<- sample(N, nrnd)
stackFits<- stack(as.data.frame(fits[, rnd]))
stackFits<- transform(stackFits, time_from_vaccination = rep(newd$time_from_vaccination, length=(length(newd$time_from_vaccination))))
stackFits<- stackFits %>%
  mutate(prob_response1=(exp(prepost_bam$family$getTheta(TRUE)[3]-values)/(1+exp(prepost_bam$family$getTheta(TRUE)[3]-values)))) %>%
  mutate(prob_some_contacts=(1-prob_response1)) %>%
  dplyr::select(-values,-prob_response1)
summary(stackFits$prob_some_contacts)

# plot simulations (include simulated CIs if computed above)
ggplot(pred, aes(x=time_from_vaccination)) +
  geom_path(data=stackFits, mapping=aes(y = prob_some_contacts, x = time_from_vaccination, group = ind), alpha=0.4, colour="grey20")

# compute slopes in periods for each draw
class(stackFits$time_from_vaccination)
slopedat<-stackFits %>%
  dplyr::filter(time_from_vaccination%in%c(-50,-14,1,13,14,50)) %>%
  mutate(period=case_when(
    time_from_vaccination %in% c(-50,-14) ~ "pre_50_14",
    time_from_vaccination %in% c(1,13) ~ "post_1_13",
    time_from_vaccination %in% c(14,50) ~ "post_14_50")) %>%
  group_by(period,ind) %>%
  mutate(slope=diff(prob_some_contacts)/abs(diff(time_from_vaccination))) %>%
  dplyr::select(-prob_some_contacts,-time_from_vaccination) %>%
  mutate(obs=paste0(ind,period,slope)) %>%
  filter(!duplicated(obs)) %>%
  dplyr::select(-obs)

# compute differences in slopes between periods for each draw
slopedat_pp_1_13 <- as.data.frame(slopedat)
slopedat_pp_1_13 <- slopedat_pp_1_13 %>%
  dplyr::filter(period!="post_14_50") %>%
  group_by(ind) %>%
  dplyr::mutate(slope_diff=diff(slope) )%>%
  dplyr::select(ind,slope,slope_diff) %>%
  filter(!duplicated(ind)) %>%
  mutate(slope=slope+slope_diff)

slopedat_pp_14_50 <- as.data.frame(slopedat)
slopedat_pp_14_50 <- slopedat_pp_14_50 %>%
  dplyr::filter(period!="post_1_13") %>%
  group_by(ind) %>%
  dplyr::mutate(slope_diff=diff(slope)) %>%
  dplyr::select(ind,slope,slope_diff) %>%
  filter(!duplicated(ind)) %>%
  mutate(slope=slope+slope_diff)

slopedat_pp_pre <- as.data.frame(slopedat)
slopedat_pp_pre <- slopedat_pp_pre %>%
  dplyr::filter(period!="post_1_13") %>%
  dplyr::filter(period!="post_14_50") 
summary(slopedat_pp_pre$slope)

# compute medians and percentiles of differences over draws for pre/post periods
slope_pp_pre <- quantile(slopedat_pp_pre$slope, c(0.025,0.5,0.975))
slope_pp_pre
slope_pp_1_13 <- quantile(slopedat_pp_1_13$slope, c(0.025,0.5,0.975))
slope_pp_1_13
slope_pp_14_50 <- quantile(slopedat_pp_14_50$slope, c(0.025,0.5,0.975))
slope_pp_14_50
percentiles_slopediff_pp_1_13 <- quantile(slopedat_pp_1_13$slope_diff, c(0.025,0.5,0.975))
percentiles_slopediff_pp_14_50 <- quantile(slopedat_pp_14_50$slope_diff, c(0.025,0.5,0.975))


# add to a table
p1<-as.data.frame((percentiles_slopediff_pp_1_13))
p1a<-as.data.frame((slope_pp_1_13))
p1<-cbind("post_1_13",p1a,p1)
p2<-as.data.frame((percentiles_slopediff_pp_14_50))
p2a<-as.data.frame((slope_pp_14_50))
p2<-cbind("post_14_50",p2a,p2)
p3<-as.data.frame((slope_pp_pre))
slope_diffs<-cbind(p1,p2,p3)





