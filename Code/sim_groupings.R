library(nlme)
library(lme4)
library(MASS)
library(survey)
library(ClusterR)

h_tot = 200 # the number of clusters

# parameter for random intercepts
mu_alpha = 0; sd_alpha = 1 
mu_beta = 0; sd_beta = 1 
mu_kappa = 0; sd_kappa = 1

# parameter for fixed effect
alpha1 = -1; alpha2 = -1; alpha3 = -0.5
alpha4 = c(-2, 0, 2)

beta1 = 1; beta2 = -1; beta3 = 0.5
beta4 = c(-2, 0,  2)

kappa1 = -0.5; kappa2 = 1; kappa3 = -1
kappa4 = c(-2, 0, 2)


## parameters to control ##
q = l = m = 1 ## q,l,m = 1,2,3

n_h = as.integer(runif(h_tot, 5, 25)) # the number of subjects within each cluster
n_tot = sum(n_h) # total number of subjects
size = n_h # cluster size 
h_index = rep(1:h_tot, size)

alpha_h = rnorm(h_tot, mean = mu_alpha, sd = sd_alpha); alpha_h = rep(alpha_h, size) 
beta_h = rnorm(h_tot, mean = mu_beta, sd = sd_beta); beta_h = rep(beta_h, size)
kappa_h = rnorm(h_tot, mean = mu_kappa, sd = sd_kappa); kappa_h = rep(kappa_h, size)

## individual-level covariates
mu_X= rep(rnorm(h_tot, mean = 0, sd = 1), size); sd_X= 1
X = rnorm(n_tot, mean = mu_X, sd = sd_X)
Xbar.tmp = rep(0, h_tot)
count = 1
for(i in 1:h_tot){
  no.temp = size[i] 
  countnew = count + no.temp - 1
  Xbar.tmp[i] = mean(X[count:countnew])
  count=countnew+1
}
Xbar = rep(Xbar.tmp, size)
Xrel = X - Xbar

## cluster-level covariates
V.h = runif(h_tot, -1, 1); V = rep(V.h, size)
U.h = runif(h_tot, -2, 2); U = rep(U.h, size)


## generate the trt indicator
am_lin = Xbar*alpha1 + Xrel*alpha2 + V*alpha3 + (U-mean(U))*alpha4[q] + alpha_h 
ps_true = exp(am_lin)/(1+exp(am_lin))
ps_true = ps_true/0.7^-1 + 0.15
trt = rbinom(n_tot, size=1, prob=ps_true)

## generate obs outcome
trteffect_lin = kappa_h + Xbar*kappa1 + Xrel*kappa2 + V*kappa3 + (U-mean(U))^2*kappa4[m]
err_y = rnorm(n_tot, mean=0, sd=1)
y = trt*(trteffect_lin) + 
  beta_h + Xbar*beta1 + Xrel*beta2 + V*beta3 + (U-mean(U))*beta4[l] + err_y # Equation 19

## true ATE
ATE = mean(trteffect_lin)

## cluster-specific values
trt.prop = ATE.cluster = true.ps = rep(0, h_tot)
count = 1
for(i in 1:h_tot){
  no.temp=size[i] 
  countnew=count+no.temp-1 # individual index
  true.ps[i] = mean(ps_true[count:countnew])
  trt.prop[i] = mean(trt[count:countnew])
  ATE.cluster[i] = mean(trteffect_lin[count:countnew])
  count=countnew+1
}

## observed data
data_obs = data.frame(y = y, h_index = h_index, 
                      trt = trt, X = X, Xbar = Xbar, Xrel = Xrel, V = V, ps_true = ps_true)

## observed cluster data
data_obs_cluster = data.frame(h_index = 1:h_tot, Xbar = Xbar.tmp, V = V.h, U = U.h,
                              trt.prop = trt.prop, true.ps = true.ps, ATE.cluster = ATE.cluster)

## (1) grouping by similar (proportions) ##
cluster_dat =  as.data.frame(data_obs_cluster$trt.prop)
cluster_dat = center_scale(cluster_dat)
cm = Cluster_Medoids(cluster_dat, clusters = 10, swap_phase = TRUE)
data_obs_cluster$cluster.dummy = cm$clusters
dummy.num = 10; 
group.ps = rep(NA, n_tot)
ate.group = w.group = sesq.group.group = rep(NA, dummy.num)
delta_h_proportion = lambda_h_proportion =  rep(NA, h_tot)

for(d in 1:dummy.num){
  tmp.data = data_obs[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d), ]
  forktrt=(tmp.data$trt==1) 
  forkcon=(tmp.data$trt==0)
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glmer(trt ~ Xbar + Xrel + V + (1|h_index), data = tmp.data, family=binomial(link="logit"))
    tmp.ps = predict(tmp.reg)
    tmp.ps = exp(tmp.ps)/(1+exp(tmp.ps))
    
    group.ps[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d)] = tmp.ps
    
    tmp.data$ps = tmp.ps
    
    if(sum(is.na(tmp.ps)) == 0 & sum(tmp.ps == 0) == 0 & sum(tmp.ps == 1) == 0){
      ate.group[d] =  sum((1/tmp.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.ps[forktrt])-
        sum((1/(1-tmp.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.ps[forkcon]))     
      w.group[d] = sum(1/tmp.ps[forktrt])+sum(1/(1-tmp.ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$y)*(sum((1/tmp.ps[forktrt])^2)/(sum(1/tmp.ps[forktrt]))^2 + 
                                               sum((1/(1-tmp.ps[forkcon]))^2)/(sum(1/(1-tmp.ps[forkcon])))^2)
      
    }
    for(i in 1:h_tot){
      dum.index = which(data_obs_cluster$cluster.dummy == d)
      # proportion of the treated 
      delta_h_proportion[dum.index] =  data_obs_cluster$trt.prop[dum.index] - mean(tmp.data$trt)
      lambda_h_proportion[dum.index] = n_h[dum.index]*(delta_h_proportion[dum.index]*(1/mean(tmp.data$trt))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))^2*kappa4[m] + 
                                                  delta_h_proportion[dum.index]*(1/mean(tmp.data$trt) + 1/(1-mean(tmp.data$trt)))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))*beta4[l])
      
    }
  }
}
proportion.group = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.proportion.group= sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
lambda_proportion = sum(lambda_h_proportion, na.rm = TRUE) / n_tot

## (2) grouping by similar (proportions + observed covariates) ## 
cluster_dat =  as.data.frame(cbind(data_obs_cluster$Xbar, data_obs_cluster$V, data_obs_cluster$trt.prop))
cluster_dat = center_scale(cluster_dat)
cm = Cluster_Medoids(cluster_dat, clusters = 10, swap_phase = TRUE)
data_obs_cluster$cluster.dummy = cm$clusters
dummy.num = 10; 
group.ps = rep(NA, n_tot)
ate.group = w.group = sesq.group.group = rep(NA, dummy.num)
delta_h_both = lambda_h_both =  rep(NA, h_tot)

for(d in 1:dummy.num){
  tmp.data = data_obs[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d), ]
  forktrt=(tmp.data$trt==1) 
  forkcon=(tmp.data$trt==0) 
  
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glmer(trt ~ Xbar + Xrel + V + (1|h_index), data = tmp.data, family=binomial(link="logit"))
    tmp.ps = predict(tmp.reg)
    tmp.ps = exp(tmp.ps)/(1+exp(tmp.ps))
    group.ps[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d)] = tmp.ps
    tmp.data$ps = tmp.ps
    
    if(sum(is.na(tmp.ps)) == 0 & sum(tmp.ps == 0) == 0 & sum(tmp.ps == 1) == 0){
      ate.group[d] =  sum((1/tmp.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.ps[forktrt])-
        sum((1/(1-tmp.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.ps[forkcon]))     
      w.group[d] = sum(1/tmp.ps[forktrt])+sum(1/(1-tmp.ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$y)*(sum((1/tmp.ps[forktrt])^2)/(sum(1/tmp.ps[forktrt]))^2 + 
                                               sum((1/(1-tmp.ps[forkcon]))^2)/(sum(1/(1-tmp.ps[forkcon])))^2)
      
    } 
    for(i in 1:h_tot){
      dum.index = which(data_obs_cluster$cluster.dummy == d)
      # proportion of the treated 
      delta_h_both[dum.index] =  data_obs_cluster$trt.prop[dum.index] - mean(tmp.data$trt)
      lambda_h_both[dum.index] = n_h[dum.index]*(delta_h_both[dum.index]*(1/mean(tmp.data$trt))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))^2*kappa4[m] + 
                                                           delta_h_both[dum.index]*(1/mean(tmp.data$trt) + 1/(1-mean(tmp.data$trt)))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))*beta4[l])

    }
  }
}
proportion.covariate.group = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.proportion.covariate.group= sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
lambda_both = sum(lambda_h_both, na.rm = TRUE) / n_tot


## (3) grouping by similar observed covariates ## 
cluster_dat =  as.data.frame(cbind(data_obs_cluster$Xbar, data_obs_cluster$V))
cluster_dat = center_scale(cluster_dat)
cm = Cluster_Medoids(cluster_dat, clusters = 10, swap_phase = TRUE)
data_obs_cluster$cluster.dummy = cm$clusters
dummy.num = 10; 
group.ps = rep(NA, n_tot)
ate.group = w.group = sesq.group.group = rep(NA, dummy.num)
delta_h_cov = lambda_h_cov =  rep(NA, h_tot)

for(d in 1:dummy.num){
  tmp.data = data_obs[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d), ]
  forktrt=(tmp.data$trt==1) 
  forkcon=(tmp.data$trt==0) 
  
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glmer(trt ~ Xbar + Xrel + V + (1|h_index), data = tmp.data, family=binomial(link="logit"))
    tmp.ps = predict(tmp.reg)
    tmp.ps = exp(tmp.ps)/(1+exp(tmp.ps))
    
    group.ps[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d)] = tmp.ps
    tmp.data$ps = tmp.ps
    
    if(sum(is.na(tmp.ps)) == 0 & sum(tmp.ps == 0) == 0 & sum(tmp.ps == 1) == 0){
      ate.group[d] =  sum((1/tmp.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.ps[forktrt])-
        sum((1/(1-tmp.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.ps[forkcon]))     
      w.group[d] = sum(1/tmp.ps[forktrt])+sum(1/(1-tmp.ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$y)*(sum((1/tmp.ps[forktrt])^2)/(sum(1/tmp.ps[forktrt]))^2 + 
                                               sum((1/(1-tmp.ps[forkcon]))^2)/(sum(1/(1-tmp.ps[forkcon])))^2)
    }
    for(i in 1:h_tot){
      dum.index = which(data_obs_cluster$cluster.dummy == d)
      # proportion of the treated 
      delta_h_cov[dum.index] =  data_obs_cluster$trt.prop[dum.index] - mean(tmp.data$trt)
      lambda_h_cov[dum.index] = n_h[dum.index]*(delta_h_cov[dum.index]*(1/mean(tmp.data$trt))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))^2*kappa4[m] + 
                                                           delta_h_cov[dum.index]*(1/mean(tmp.data$trt) + 1/(1-mean(tmp.data$trt)))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))*beta4[l])

    }
    
  }
}
covariate.group = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.covariate.group= sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
lambda_cov = sum(lambda_h_cov, na.rm = TRUE) / n_tot

## (4) grouping randomly ## 
data_obs_cluster$cluster.dummy = sample(c(1:10), nrow(data_obs_cluster), replace = TRUE)
dummy.num = 10; 
group.ps = rep(NA, n_tot)
ate.group = w.group = sesq.group.group = rep(NA, dummy.num)
delta_h_random = lambda_h_random =  rep(NA, h_tot)

for(d in 1:dummy.num){
  tmp.data = data_obs[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d), ]
  forktrt=(tmp.data$trt==1) 
  forkcon=(tmp.data$trt==0) 
  
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glmer(trt ~ Xbar + Xrel + V + (1|h_index), data = tmp.data, family=binomial(link="logit"))
    tmp.ps = predict(tmp.reg)
    tmp.ps = exp(tmp.ps)/(1+exp(tmp.ps))
    group.ps[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d)] = tmp.ps
    tmp.data$ps = tmp.ps
    
    if(sum(is.na(tmp.ps)) == 0 & sum(tmp.ps == 0) == 0 & sum(tmp.ps == 1) == 0){
      ate.group[d] =  sum((1/tmp.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.ps[forktrt])-
        sum((1/(1-tmp.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.ps[forkcon]))     
      w.group[d] = sum(1/tmp.ps[forktrt])+sum(1/(1-tmp.ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$y)*(sum((1/tmp.ps[forktrt])^2)/(sum(1/tmp.ps[forktrt]))^2 + 
                                               sum((1/(1-tmp.ps[forkcon]))^2)/(sum(1/(1-tmp.ps[forkcon])))^2)
      
    } 
    for(i in 1:h_tot){
      dum.index = which(data_obs_cluster$cluster.dummy == d)
      # proportion of the treated 
      delta_h_random[dum.index] =  data_obs_cluster$trt.prop[dum.index] - mean(tmp.data$trt)
      lambda_h_random[dum.index] = n_h[dum.index]*(delta_h_random[dum.index]*(1/mean(tmp.data$trt))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))^2*kappa4[m] + 
                                                           delta_h_random[dum.index]*(1/mean(tmp.data$trt) + 1/(1-mean(tmp.data$trt)))*(data_obs_cluster$U[dum.index]-mean(data_obs_cluster$U))*beta4[l])

    }
  }
}
random.group = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.random.group= sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
lambda_random = sum(lambda_h_random, na.rm = TRUE) / n_tot

print(list(ATE = ATE, 
            # estimate
            random.group = random.group, proportion.group = proportion.group,
            proportion.covariate.group = proportion.covariate.group, covariate.group = covariate.group,
            # standard error
            se.random.group = se.random.group, se.proportion.group = se.proportion.group,
            se.proportion.covariate.group = se.proportion.covariate.group, se.covariate.group = se.covariate.group,
            # lambda (bias due to U)
            lambda_proportion = lambda_proportion, lambda_both = lambda_both,
            lambda_random = lambda_random, lambda_cov = lambda_cov
))