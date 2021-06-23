library(nlme)
library(lme4)
library(MASS)
library(flexmix)
library(ClusterR)

h_tot = 200 # the number of clusters

# parameter for random intercepts
mu_alpha = 0; sd_alpha = 1 
mu_beta = 0; sd_beta = 1 
mu_kappa = 0; sd_kappa = 1

# parameter for fixed effect
alpha1 = -1; alpha2 =  c(-2, -1, 0, 1, 2); alpha3 = c(-2, -1, 0, 1, 2);
alpha4 = 0.5

beta1 = 1; beta2 = -1; beta3 = 0.5
beta4 = 0.5

kappa1 = -0.5; kappa2 = 1; kappa3 = -1
kappa4 = 0.5

ATE = overlap.cl.count = ate.pooled.cluster = ate.pooled.group = ate.pooled.pooled =
ate.group.cluster= ate.group.pooled = ate.group.group = se.pooled.group = 
se.pooled.cluster = se.pooled.pooled = se.group.pooled= se.group.group = se.group.cluster = ate.pooled.pooled.fixed= 
se.pooled.pooled.fixed = c()


for(q in 1:5){
  
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
  
  ## cluster-level covariates (W and U are unobserved)
  V.h = runif(h_tot, -1, 1); V = rep(V.h, size)
  U.h = runif(h_tot, -2, 2); U = rep(U.h, size)
  

  ## generate the trt indicator
  am_lin = Xbar*alpha1 + Xrel*alpha2[q] + V*alpha3[q] + (U-mean(U))*alpha4 + alpha_h 
  ps_true = exp(am_lin)/(1+exp(am_lin))
  ps_true = ps_true/0.7^-1 + 0.15
  trt = rbinom(n_tot, size=1, prob=ps_true)
  
  ## generate obs outcome
  trteffect_lin = kappa_h + Xbar*kappa1 + Xrel*kappa2 + V*kappa3 + (U-mean(U))^2*kappa4
  err_y = rnorm(n_tot, mean=0, sd=1)
  y = trt*(trteffect_lin) + 
    beta_h + Xbar*beta1 + Xrel*beta2 + V*beta3 + (U-mean(U))*beta4 + err_y # Equation 19
  
  ## true ATE
  ATE[q] = mean(trteffect_lin)
  
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
  
  
  ## fully pooled propensity scores (random effect only)
  reg.pooled = tryCatch(glmer(trt ~ Xbar + Xrel + V + (1|h_index), data=data_obs, family=binomial), error=function(err) NA)
  
  if(!is.na(reg.pooled)){
    pooled.ps = predict(reg.pooled)
    pooled.ps = exp(pooled.ps)/(1+exp(pooled.ps))
    
    data_obs$pooled.ps = pooled.ps
    data_obs_cluster$pooled.ps = aggregate(pooled.ps~h_index, data = data_obs, mean)$pooled.ps
    
    ate.pooled.pooled[q] = sum(data_obs$trt*data_obs$y/data_obs$pooled.ps)/sum(data_obs$trt/data_obs$pooled.ps) - 
      sum((1-data_obs$trt)*data_obs$y/(1-data_obs$pooled.ps))/sum((1-data_obs$trt)/(1-data_obs$pooled.ps))
    
    sigsq1 = var(data_obs$y[data_obs$trt==1])
    sigsq0 = var(data_obs$y[data_obs$trt==0])
    se.pooled.pooled[q] = sqrt(sigsq1*sum((1/data_obs$pooled.ps[data_obs$trt==1])^2)/(sum(1/data_obs$pooled.ps[data_obs$trt==1]))^2 + 
                                 sigsq0*sum((1/data_obs$pooled.ps[data_obs$trt==0])^2)/(sum(1/data_obs$pooled.ps[data_obs$trt==0]))^2)
    
  }
  ## fully pooled propensity scores (fixed effect)
  fixed.pooled = tryCatch(glm(trt ~ Xbar + Xrel + V + as.factor(h_index), data=data_obs, family=binomial), error=function(err) NA); 
  
  if(sum(is.na(predict(fixed.pooled))) > 0 | is.na(fixed.pooled)){
    ate.pooled.pooled.fixed[q] = NA
    se.pooled.pooled.fixed[q] = NA
  }else{
    fixed.ps = predict(fixed.pooled)
    fixed.ps = exp(fixed.ps)/(1+exp(fixed.ps))
    
    data_obs$fixed.ps = fixed.ps
    #data_obs_cluster$fixed.ps = aggregate(fixed.ps~h_index, data = data_obs, mean)$fixed.ps
    
    ate.pooled.pooled.fixed[q] = sum(data_obs$trt*data_obs$y/data_obs$fixed.ps)/sum(data_obs$trt/data_obs$fixed.ps) - 
      sum((1-data_obs$trt)*data_obs$y/(1-data_obs$fixed.ps))/sum((1-data_obs$trt)/(1-data_obs$fixed.ps))
    
    sigsq1 = var(data_obs$y[data_obs$trt==1])
    sigsq0 = var(data_obs$y[data_obs$trt==0])
    se.pooled.pooled.fixed[q] = sqrt(sigsq1*sum((1/data_obs$fixed.ps[data_obs$trt==1])^2)/(sum(1/data_obs$fixed.ps[data_obs$trt==1]))^2 + 
                                       sigsq0*sum((1/data_obs$fixed.ps[data_obs$trt==0])^2)/(sum(1/data_obs$fixed.ps[data_obs$trt==0]))^2)
  }
  
  ## grouping
  #cluster_dat =  as.data.frame(data_obs_cluster$trt.prop)
  #cluster_dat = center_scale(cluster_dat)
  cm = kmeans(data_obs_cluster$trt.prop, 10)
  data_obs_cluster$cluster.dummy = cm$cluster
  dummy.num = 10
  ate.cluster.h = sesq.group.cluster = w.cluster.h = rep(NA, h_tot)
  ate.h = sesq.pooled.cluster = w.h = rep(NA, h_tot)
  
  data_obs$group.ps = rep(NA, nrow(data_obs))
  ate.group = w.group = rep(NA, dummy.num)
  ate.pool = w.pool = rep(NA, dummy.num)
  sesq.pooled.group = sesq.group.pooled = sesq.group.group = rep(NA, dummy.num)
  
  for(d in 1:dummy.num){
    tmp.data = data_obs[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d), ]
    forktrt=(tmp.data$trt==1) 
    forkcon=(tmp.data$trt==0) 
    
    if(sum(forktrt)!=0 & sum(forkcon)!=0){
      ## (full, group) : uses fully pooled propensity scores
      ate.pool[d] =  sum((1/tmp.data$pooled.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.data$pooled.ps[forktrt])-
        sum((1/(1-tmp.data$pooled.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.data$pooled.ps[forkcon]))   
      w.pool[d] =  sum(1/tmp.data$pooled.ps[forktrt])+sum(1/(1-tmp.data$pooled.ps[forkcon]))
      sesq.pooled.group[d] = var(tmp.data$y)*(sum((1/tmp.data$pooled.ps[forktrt])^2)/(sum(1/tmp.data$pooled.ps[forktrt]))^2 + 
                                                sum((1/(1-tmp.data$pooled.ps[forkcon]))^2)/(sum(1/(1-tmp.data$pooled.ps[forkcon])))^2)
      
      ## estimate partially pooled propensity scores
      tmp.reg = tryCatch(glmer(trt ~ Xbar + Xrel + V + (1|h_index), data = tmp.data, family=binomial(link="logit")), error=function(err) NA); 
      
      if(!is.na(tmp.reg)){
        tmp.ps = predict(tmp.reg)
        tmp.ps = exp(tmp.ps)/(1+exp(tmp.ps))
        data_obs$group.ps[data_obs$h_index %in% which(data_obs_cluster$cluster.dummy == d)] = tmp.ps         
        tmp.data$group.ps = tmp.ps
        
        if(sum(is.na(tmp.ps)) == 0 & sum(tmp.ps == 0) == 0 & sum(tmp.ps == 1) == 0){
          
          ate.group[d] =  sum((1/tmp.ps[forktrt])*tmp.data$y[forktrt])/sum(1/tmp.ps[forktrt])-
            sum((1/(1-tmp.ps[forkcon]))*tmp.data$y[forkcon])/sum(1/(1-tmp.ps[forkcon]))     
          w.group[d] = sum(1/tmp.ps[forktrt])+sum(1/(1-tmp.ps[forkcon]))
          sesq.group.group[d] = var(tmp.data$y)*(sum((1/tmp.ps[forktrt])^2)/(sum(1/tmp.ps[forktrt]))^2 + 
                                                   sum((1/(1-tmp.ps[forkcon]))^2)/(sum(1/(1-tmp.ps[forkcon])))^2)
          num.cluster = length(tmp.data$h_index)
          
          # for each cluster in the group
          for(k in unique(tmp.data$h_index)){
            tmp.cluster.data = tmp.data[tmp.data$h_index == k, ]
            tmp.forktrt=(tmp.cluster.data$trt==1) # index of treated individual within cluster
            tmp.forkcon=(tmp.cluster.data$trt==0) # index of control individual within cluster
            
            if(sum(tmp.forktrt)!=0 & sum(tmp.forkcon)!=0) {
              y.temp = tmp.cluster.data$y
              ## with fully pooled propensity scores
              w.temp.re = tmp.cluster.data$pooled.ps     
              ate.h[k] =  sum((1/w.temp.re[tmp.forktrt])*y.temp[tmp.forktrt])/sum(1/w.temp.re[tmp.forktrt])-sum((1/(1-w.temp.re[tmp.forkcon]))*y.temp[tmp.forkcon])/sum(1/(1-w.temp.re[tmp.forkcon]))     
              w.h[k] = sum(1/w.temp.re[tmp.forktrt])+sum(1/(1-w.temp.re[tmp.forkcon]))
              sesq.pooled.cluster[k] = var(y.temp)*(sum((1/w.temp.re[tmp.forktrt])^2)/(sum(1/w.temp.re[tmp.forktrt]))^2 + 
                                                      sum((1/(1-w.temp.re[tmp.forkcon]))^2)/(sum(1/(1-w.temp.re[tmp.forkcon])))^2)
              
              ## with partially pooled propensity scores
              group.ps = tmp.cluster.data$group.ps
              ate.cluster.h[k] = sum((1/group.ps[tmp.forktrt])*y.temp[tmp.forktrt])/sum(1/group.ps[tmp.forktrt]) -
                sum((1/(1-group.ps[tmp.forkcon]))*y.temp[tmp.forkcon])/sum(1/(1-group.ps[tmp.forkcon]))
              sesq.group.cluster[k] = var(y.temp)*(sum((1/group.ps[tmp.forktrt])^2)/(sum(1/group.ps[tmp.forktrt]))^2+
                                                     sum((1/(1-group.ps[tmp.forkcon]))^2)/(sum(1/(1-group.ps[tmp.forkcon])))^2)
              w.cluster.h[k] = sum(1/group.ps[tmp.forktrt])+sum(1/(1-group.ps[tmp.forkcon]))
              
            }
          }
        }
        
      }
    }
  }
  
  ## 2. (full, group)
  ate.pooled.group[q] = sum(ate.pool*w.pool, na.rm = TRUE)/sum(w.pool, na.rm = TRUE)
  se.pooled.group[q] = sqrt(sum(w.pool^2*sesq.pooled.group, na.rm = TRUE)) / sum(w.pool, na.rm = TRUE)
  
  ## 3. (full, cluster)
  overlap.cl.count[q] = sum(!is.na(ate.h)) # number of overlapping clusters in simulation.
  ate.pooled.cluster[q] = sum(ate.h*w.h, na.rm = TRUE)/sum(w.h, na.rm = TRUE)
  se.pooled.cluster[q] = sqrt(sum(w.h^2*sesq.pooled.cluster, na.rm = TRUE)) / sum(w.h, na.rm = TRUE)
  
  ## 4. (group, full)
  ate.group.pooled[q] = sum(data_obs$trt*data_obs$y/data_obs$group.ps, na.rm = TRUE)/sum(data_obs$trt/data_obs$group.ps, na.rm = TRUE) - 
    sum((1-data_obs$trt)*data_obs$y/(1-data_obs$group.ps), na.rm = TRUE)/sum((1-data_obs$trt)/(1-data_obs$group.ps), na.rm = TRUE)
  se.group.pooled[q] = sqrt(sigsq1*sum((1/data_obs$group.ps[data_obs$trt==1])^2, na.rm = TRUE)/(sum(1/data_obs$group.ps[data_obs$trt==1], na.rm = TRUE))^2 + 
                              sigsq0*sum((1/data_obs$group.ps[data_obs$trt==0])^2, na.rm = TRUE)/(sum(1/data_obs$group.ps[data_obs$trt==0], na.rm = TRUE))^2)
  
  ## 5. (group, group)
  ate.group.group[q] = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
  se.group.group[q] = sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
  
  ## 6. (group, cluster)
  ate.group.cluster[q] = sum(ate.cluster.h*w.cluster.h, na.rm = TRUE)/sum(w.cluster.h, na.rm = TRUE)
  se.group.cluster[q] = sqrt(sum(w.cluster.h^2*sesq.group.cluster, na.rm = TRUE)) / sum(w.cluster.h, na.rm = TRUE)
  
}


results = list(ATE = ATE, overlap.cl.count = overlap.cl.count,
            ate.pooled.cluster = ate.pooled.cluster, ate.pooled.group = ate.pooled.group, ate.pooled.pooled = ate.pooled.pooled, 
            ate.group.cluster = ate.group.cluster, ate.group.pooled = ate.group.pooled, ate.group.group = ate.group.group, 
            se.pooled.cluster = se.pooled.cluster, se.pooled.group = se.pooled.group, se.pooled.pooled = se.pooled.pooled, 
            se.group.pooled = se.group.pooled, se.group.group = se.group.group, se.group.cluster = se.group.cluster,
            ate.pooled.pooled.fixed = ate.pooled.pooled.fixed, se.pooled.pooled.fixed = se.pooled.pooled.fixed)

