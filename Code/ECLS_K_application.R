## download libraries and data
library(readstata13)
library(lme4)
library(foreign)
library(MatchIt)
library(survey)
library(data.table)
library(Hmisc)
library(xtable)
library(tableone)
library(Matching)
library(reshape2)
library(ggplot2)
library(ClusterR)
library(plotrix)
dat = read.dta13("~/Dropbox/2019 postdoc/ECLS/Data/childK8p.dta")
# data is publicly available at https://nces.ed.gov/ecls/dataproducts.asp.

## auxiliary functions
findU = function(x) quantile(x,probs=0.975, na.rm = TRUE)
findL = function(x) quantile(x,probs=0.025, na.rm = TRUE)

## extract the variable of interest from the raw data
data.1 = dat
VarsToKeep = c(
  # ID for school and child
  "CHILDID",
  "S1_ID"  , 
  "S2_ID"  ,
  
  # treatment
  "P1PRIMPK", # P1: primary type nonparental care pre-K 
  
  # outcome
  "C1R4MSCL", # Kindergarten math score
  
  # individual-level covariates
  "GENDER" , # male 
  "WKBLACK", # black
  "WKHISP", # hispanic
  "C1HEIGHT",  # height (inches)
  "C1WEIGHT", # weight (pounds)
  "P1AGEENT",  # age at the kindergarten entry
  "C2SCREEN", # speak non-English language at home  
  "P1HPARNT", # biological mother
  "P1HFAMIL", # Family type
  "C1CMOTOR", # motor skill
  "R1_KAGE", # child assessment age (month)
  "WKSESL" ,# socioeconomic status
  "P1HIG_1", # paraent education - continuous
  "WKINCOME", # weekly income (impuated)
  
  # cluster-level covariates
  "CREGION", # census region in sample frame
  "KURBAN_R" # location type in base year sample frame
) 


data.t        = subset(data.1, select=VarsToKeep)
# only consider center-based program vs. parental case
data.t$P1PRIMPK = ifelse(data.t$P1PRIMPK == 6, 1, ifelse(data.t$P1PRIMPK == 0, 0, NA))
data.t$sample = complete.cases(data.t)
data.t        = subset(data.t, sample==T)
data.t        = subset(data.t, data.t$C1R4MSCL >=0)
data.t        = subset(data.t, data.t$WKBLACK >=0)
data.t        = subset(data.t, data.t$WKHISP >=0)
data.t        = subset(data.t, data.t$C1HEIGHT >=0)
data.t        = subset(data.t, data.t$C1WEIGHT >=0)
data.t        = subset(data.t, data.t$P1AGEENT >=0)
data.t        = subset(data.t, data.t$P1HFAMIL >=0)
data.t        = subset(data.t, data.t$C1CMOTOR >=0)
data.t        = subset(data.t, data.t$R1_KAGE >=0)
data.t        = subset(data.t, data.t$P1HIG_1 >=0)

# define a binary treatment assignment
data.t$TREAT  = ifelse(data.t$P1PRIMPK == 1,1,0)

data.t$FEMALE = ifelse(data.t$GENDER   == 2,1,0)
data.t$BLACK  = ifelse(data.t$WKBLACK  == 1,1,0)
data.t$HISPANIC  = ifelse(data.t$WKHISP  == 1,1,0)
data.t$ENGLISHHOME  = ifelse(data.t$C2SCREEN  == 2,1,0)
data.t$TWOPAR = ifelse(data.t$P1HFAMIL==1|data.t$P1HFAMIL==2,1,0) # two parents
data.t$FAMTYPE = as.factor(data.t$P1HFAMIL)
data.t$BIOMATHER = ifelse(data.t$P1HPARNT==1 | data.t$P1HPARNT== 2 | 
                            data.t$P1HPARNT==4 ,1,0) # at linest one biological
data.t$CREGION = as.factor(data.t$CREGION)
data.t$KURBAN_R = as.factor(data.t$KURBAN_R)

# propensity score model at Dong et al. 
formula_PS_Dong =  TREAT  ~ FEMALE + BLACK + HISPANIC + CREGION + 
  TWOPAR + FAMTYPE + BIOMATHER + 
  ENGLISHHOME + C1CMOTOR + C1HEIGHT + C1WEIGHT + R1_KAGE + 
  WKSESL + P1HIG_1 + WKINCOME + KURBAN_R 

# propensity score with three individual-level covariates omitted: WKSESL, P1HIG_1, WKINCOME
formula_PS =  TREAT  ~ FEMALE + BLACK + HISPANIC + CREGION + 
  TWOPAR + FAMTYPE + BIOMATHER + 
  ENGLISHHOME + C1CMOTOR + C1HEIGHT + C1WEIGHT + R1_KAGE  + 
  KURBAN_R 

# propensity score without cluster-level covariates:  CREGION, KURBAN_R, WKSESL, P1HIG_1, WKINCOME
formula_PS_noregion =  TREAT  ~ FEMALE + BLACK + HISPANIC +  
  TWOPAR + FAMTYPE + BIOMATHER + 
  ENGLISHHOME + C1CMOTOR + C1HEIGHT + C1WEIGHT + R1_KAGE 

# consider clusters having at least one treated and control
subindex = aggregate(TREAT ~ S2_ID, data = data.t, mean)$S2_ID[aggregate(TREAT ~ S2_ID, data = data.t, mean)$TREAT > 0 &
                                                                 aggregate(TREAT ~ S2_ID, data = data.t, mean)$TREAT < 1]
subdat = data.t[data.t$S2_ID %in% subindex, ]
subdat = subdat[subdat$S2_ID != "",]

# proportion of treated units 
proportion = aggregate(TREAT ~ S2_ID, data = subdat, mean)
cluster_dat =  data.frame(proportion$TREAT)
cluster_dat = center_scale(cluster_dat)
# use PAM algorithm 
cm = Cluster_Medoids(cluster_dat, clusters = 10, distance_metric = 'mahalanobis', 
                     swap_phase = TRUE)
table(cm$clusters)
proportion$class = cm$clusters; dummy = rep(NA, nrow(subdat))
for(i in 1:nrow(proportion)){
  dummy[subdat$S2_ID %in% proportion$S2_ID[i]] = proportion$class[i]
}
subdat$dummy = dummy
pdf("Figure/region.pdf",  width = 12, height = 6)
par(mfrow = c(1,1),   mar = c(7,7,5,5),  cex.lab = 2, 
    cex.main = 2, cex.axis = 1.5, tcl = 0.5,  omi = c(0,0,0,0))
plot(mat$Proportion, mat$Northeast, ylim = c(0,0.6),  xlim = c(0,1),
     pch = 18, col = "dodgerblue", mgp = c(5,2,0), type = "b", cex = 1.5,
     ylab = "Proportion of census region", xlab = "Average proportion of treated units in group (G = 10)", lwd = 2)
points(mat$Proportion, mat$Midwest,  type = "b", pch = 19, col = "brown",  lwd = 2)
points(mat$Proportion, mat$South,  type = "b", pch = 17, col = "darkgoldenrod",  lwd = 2)
points(mat$Proportion, mat$West,  type = "b", pch = 15, col = "forestgreen",  lwd = 2)
legend("topright", c("Northeast", "Midwest", "South", "West"), bty = "n",
       pch = c(18, 19, 17, 15), col = c("dodgerblue", "brown", "darkgoldenrod", "forestgreen"),
       pt.cex = c(1.5, 1, 1, 1))
dev.off()

## generate bootstrap samples
#data_bs = list()
#for(b in 1:1000){
#  tmp = NULL
#  # sampling with replacement of clusters
#  h_id_bs=sample(1:length(unique(subdat$S2_ID)), length(unique(subdat$S2_ID)), replace=T) 
#  for(i in 1:length(unique(subdat$S2_ID))){
#    tmp = rbind(tmp, subdat[subdat$S2_ID %in% unique(subdat$S2_ID)[h_id_bs[i]], ])
#  }
#  data_bs[[b]] = tmp
#}

## unweighted IPW 
ate.pooled.pooled.noweight = mean(subdat$C1R4MSCL[subdat$TREAT == 1])- mean(subdat$C1R4MSCL[subdat$TREAT == 0])  
ate.pooled.pooled.boot.noweight = c()
sigsq1 = var(subdat$C1R4MSCL[subdat$TREAT==1])
sigsq0 = var(subdat$C1R4MSCL[subdat$TREAT==0])
se.noweight =  sqrt(sigsq1/sum(subdat$TREAT==1)) + sqrt(sigsq0/sum(subdat$TREAT==0))

## fully pooled propensity scores (with U)
reg = glm(formula_PS, data = subdat, family=binomial(link="logit"))
ps = predict(reg)
ps = exp(ps)/(1+exp(ps))
forktrt=(subdat$TREAT==1) 
forkcon=(subdat$TREAT==0)
ate.pooled.pooled =  sum((1/ps[forktrt])*subdat$C1R4MSCL[forktrt])/sum(1/ps[forktrt])-
  sum((1/(1-ps[forkcon]))*subdat$C1R4MSCL[forkcon])/sum(1/(1-ps[forkcon]))
poolweight = (subdat$TREAT == 1)/ps + (subdat$TREAT == 0)/(1-ps) 
poolweight.subdat = svydesign(ids = ~ 1, data = subdat, weights = ~ poolweight)

## fully pooled propensity scores (without U)
reg = glm(formula_PS_noregion, data = subdat, family=binomial(link="logit"))
ps = predict(reg)
ps = exp(ps)/(1+exp(ps))
forktrt=(subdat$TREAT==1) 
forkcon=(subdat$TREAT==0)
ate.pooled.pooled.noregion =  sum((1/ps[forktrt])*subdat$C1R4MSCL[forktrt])/sum(1/ps[forktrt])-
  sum((1/(1-ps[forkcon]))*subdat$C1R4MSCL[forkcon])/sum(1/(1-ps[forkcon]))     
poolweight.noregion = (subdat$TREAT == 1)/ps + (subdat$TREAT == 0)/(1-ps) 
poolweight.subdat.noregion = svydesign(ids = ~ 1, data = subdat, weights = ~ poolweight.noregion)


## partially pooled propensity scores (with U)
dummy.num = 10 # the number of groups
ate.group = w.group = rep(NA, dummy.num)
group.ps = rep(NA, nrow(subdat)) # partially pooled propensity scores
ate.group = sesq.group.group = rep(NA, dummy.num) 
for(d in 1:dummy.num){
  tmp.data = subdat[subdat$dummy == d, ]
  forktrt=(tmp.data$TREAT==1) # index of treated individual within cluster
  forkcon=(tmp.data$TREAT==0) # index of control individual within cluster
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glm(formula_PS, data = tmp.data, family=binomial(link="logit"))
    ps=predict(tmp.reg)
    ps=exp(ps)/(1+exp(ps))
    group.ps[subdat$dummy == d] = ps
    if(sum(is.na(ps)) == 0 & sum(ps == 0) == 0 & sum(ps == 1) == 0){
      ate.group[d] =  sum((1/ps[forktrt])*tmp.data$C1R4MSCL[forktrt])/sum(1/ps[forktrt])-
        sum((1/(1-ps[forkcon]))*tmp.data$C1R4MSCL[forkcon])/sum(1/(1-ps[forkcon]))     
      w.group[d] = sum(1/ps[forktrt])+sum(1/(1-ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$C1R4MSCL)*(sum((1/ps[forktrt])^2)/(sum(1/ps[forktrt]))^2 + 
                                                      sum((1/(1-ps[forkcon]))^2)/(sum(1/(1-ps[forkcon])))^2)
    } 
  }
}
ate.group.group = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.group.group = sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
groupweight = (subdat$TREAT == 1)/group.ps + (subdat$TREAT == 0)/(1-group.ps) 
groupweight.subdat = svydesign(ids = ~ 1, data = subdat, weights = ~ groupweight)


## partially pooled propensity scores without two cluster-level covariates (without U)
dummy.num = 10
ate.group = w.group = rep(NA, dummy.num)
group.ps = rep(NA, nrow(subdat))
ate.group = sesq.group.group = rep(NA, dummy.num)
for(d in 1:dummy.num){
  tmp.data = subdat[subdat$dummy == d, ]
  forktrt=(tmp.data$TREAT==1) 
  forkcon=(tmp.data$TREAT==0) 
  if(sum(forktrt)!=0 & sum(forkcon)!=0){
    tmp.reg = glm(formula_PS_noregion, data = tmp.data, family=binomial(link="logit"))
    ps=predict(tmp.reg)
    ps=exp(ps)/(1+exp(ps))
    group.ps[subdat$dummy == d] = ps
    if(sum(is.na(ps)) == 0 & sum(ps == 0) == 0 & sum(ps == 1) == 0){
      ate.group[d] =  sum((1/ps[forktrt])*tmp.data$C1R4MSCL[forktrt])/sum(1/ps[forktrt])-
        sum((1/(1-ps[forkcon]))*tmp.data$C1R4MSCL[forkcon])/sum(1/(1-ps[forkcon]))     
      w.group[d] = sum(1/ps[forktrt])+sum(1/(1-ps[forkcon]))
      sesq.group.group[d] = var(tmp.data$C1R4MSCL)*(sum((1/ps[forktrt])^2)/(sum(1/ps[forktrt]))^2 + 
                                                      sum((1/(1-ps[forkcon]))^2)/(sum(1/(1-ps[forkcon])))^2)
    } 
  }
}
ate.group.group.noregion = sum(ate.group*w.group, na.rm = TRUE)/sum(w.group, na.rm = TRUE)
se.group.group.noregion = sqrt(sum(w.group^2*sesq.group.group, na.rm = TRUE)) / sum(w.group, na.rm = TRUE)
groupweight.noregion = (subdat$TREAT == 1)/group.ps + (subdat$TREAT == 0)/(1-group.ps) 
groupweight.subdat.noregion = svydesign(ids = ~ 1, data = subdat, weights = ~ groupweight.noregion)


## You can replicate above procedures with 1000 bootstraped samples, data_bs.
## all the bootstrapped results are saved at "Data/bootresult.RData"
load("Data/bootresult.RData")
ate.pooled.pooled.boot.noweight = bootresult$ate.pooled.pooled.boot.noweight
ate.group.group.boot = bootresult$ate.group.group.boot
ate.pooled.pooled.boot = bootresult$ate.pooled.pooled.boot
ate.group.group.boot.noregion = bootresult$ate.group.group.boot.noregion
ate.pooled.pooled.boot.noregion = bootresult$ate.pooled.pooled.boot.noregion
# calculate empirical confidence intervals
noweight.ci = c(findL(ate.pooled.pooled.boot.noweight), findU(ate.pooled.pooled.boot.noweight))
group.ci = c(findL(ate.group.group.boot), findU(ate.group.group.boot))
pooled.ci = c(findL(ate.pooled.pooled.boot), findU(ate.pooled.pooled.boot))
group.noregion.ci = c(findL(ate.group.group.boot.noregion), findU(ate.group.group.boot.noregion))
pooled.noregion.ci = c(findL(ate.pooled.pooled.boot.noregion), findU(ate.pooled.pooled.boot.noregion))

## plot for estimation results
mat = matrix(NA, 5, 4)
rownames(mat) = c("Not weighted", "(full, full)", "(group, group)",  "(full, full)", "(group, group)")
colnames(mat) = c("Estimate", "SE", "95% CI (bootstrap)", "95% CI (bootstrap)")
mat[,1] =  c(ate.pooled.pooled.noweight, ate.pooled.pooled,
             ate.group.group, ate.pooled.pooled.noregion,
             ate.group.group.noregion)
mat[,2] = c(se.noweight, se.pooled.pooled, se.group.group, se.pooled.pooled.noregion, se.group.group.noregion)
mat[1,c(3,4)] = noweight.ci 
mat[2,c(3,4)] = pooled.ci
mat[3,c(3,4)] = group.ci
mat[4,c(3,4)] = pooled.noregion.ci
mat[5,c(3,4)] = group.noregion.ci
print(xtable(mat, digits = 2))
tmp.result = c(ate.pooled.pooled.boot.noweight, ate.pooled.pooled.boot,
               ate.group.group.boot, ate.pooled.pooled.boot.noregion, ate.group.group.boot.noregion)
types = rep(c("Unweighted effect", "Fully pooled PS", "Partially PS", "Fully pooled PS (without U)", "Fully pooled PS (without U)") , each = 1000)
types = factor(types, levels = c("Unweighted effect", "Pooled PS", "Group-sourced PS", "Pooled PS (without U)", "Group-sourced PS (without U)"))

tmp = data.frame(tmp.result = tmp.result, types = types)
pdf("Figure/results.pdf",  width = 22, height = 10)
par(mfrow = c(1,1),   mar = c(5,7,5,5),  cex.lab = 3, 
    cex.main = 3.0, cex.axis = 1.5, tcl = 0.5,  omi = c(0,0.5,0,0))
plotCI(c(1:5), c(ate.pooled.pooled.noweight, ate.pooled.pooled,
                 ate.group.group, ate.pooled.pooled.noregion,
                 ate.group.group.noregion),
       li = c(noweight.ci[1], pooled.ci[1], group.ci[1], pooled.noregion.ci[1], group.noregion.ci[1]),
       ui = c(noweight.ci[2], pooled.ci[2], group.ci[2], pooled.noregion.ci[2], group.noregion.ci[2]),
       xlab = "", ylab = "Estimated ATE", 
       col= c(rep("black",3),rep("dodgerblue",2)),
       main= "", pch = 19, cex= 3,
       mgp = c(5,2,0),
       xaxt='n', yaxt = 'n',  ylim = c(-0.5, 5), lwd = 3,
       xlim = c(0.5, 5.5))
text(x = c(1:5)+0.2, y = c(ate.pooled.pooled.noweight, ate.pooled.pooled,
                           ate.group.group, ate.pooled.pooled.noregion,
                           ate.group.group.noregion), labels = formatC(c(ate.pooled.pooled.noweight, ate.pooled.pooled,
                                                                         ate.group.group, ate.pooled.pooled.noregion,
                                                                         ate.group.group.noregion), format = "f", digits = 2),
     cex = 2)
axis(2, at = seq(-2, 5, 1), 
     labels = seq(-2, 5, 1), las = 2,
     col = "black", cex.axis = 2)
axis(1, at = c(1:5),  labels =  c("Unweighted effect", "Fully pooled PS", "Partially pooled PS", "Fully pooled PS (without U)", "Partially pooled PS (without U)"),
     col = "black", cex = 1.3)
abline(h = 0, col = "red")
dev.off()


## covariate balance when unweighted
vars =  c("C1R4MSCL", "FEMALE", "BLACK", "HISPANIC", 
                   "TWOPAR", "FAMTYPE", "BIOMATHER","ENGLISHHOME", 
                   "C1CMOTOR", "C1HEIGHT", "C1WEIGHT", "R1_KAGE", 
                   "CREGION", "KURBAN_R")

# construct a table
tabUnmatched <- CreateTableOne(vars = vars, strata = "P1PRIMPK", data = subdat, test = FALSE)
# print the table with SMD
print(tabUnmatched, smd = TRUE)

## covariate balance when weighted by fully pooled propensity scores
vars =  c("C1R4MSCL", "FEMALE", "BLACK", "HISPANIC", 
          "TWOPAR", "FAMTYPE", "BIOMATHER","ENGLISHHOME", 
          "C1CMOTOR", "C1HEIGHT", "C1WEIGHT", "R1_KAGE", 
          "CREGION", "KURBAN_R")
# construct a table
tabPooled = svyCreateTableOne(vars = vars, strata = "TREAT", data = poolweight.subdat.noregion, test = FALSE)
# print the table with SMD
print(tabPooled, smd = TRUE)

## covariate balance when weighted by partially pooled propensity scores
vars <-  c("C1R4MSCL", "FEMALE", "BLACK", "HISPANIC", 
           "TWOPAR", "FAMTYPE", "BIOMATHER",  
           "ENGLISHHOME", "C1CMOTOR", "C1HEIGHT", "C1WEIGHT", "R1_KAGE", 
           "CREGION", "KURBAN_R")
# construct a table
tabGroup = svyCreateTableOne(vars = vars, strata = "TREAT", data = groupweight.subdat.noregion, test = FALSE)
# print the table with SMD
print(tabGroup, smd = TRUE)

## combine tables (unweighted/weighted by fully pooled PS/weighted by partially pooled PS)
resCombo = cbind(print(tabUnmatched, printToggle = FALSE, smd = TRUE),
                  print(tabPooled,  printToggle = FALSE, smd = TRUE),
                  print(tabGroup,   printToggle = FALSE, smd = TRUE))
resCombo = resCombo
colnames(resCombo) = c("Unweighted", "Pooled PS", "Group PS")
print(resCombo, quote = FALSE)
print(xtable(resCombo))