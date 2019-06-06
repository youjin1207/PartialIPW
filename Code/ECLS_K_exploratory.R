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
dat <- read.dta13("~/Dropbox/2019 postdoc/ECLS/Data/childK8p.dta")
# data is publicly available at https://nces.ed.gov/ecls/dataproducts.asp.

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

# consider clusters having at least one treated and control
subindex = aggregate(TREAT ~ S2_ID, data = data.t, mean)$S2_ID[aggregate(TREAT ~ S2_ID, data = data.t, mean)$TREAT > 0 &
                                                                 aggregate(TREAT ~ S2_ID, data = data.t, mean)$TREAT < 1]
subdat = data.t[data.t$S2_ID %in% subindex, ]
subdat = subdat[subdat$S2_ID != "",]

##
region.mat = cbind(aggregate(TREAT ~ CREGION, data = subdat, mean), aggregate(C1R4MSCL ~ CREGION, data = subdat, mean))
rural.mat = cbind(aggregate(TREAT ~ KURBAN_R, data = subdat, mean), aggregate(C1R4MSCL ~ KURBAN_R, data = subdat, mean))
print(xtable(region.mat, digits = 2, row.names = FALSE))

prop = aggregate(TREAT ~ S2_ID, data = subdat, mean)$TREAT
prop = prop[-1]
region = c()
for(i in 1:length(prop)){
  region[i] = subdat$CREGION[subdat$S2_ID == aggregate(TREAT ~ S2_ID, data = subdat, mean)$S2_ID[i]][1]
}

prop.cate = ifelse(prop <= 0.1, 1, ifelse(prop <= 0.2, 2,
                                          ifelse(prop <= 0.3, 3, ifelse(prop <= 0.4, 4, ifelse(prop <= 0.5, 5, ifelse(prop <= 0.6, 6, ifelse(prop <= 0.7, 7, ifelse(prop <= 0.8, 8, ifelse(prop <= 0.9, 9, 10)))))))))

counts = table(region, prop.cate)
adj.counts = counts / matrix(rep(colSums(counts),each = 4), 4, 10)
colnames(counts) = seq(0.1, 1.0, 0.1)
colnames(adj.counts) = seq(0.1, 1, 0.1)
pdf("Figure/preliminary.pdf",  width = 17, height = 12)
par(mfrow = c(1,1),   mar = c(7,7,7,15),  cex.lab = 2, 
    cex.main = 2.5, cex.axis = 2, tcl = 0.5,  omi = c(0.5,0,0,0), xpd = TRUE)
barplot(adj.counts, main="",
        col=c("dodgerblue","brown", "darkgoldenrod", "forestgreen"), xlab = "",
        names.arg = c("< 0.1", "< 0.2", "< 0.3", "< 0.4",
                      "< 0.5", "< 0.6", "< 0.7", "< 0.8", "< 0.9", "< 1.0"))
legend(12, 0.5, c("Northeast", "Midwest", "South", "West"), 
       col = c("dodgerblue","brown", "darkgoldenrod", "forestgreen"), pch = 15, cex = 2.5, bty = "n",
       pt.cex = 3)
text(x = 6, y = -0.1, "Proportion of students having pre-school center-based program", cex = 2,
     xpd = TRUE)
text(x = 7,  y = 1.1, "Distribution of census region at different proportions of treated units", cex = 3,
     xpd = TRUE)
dev.off()