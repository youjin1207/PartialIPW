# Partially pooled propensity score model for average treatment effect estimation with multilevel data

Youjin Lee, Trang Q. Nguyen, and Elizabeth Stuart

## Overview

In an observational study, subjects are clustered together, sharing environment and resources. These kinds of the contextual environment are often unmeasured while they are likely to confound and modify the causal relationship which we are interested in. In our paper, we discuss the propensity score weighting method under multilevel setting where units in the cluster share unmeasured confounders or effect modifiers. Furthermore, we apply our proposed method to evaluate the effectiveness of a center-based pre-school program on the child's achievement at kindergarten. Here we provide `R` codes to reproduce our simulation studies and replicate our data analysis on Early Childhood Longitudinal Study, Kindergarten (ECLS-K) data. 

## Data Availability

The ECLS-K dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. A list of variables used for our analysis is as follows:
```
VarsToKeep = c("CHILDID", "S1_ID", "S2_ID", "P1PRIMPK",  
  "C1R4MSCL", "GENDER", "WKBLACK", "WKHISP", 
  "C1HEIGHT", "C1WEIGHT", "P1AGEENT", "C2SCREEN",
  "P1HPARNT", "P1HFAMIL", "C1CMOTOR", "R1_KAGE", 
  "WKSESL", "P1HIG_1", "WKINCOME", "CREGION",
  "KURBAN_R")
```

`ECLS_K_application.R` and `ECLS_K_exploratory.R` can be used to replicate our data analyses and generate the tables and figures presented in the article.

`Data/bootresult.RData` is provided to construct empirical confidence intervals of the estimates, and how to generate the bootstrap samples is described in `ECLS_K_application.R`.

## Code for Reproducibility

* `sim_groupings.R`  
 This `R` file can be used to compare four different grouping methods: grouping by minimizing the within-group distance of (1) the proportion of the treated units, (2) the observed covariates, (3) the observed covariates and the proportion, and (4) grouping randomly. This code also includes the derivation of approximated bias $\Lambda$ due to unmeasured characteristics. We can control the amount of confounding and effect modification from the unmeasured covariates by adjusting these three parameters:
``` 
alpha4 = beta4 = kappa4 = c(-2, 0, 2)
 ## parameters to control ##
q = l = m = 1 ## q,l,m = 1,2,3
```
At the final step, the code prints out true ATE and (ATE estimate, estimated standard error, $\Lambda$) under four different grouping methods.

* `sim_estimators.R`

This `R` file is for deriving different causal estimators with two choices: (i) propensity score using fully pooled observations and partially pooled observations and (ii) marginal IPW, group-weighted IPW, and cluster-weighted IPW. Here we partially pool the cluster groups having a similar proportion of the treated in the cluster. 

We can also control the amount of confounding and effect modification from the unmeasured covariates by adjusting these three parameters `q,l,m`.
This code prints out true ATE and (ATE estimate, estimated standard error) under $2 \times 3 = 6$ different causal estimation methods. 

* `sim_collider.R` generates the case where conditioning on the observed covariates creates a collider bias. 

* `sim_bias_amplification.R` generates the simulation scenario with the observed covariates that are strongly associated with the treatment but weakly associated with the outcome. 

* `sim_randomslope.R` generates the simulation data with random slopes to create center-specific effect heterogeneity of the cluster-level unmeasured confounders. 



