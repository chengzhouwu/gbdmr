gbdmr
=======

An R package for identification of differentially methylated regions (DMRs) from high density chip data. 

Algorithm
--------
1. Order the CpG sites based on the chromosome number
and chromosome coordinates. Calculate the adjacent CpG sites’
correlation in the same chromosome one by one. Furthermore,
get the clustered CpG sites that share similar information based
on the correlation and treat these clusters as candidate DMRs.
2. Generalized beta regression is used to model the CpG
sites in the same cluster; parameter estimation is performed
via maximum likelihood function (ML).
3. Conduct a log likelihood ratio test to the interesting
phenotype after adjusting for the covariables, then summarize
the result and report the DMRs

Installation
------------

```s
library(devtools)
install_github("chengzhouwu/gbdmr")
```
To start using the package just load it as any other package.

```s
library(gbdmr)
```

Usage
-----
For running the algorithm one has to have three objects:

* a matrix with methylation values;
* a phenotype dataset (the IDs of the two data sets need to correspond one-to-one); 
* specific correlation value from 0 to 1.  

All the work is done by one command `gbdmr` and for more information use `?gbdmr` in R.

Example
-------

An example dataset `beta` and `phenotype` is included in the package. It contains 1,000 CpG sites for 506 samples. Finding regions in this data can be accomplished using commands:

```s
data(beta)
data(phenotype)
DMR_summary_list =  gbdmr(beta,phenotype,rho,debug = F)
```
For the new finding CpG sites, mediation analysis can be doone by:

```s
Y <- new_variable$ASTHMA_18
X <- new_variable$BMI_CLUSTER.x
M <- new_beta[rownames(new_beta) == "cg09705784",]
M <- as.numeric(M[1,])

model.0 <- glm(Y~X)
summary(model.0)

model.M <- glm(M ~ X)
summary(model.M)

model.Y <- glm(Y ~ X + M)
summary(model.Y)

# ACME value
library(mediation)
results <- mediate(model.M, model.Y, treat="X", mediator = "M", boot = T, sims=500)
summary(results)
```

Assocaition study can be done:

```s
Y <- new_variable$ASTHMA_18
X1 <- new_variable$POLLEN_18
X2 <- new_beta[rownames(new_beta) == "cg22024479",]
X2 <- as.numeric(X2[1,])
X3 <- new_variable$SEX_0.x

model1 <- glm(Y ~ X1 + X2 + X3 + X1*X2 + X1*X3 + X2*X3)
summary(model1)
model2 <- glm(Y ~ X1 + X2)
summary(model2)
model3 <- glm(Y ~ X1 + X2 + X3 )
summary(model3)
coef <- model2$coefficients[2]
se <- summary(model2)$coefficients[2,2]
OR <- exp(round(coef,2))
OR
round(exp(round(coef,2)+ 1.96*se),2)
round(exp(round(coef,2)-1.96*se),2)
```
Extended application
-------

A [simple example](https://rpubs.com/chengzhouwu/973206)
showing how to apply *gbdmr* to a microarray dataset  and cluster the CpG sites by correlation.




