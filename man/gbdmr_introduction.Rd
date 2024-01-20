\name{gbdmr}
\alias{gbdmr}
\title{gbdmr: Identifying differentially methylated CpG
regions in the human genome via generalized beta
regression method}
\Description{Identifying differentially methylated CpG
regions efficiently with generalized beta regression}
\usage{
gbdmr(beta,phenotype,rho,ntask,debug = F)
}
\arguments{

\item{beta}{A matrix where columns are samples and rows correspond to the sites}
\item{phenotype}{A matrix where columns are variables and rows correspond to the samples}
\item{rho}{A correlation value ranges from 0 to 1}
\item{ntask}{parallel computing task number, default the maximum number of parallel tasks available on the machine}
}
\value{
A list containing all candidate regions and their summary statistics.
}
\details{
gbdmr is used to identify the differential methylated CpG regions for the microarray data. The method can be smoothly conducted in three steps. First, order the CpG sites based on the chromosome number and chromosome coordinates; calculate the adjacent CpG sites’ correlation in the same chromosome one by one; get the clustered CpG sites that share similar information based on the correlation and treat these clusters as candidate DMRs. Then generalized beta regression is used to model the CpG sites in the same cluster; parameter estimation is performed via the maximum likelihood function (ML). Finally, conduct a log likelihood ratio test to the interesting phenotype after adjusting for the confounders, then summarize the result and report the DMRs. The method can be used to find the risk CpG sites related to any disease or CpG sites that can be potential mediators in the causal inference.
}

\author{Chengzhou Wu, Xichen Mou, Hongmei Zhang}
\examples{
##############################################################################
#DMR finding
data(beta)
data(phenotype)
rho <- seq(0,1,0.1)
i <- 1:length()

# Find regions
gbdmr(beta,phenotype,rho, ntask = NULL)
# rho should be a specific value eg. 0.5
# the result will return three list: one all test statistics; one positive CpG sites test statistics; one positive CpG names.

# find the trait of the CpG sites:
https://ngdc.cncb.ac.cn/ewas

##############################################################################
##mediation analysis
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
##############################################################################
##association study
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

}
