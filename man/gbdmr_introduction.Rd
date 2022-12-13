\name{gbdmr}
\alias{gbdmr}
\title{gbdmr: Identifying differentially methylated CpG
regions in the human genome via generalized beta
regression method}
\Description{Identifying differentially methylated CpG
regions efficiently with generalized beta regression}
\usage{
gbdmr(beta,phenotype,rho,debug = F)
}
\arguments{

\item{beta}{A matrix where columns are samples and rows correspond to the sites}
\item{phenotype}{A matrix where columns are variables and rows correspond to the samples}
\item{rho}{A correlation value ranges from 0 to 1}
}

\details{
gbdmr is used to identify the differential methylated CpG regions for the microarray data. The method can be smoothly conducted in three steps. First, order the CpG sites based on the chromosome number and chromosome coordinates; calculate the adjacent CpG sites’ correlation in the same chromosome one by one; get the clustered CpG sites that share similar information based on the correlation and treat these clusters as candidate DMRs. Then generalized beta regression is used to model the CpG sites in the same cluster; parameter estimation is performed via the maximum likelihood function (ML). Finally, conduct a log likelihood ratio test to the interesting phenotype after adjusting for the confounders, then summarize the result and report the DMRs. The method can be used to find the risk CpG sites related to any disease or CpG sites that can be potential mediators in the causal inference.
}

\author{ChengzhouWu, Xichen Mou, Hongmei Zhang}
\examples{
data(beta)
data(phenotype)
rho <- seq(0,1,0.1)
i <- 1:length()

# Find regions
gbdmr(beta,phenotype,rho)
# rho should be a specific value eg. 0.5
# the result will return three list: one all test statistics; one positive CpG sites test statistics; one positive CpG names.

# find the trait of the CpG sites:
https://ngdc.cncb.ac.cn/ewas

##mediation analysis



}
