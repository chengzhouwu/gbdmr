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
library(maxLik)
library(gbdmr)
```

Usage
-----
For running the algorithm one has to have three objects:

* a matrix with methylation values;
* a phenotype dataset (the IDs of the two data sets need to correspond one-to-one); 
* specific correlation value from 0 to 1.  

All the work is done by one command `gbdmr` and for more information use `?gbdmr` in R.

Tutorial
-------
A comprehensive usage protocol for gbdmr has been documented and is available on RPubs: 
This detailed guide provides step-by-step instructions for leveraging the functionalities of gbdmr effectively.
The output generated by gbdmr is versatile, allowing users to conduct various downstream analyses. Some of the notable analyses include enrichment analysis, mediation analysis, association studies, and prediction studies. The flexibility of gbdmr output empowers researchers to delve deeper into the biological and statistical aspects of their data, making it a valuable tool for comprehensive data exploration and interpretation.




