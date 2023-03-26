## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("HTRX")

## ----setup--------------------------------------------------------------------
library(HTRX)

## -----------------------------------------------------------------------------
## load the data
data(example_data_nosnp)
data(example_hap1)
data(example_hap2)

## -----------------------------------------------------------------------------
## 
example_data_nosnp[41:43,1:6]

## -----------------------------------------------------------------------------
head(example_hap1,3)

## -----------------------------------------------------------------------------
## create haplotype data removing haplotypes rarer than 1%
HTRX_matrix_rmrare = make_htrx(hap1=example_hap1[1:2000,1:4],
                               hap2=example_hap2[1:2000,1:4],
                               rareremove=TRUE,rare_threshold=0.01)

## display the created haplotype data
HTRX_matrix_rmrare[1:3,47:50]

## create haplotype data without removing any haplotypes
HTRX_matrix_allhaps = make_htrx(hap1=example_hap1[1:2000,1:4],
                                hap2=example_hap2[1:2000,1:4])

## create haplotype data while at maximum 3 SNPs can interact
HTRX_matrix_3snphaps = make_htrx(hap1=example_hap1[1:2000,1:4],
                                 hap2=example_hap2[1:2000,1:4],max_int=3)

## compare the numbers of haplotypes created by setting different 'mat_int'
cat(ncol(HTRX_matrix_rmrare),
    ncol(HTRX_matrix_allhaps),
    ncol(HTRX_matrix_3snphaps))


## ----warning=FALSE------------------------------------------------------------
## selecting the best haplotype model using "AIC" from all the haplotypes
CV_results_nocovar <- do_cv(data_nosnp=example_data_nosnp[1:2000,1,drop=FALSE],
                            featuredata=HTRX_matrix_rmrare,
                            sim_times=2,featurecap=4,usebinary=1,
                            method="simple",criteria="BIC",gain=FALSE)

cat('The selected features', as.character(CV_results_nocovar[[2]]),
    'explains \n',mean(CV_results_nocovar[[1]])*100,
    '% average out-of-sample variance')

## ----warning=FALSE------------------------------------------------------------
## selecting the best haplotype model using "BIC" from all the haplotypes
## here we include the sex and age as fixed covariates
CV_results_withcovar <- do_cv(data_nosnp=example_data_nosnp[1:2000,1:3],
                              featuredata=HTRX_matrix_rmrare,
                              sim_times=2,featurecap=8,usebinary=1,
                              method="stratified",criteria="AIC",gain=TRUE)

cat('The selected features', as.character(CV_results_withcovar[[2]]),
    'explains \n', mean(CV_results_withcovar[[1]])*100,
    '% extra average out-of-sample variance')

## ----warning=FALSE------------------------------------------------------------
## selecting the best haplotype model using "BIC"
## we include all the 8 SNPs, but specify at most 4 SNPs can interact
## we also include the sex and age as fixed covariates
cumu_CV_results <- do_cumulative_htrx(data_nosnp=example_data_nosnp[1:2000,1:3],
                                      hap1=example_hap1[1:2000,],
                                      hap2=example_hap2[1:2000,],
                                      sim_times=1,featurecap=8,usebinary=1,
                                      method="stratified",criteria="AIC",
                                      gain=TRUE,max_int=4)

cat('The selected features', as.character(cumu_CV_results[[2]]),
    'explains \n',mean(cumu_CV_results[[1]])*100,
    '% average out-of-sample variance')

