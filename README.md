# AN-EXTROPY-BASED-DIVERGENCE-MEASURE
NONPARAMETRIC INFERENCE FOR AN EXTROPY-BASED DIVERGENCE MEASURE

# Abstract. 
Survival extropy, which quantifies the uncertainty associated with the remaining lifetime distribution, provides an information-theoretic perspective on survival behavior. We consider a divergence measure based on survial extropy and derive its nonparametric estimators based on U statistics, empirical distribution functions, and kernel density. Further, we construct confidence intervals for the divergence measure using the jackknife empirical likelihood (JEL) method and the normal approximation method with a jackknife pseudo value based variance estimator. A comprehensive Monte Carlo simulation study is conducted to compare the performance of the measure with existing divergence measures. In addition, we evaluate the finite sample performance of various estimators for the proposed measure. The findings highlight the effectiveness of the divergence measure and its estimators in practical applications. Finally, we show how the proposed divergence measure is used to detect small differences between images in image datasets, which is common in biomedical studies.

Keywords: Extropy; Jackknife empirical likelihood; Measure of divergence; U statistics.

# Nonparametric Inference for an Extropy-Based Divergence Measure

This repository contains the R code used in the paper:


## Repository:

```text
Section 5:
   Section 5.1:
       Sim_RelativeMSE_comparison_Exponential_dist.R
       Sim_RelativeMSE_comparison_Weibull_dist.R
   Section 5.2:
       Sim_MSE_Estimators_Exponential_dist.R
       Sim_MSE_Estimators_Weibull_dist.R
   Section 5.3:
       Sim_Confidence_Intervals_Exponential_dist.R
       Sim_Confidence_Intervals_Weibull_dist.R
   Section 5.4:
       Censored Simulations/
           Point_estimation_right_censoring.R
           Confidence_Intervals_right_censoring.R

Section 6:
   Real Data Analysis/
       Censored_Real_Data_Analysis.R
       Image_based_Real_Data_Analysis.R
       Images/

