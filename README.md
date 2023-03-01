# VF-NLM
Code for variance factored non local means (VF-NLM)
Based on the non local means (NLM) algorithm.
VF-NLM is a modified version that accounts for spatially varying noise.

For choosing parameters: code for calculating and optimizing Stein's Unbiased Risk Estimate (SURE) of the MSE

## Papers
Non local means algrithm (Buades et al): https://hal.science/hal-00271141/document
Stein's Unbiased Risk Estimate (Stein)



## Functions

* Non local means (NLM):
    *   file: nlm.cpp
    *   call: nlf(x, h, alpha, pad)
    *   description: Apply NLM to noisy signal
    
* Varince factoered non local means (VF-NLM):
    *   file: vfNlm.cpp
    *   call: vf_nlf(x, h, alpha, sig, pad)
    *   description: Apply VF-NLM to noisy signal
    
* VF-NLM with Stein's Unbiased Risk Estimate (SURE):
   * file: vfNlmSure.cpp
   * call: vf_nlf_sure(x, h, alpha, sig, pad)
   * description: Apply VF-NLM to noisy signal and calculate SURE of the MSE for the estimate

