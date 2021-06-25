# PredPreyWaveletAnalysis
The goal of this code is to provide reliable surrogate data testing for time series. 
Three different methods to generate surrogate data are presented with different null hypotheses.
H_0 scrambling method: The observed data is a realisation of independently and identically
distributed random variables (white noise).
H_0 AAFT method: The observed data is a realisation of an underlying stationary linear Gaussian process measured by an invertible time-independent instantaneous measurement function.
H_0 IAAFT_td method: The observed data is a realisation of an underlying linear Gaussian process with possibly a deterministic trend and slowly changing parameters.
To provide testing with a high specificity the method has to be chosen according to the observed process
