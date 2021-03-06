# Testing for Interactions in Predator-Prey Time Series
The goal of this code is to provide significance testing for interactions between time series by using surrogate data methods. The wavelet coherence of the time series is estimated and then its significance is tested.
Three different methods to generate surrogate data are presented. Depending on the randomization rule each method tests for a different null hypothesis.  
H_0 scrambling method: The observed data is a realisation of independently and identically
distributed random variables (white noise).  
H_0 AAFT method: The observed data is a realisation of an underlying stationary linear Gaussian process measured by an invertible time-independent instantaneous measurement function.  
H_0 IAAFT_td method: The observed data is a realisation of an underlying linear Gaussian process with possibly a deterministic trend and slowly changing parameters.  
To provide testing with a high specificity, the method has to be chosen according to the observed process.
The predator-prey data investigated here was published by:  
Blasius, B., Rudolf, L., Weithoff, G., Gaedke, U., and Fussmann, G. F. Long-term cyclic persistence in an experimental predator–prey system. Nature (2020), 577(7789):226–230.
