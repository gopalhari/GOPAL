# GOPAL
Scripts used for parameterizing GOPAL (Gibbs Optimized Potential for Atomistic simulation of Liquid water)

getexp.py : code for IAPWS95 implementation 
getcorr.py : code for implementing all the ideal gas corrections
getsimprop.py : code for estimating all thermodynamic properties (except the mixed derivatives) using reduced internal energy matrix uklt(k,l,t). Uncertainties are not estimated. This code is used during the optimization
getsimprop2.py : code for estimating all thermodynamic properties using reduced internal energy matrix uklt(k,l,t). Uncertainties are also estimated. This code is used during for generating the property estimates along with uncertainty estimates.

