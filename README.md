# GOPAL
Scripts used for parameterizing GOPAL (Gibbs Optimized Potential for Atomistic simulation of Liquid water)

The most important scripts developed for carrying out the water model parameterization have been shared. Note the scripts are written in python 2.X and for Some of the scripts like the IAPWS95 implementation (in getexp.py) and ideal gas corrections (in get_corr.py) can be used directly.  

* krsna_t4_rc_final.py : Code used for running optimization 
* getexp.py : code for IAPWS95 implementation  
* get_corr.py : code for implementing all the ideal gas corrections 
* getsimprop.py : code for estimating all thermodynamic properties (except the mixed derivatives) using reduced internal energy matrix uklt(k,l,t) of samples generated using molecular simulations. Uncertainties are not estimated. This code is used during the optimization 
* getsimprop2.py : code for estimating all thermodynamic properties using reduced internal energy matrix uklt(k,l,t). Uncertainties are also estimated. This code is used for generating the property estimates along with uncertainty estimates for final report. 
* gop2_c16_surf_t : Nomads input file 
* genPvec_v2.py : code that generates input vector for a thermodynamic state defined by [ T, P, roh, rhh, rom, q_o, sigma_o, epsilon_o]
* getenev3bnew_pp_a3.py : code that generates energies u[k,l,t] by reading and mapping (if required) trajectories 
* gentrr.py : code used to generate trajectories 





