# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:15:40 2022

@author: admin
"""

import numpy as np

import getexp

import scipy.integrate as integrate


# Fundamental vibrational frequencies ν (cm−1); 

# intra-molecular(liquid)
intral = 29979245800*np.array([3490,3280,1645]) # s-1 # since 1cm-1 = 29979245800 Hz
#inter-molecular(liquid)
interl = 29979245800*np.array([800,500,200,50]) # s-1

#kB = 1.381*6.02214/4.184 ## cal/mol/K

kB = 1.380649E-23 *6.02214E23 # m2 kg s-2 K-1 mol-1 # J/mol/K

h = 6.62607015E-34 *6.02214E23 # Js mol-1

def getE_CM_vib(ni,T):
    
    E = ni*kB*T
    return(E)

def getE_QM_vib(T):
    
    Evib = 0
    
    for newi in intral:
        
        thetaj = h*newi/kB
        Evib = Evib + thetaj/2 + thetaj*np.exp(-thetaj/T)/(1-np.exp(-thetaj/T))
        
    for newi in interl:
        thetaj = h*newi/kB
        Evib = Evib + thetaj/2 + thetaj*np.exp(-thetaj/T)/(1-np.exp(-thetaj/T))

    return(Evib*kB) 
    

#Correction in kJ/mol

def getEcorr(T):
    
    modes = len(interl) #+ len(intral)
    Ecorr =  getE_QM_vib(T)  - getE_CM_vib(modes,T)      
    return(Ecorr/1000)

def getCv_vibQM(T):
    
    Cvvib = 0
    
    for newi in intral:
        
        thetaj = h*newi/kB
        Cvvib = Cvvib + (thetaj/T)**2*np.exp(-thetaj/T)/(1-np.exp(-thetaj/T))**2
        
    for newi in interl:
            
        thetaj = h*newi/kB
        Cvvib = Cvvib + (thetaj/T)**2*np.exp(-thetaj/T)/(1-np.exp(-thetaj/T))**2

    return(Cvvib*kB)   

def getCv_vibCM(T):
    modes = len(interl) #+ len(intral)
    Cvcm = modes*kB
    return(Cvcm)
    

def getCv_corr(T):
    return(getCv_vibQM(T)-getCv_vibCM(T))


def getU_corr(T0,T): # kJ/mol
    res = integrate.quad(getCv_corr,T0,T)
    return(res[0]/1000)


def get_corrected_Cpid(T):
    #cpid = getexp.getcpidTPig(T,P)
    cpidc = 6/2*kB + getCv_corr(T) + kB
    
    return(cpidc)

def get_CpidbT_f(T):
    return( (6/2*kB + getCv_corr(T) + kB)/T)


def get_corrected_Sid(T0,T): # J/mol/K
    res = integrate.quad(get_CpidbT_f,T0,T)
    return(res[0])
    
def get_Hid_corr(T0,T):
    res = integrate.quad(get_corrected_Cpid,T0,T)
    return(res[0]/1000)


T = 278.15
T0 = 274.15
P = 1.01325
print(getCv_corr(T), getEcorr(T)-getEcorr(T0),getU_corr(T0,T) )

print(get_corrected_Cpid(T))
