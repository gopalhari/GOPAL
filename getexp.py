#!/usr/bin/python 

import numpy 
from scipy import optimize
import pdb
#from freesteam import *

R = 0.46151805*1000 # J/kg/K
Molwt = 18.015268 # g/mol
kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)

coid95 = numpy.array([
[1,0*-8.32044648201,0],\
[2,0*6.68321052680,0],\
[3, 3.00632000000,0],\
[4, 0.01243600000,1.28728967],\
[5, 0.97315000000,3.53734222],\
[6, 1.27950000000,7.74073708],\
[7, 0.96956000000,9.24437796],\
[8, 0.24873000000,27.5075105]], numpy.float64)

co1 = numpy.array([
[1,  0, 1 , -0.5   ,	 0.12533547935523E-1],\
[2,  0, 1 ,  0.875 ,	 0.78957634722828E1],\
[3,  0, 1 ,  1 	  ,     -0.87803203303561E1],\
[4,  0, 2 ,  0.5   ,     0.31802509345418],\
[5,  0, 2 ,  0.75  ,    -0.26145533859358],\
[6,  0, 3 ,  0.375 ,	-0.78199751687981E-2],\
[7,  0, 4 ,  1 	,	 0.88089493102134E-2],\
[8,  1, 1 ,  4 	,	-0.66856572307965],\
[9,  1, 1 ,  6 	,	 0.20433810950965],\
[10, 1, 1 ,  12	,	-0.66212605039687E-4],\
[11, 1, 2 ,  1 	,	-0.19232721156002],\
[12, 1, 2 ,  5 	,	-0.25709043003438],\
[13, 1, 3 ,  4 	,	 0.16074868486251],\
[14, 1, 4 ,  2 	,	-0.40092828925807E-1],\
[15, 1, 4 ,  13 ,	 0.39343422603254E-6],\
[16, 1, 5 ,  9 	,	-0.75941377088144E-5],\
[17, 1, 7 ,  3 	,	 0.56250979351888E-3],\
[18, 1, 9 ,  4 	,	-0.15608652257135E-4],\
[19, 1, 10,  11 ,	 0.11537996422951E-8],\
[20, 1, 11,  4 	,	 0.36582165144204E-6],\
[21, 1, 13,  13 ,	-0.13251180074668E-11],\
[22, 1, 15,  1 	,	-0.62639586912454E-9],\
[23, 2, 1 ,  7 	,	-0.10793600908932],\
[24, 2, 2 ,  1 	,	 0.17611491008752E-1],\
[25, 2, 2 ,  9 	,	 0.22132295167546],\
[26, 2, 2 ,  10 ,	-0.40247669763528],\
[27, 2, 3 ,  10 ,	 0.58083399985759],\
[28, 2, 4 ,  3 	,	 0.49969146990806E-2],\
[29, 2, 4 ,  7 	,	-0.31358700712549E-1],\
[30, 2, 4 ,  10 ,	-0.74315929710341],\
[31, 2, 5 ,  10 ,	 0.47807329915480],\
[32, 2, 6 ,  6 	,	 0.20527940895948E-1],\
[33, 2, 6 ,  10 ,	-0.13636435110343],\
[34, 2, 7 ,  10 ,	 0.14180634400617E-1],\
[35, 2, 9 ,  1 	,	 0.83326504880713E-2],\
[36, 2, 9 ,  2 	,	-0.29052336009585E-1],\
[37, 2, 9 ,  3 	,	 0.38615085574206E-1],\
[38, 2, 9 ,  4 	,	-0.20393486513704E-1],\
[39, 2, 9 ,  8 	,	-0.16554050063734E-2],\
[40, 2, 10,  6 	,	 0.19955571979541E-2],\
[41, 2, 10,  9 	,	 0.15870308324157E-3],\
[42, 2, 12,  8 	,	-0.16388568342530E-4],\
[43, 3, 3 ,  16 ,	 0.43613615723811E-1],\
[44, 3, 4 ,  22 ,	 0.34994005463765E-1],\
[45, 3, 4 ,  23 ,	-0.76788197844621E-1],\
[46, 3, 5 ,  23 ,	 0.22446277332006E-1],\
[47, 4, 14,  10 ,	-0.62689710414685E-4],\
[48, 6, 3 ,  50 ,	-0.55711118565645E-9],\
[49, 6, 6 ,  44 ,	-0.19905718354408],\
[50, 6, 6 ,  46	,	 0.31777497330738],\
[51, 6, 6 ,  50 ,	-0.11841182425981]], numpy.float64)

#i ,ci, di, ti, ni, alpha , beta, gamma, epsilon

co2 = numpy.array([
[52, 0,3.0, 0.0, -0.31306260323435E2, 20.0, 150.0, 1.21, 1.0],\
[53, 0,3.0, 1.0,  0.31546140237781E2, 20.0, 150.0, 1.21, 1.0],\
[54, 0,3.0, 4.0, -0.25213154341695E4, 20.0, 250.0, 1.25, 1.0]], numpy.float64)


co3 = numpy.array([
[55, 3.5, 0.85, 0.2, -0.14874640856724, 28, 700, 0.32, 0.3],\
[56, 3.5, 0.95, 0.2,  0.31806110878444, 32, 800, 0.32, 0.3]], numpy.float64)

Tc = 647.096 # K
rhoc = 322.0 # kg/m3 

def getAres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sum1 = 0  
	for i in range(0,7):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum1 = sum1 + ni*rrho**di*tau**ti 

	sum2 = 0  
	for i in range(7,51):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum2 = sum2 + ni*rrho**di*tau**ti*numpy.exp(-rrho**ci)

	sum3 = 0  
	#for i in range(51,54):
	for i in range(0,3):
		ci       = co2[i][1]
		di       = co2[i][2]
		ti       = co2[i][3]
		ni       = co2[i][4]  
		alphai   = co2[i][5]  
		betai    = co2[i][6]  
		gammai   = co2[i][7]  
		epsiloni = co2[i][8]  
		exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
		sum3 = sum3 + ni*rrho**di*tau**ti*numpy.exp(exponent)


	sum4 = 0  
	#for i in range(51,54):
	for i in range(0,2):
		ai    = co3[i][1]
		bi    = co3[i][2]
		Bi    = co3[i][3]
		ni    = co3[i][4]  
		Ci    = co3[i][5]  
		Di    = co3[i][6]  
		Ai    = co3[i][7]  
		betai = co3[i][8]  
		exponent = -Ci*(rrho-1)**2-Di*(tau-1)**2
		psi = numpy.exp(exponent)
		theta = 1-tau + Ai*((rrho-1)**2)**(1/(2*betai))
		delta = theta**2+Bi*((rrho-1)**2)**ai
		sum4 = sum4 + ni*delta**bi*rrho*psi
	
	sum = sum1+sum2+sum3+sum4
	return(sum)


 
def getApres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sum1 = 0  
	for i in range(0,7):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum1 = sum1 + ni*di*rrho**(di-1)*tau**ti 

	sum2 = 0  
	for i in range(7,51):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum2 = sum2 + ni*numpy.exp(-rrho**ci)*(rrho**(di-1)*tau**ti*(di-ci*rrho**ci))

	sum3 = 0  
	#for i in range(51,54):
	for i in range(0,3):
		ci       = co2[i][1]
		di       = co2[i][2]
		ti       = co2[i][3]
		ni       = co2[i][4]  
		alphai   = co2[i][5]  
		betai    = co2[i][6]  
		gammai   = co2[i][7]  
		epsiloni = co2[i][8]  
		exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
		sum3 = sum3 + ni*rrho**di*tau**ti*numpy.exp(exponent)*(di/rrho-2*alphai*(rrho-epsiloni))


	sum4 = 0  
	#for i in range(51,54):
	for i in range(0,2):
		ai    = co3[i][1]
		bi    = co3[i][2]
		Bi    = co3[i][3]
		ni    = co3[i][4]  
		Ci    = co3[i][5]  
		Di    = co3[i][6]  
		Ai    = co3[i][7]  
		betai = co3[i][8]  
		exponent = -Ci*(rrho-1)**2-Di*(tau-1)**2
		psi = numpy.exp(exponent)
		theta = 1-tau + Ai*((rrho-1)**2)**(1/(2*betai))
		delta = theta**2+Bi*((rrho-1)**2)**ai
		dpsiBrho = -2*Ci*(rrho-1)*psi
		ddeltaBrho = (rrho-1)*(Ai*theta*2/betai*((rrho-1)**2)**(1/(2*betai)-1) + 2*Bi*ai*((rrho-1)**2)**(ai-1))
		ddeltabiBrho = bi*delta**(bi-1)*ddeltaBrho
		sum4 = sum4 + ni*(delta**bi*(psi+rrho*dpsiBrho)+(ddeltabiBrho*rrho*psi) )
	
	sum = sum1+sum2+sum3+sum4
	return(sum)

def getAid95(rho,T):
	#st = steam_pT(P,T)
	#rho = st.rho
	#print rho	
	rrho = rho/rhoc 
	tau = Tc/T
	gamma = 0
	n10 = coid95[0][1]
	n20 = coid95[1][1]
	n30 = coid95[2][1]
	for i in range(3,8):
		ni = coid95[i][1]
		gi = coid95[i][2]
		gamma= gamma+ni*numpy.log(1-numpy.exp(-gi*tau))
	gamma = numpy.log(rrho)+n10+n20*tau+n30*numpy.log(tau)+gamma
	return(gamma)

def getApid95(rho,T):
	rrho = rho/rhoc
	gamma = 1.0/rrho	 
	return(gamma)	

def getAppid95(rho,T):
	rrho = rho/rhoc
	gamma = -1.0/rrho**2	 
	return(gamma)	

def getAptid95(rho,T):
	return(0.0)
	
def getAtid95(rho,T):
	#st = steam_pT(P,T)
	#rho = st.rho
	#print rho	
	rrho = rho/rhoc 
	tau = Tc/T
	gamma = 0
	n10 = coid95[0][1]
	n20 = coid95[1][1]
	n30 = coid95[2][1]
	for i in range(3,8):
		ni = coid95[i][1]
		gi = coid95[i][2]
		gamma= gamma+ni*gi*((1-numpy.exp(-gi*tau))**(-1)-1)
	gamma = n20+n30/tau+gamma
	return(gamma)

def getAttid95(rho,T):
	#st = steam_pT(P,T)
	#rho = st.rho
	#print rho	
	rrho = rho/rhoc 
	tau = Tc/T
	gamma = 0
	n10 = coid95[0][1]
	n20 = coid95[1][1]
	n30 = coid95[2][1]
	for i in range(3,8):
		ni = coid95[i][1]
		gi = coid95[i][2]
		gamma= gamma+ni*gi**2*numpy.exp(-gi*tau)*(1-numpy.exp(-gi*tau))**(-2)
	gamma = -n30/tau**2-gamma
	return(gamma)

def getAppres95(rho,T):
        #st = steam_pT(P*1E5,T)
        #rho = st.rho
        tau = Tc/T
        rrho = rho/rhoc
        sum1 = 0
        for i in range(0,7):
                ci = co1[i][1]
                di = co1[i][2]
                ti = co1[i][3]
                ni = co1[i][4]
                #print co1[i][0],ci,di,ti,ni
                sum1 = sum1 + ni*di*(di-1)*rrho**(di-2)*tau**ti

        sum2 = 0
        for i in range(7,51):
                ci = co1[i][1]
                di = co1[i][2]
                ti = co1[i][3]
                ni = co1[i][4]
                #print co1[i][0],ci,di,ti,ni
                sum2 = sum2 + ni*numpy.exp(-rrho**ci)*(rrho**(di-2)*tau**ti*((di-ci*rrho**ci)*(di-1-ci*rrho**ci)-ci**2*rrho**ci))

        sum3 = 0
        #for i in range(51,54):
        for i in range(0,3):
                ci       = co2[i][1]
                di       = co2[i][2]
                ti       = co2[i][3]
                ni       = co2[i][4]
                alphai   = co2[i][5]
                betai    = co2[i][6]
                gammai   = co2[i][7]
                epsiloni = co2[i][8]
                exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
                sum3 = sum3 + ni*tau**ti*numpy.exp(exponent)*(-2*alphai*rrho**di+4*alphai**2*rrho**di*(rrho-epsiloni)**2-4*di*alphai*rrho**(di-1)*(rrho-epsiloni)+di*(di-1)*rrho**(di-2))


        sum4 = 0
        #for i in range(51,54):
        for i in range(0,2):
                ai    = co3[i][1]
                bi    = co3[i][2]
                Bi    = co3[i][3]
                ni    = co3[i][4]
                Ci    = co3[i][5]
                Di    = co3[i][6]
                Ai    = co3[i][7]
                betai = co3[i][8]
                exponent = -Ci*(rrho-1)**2.0-Di*(tau-1)**2
                psi = numpy.exp(exponent)
                theta = 1.0-tau + Ai*((rrho-1)**2)**(1.0/(2*betai))
                delta = theta**2+Bi*((rrho-1)**2)**ai
                dpsiBrho = -2.0*Ci*(rrho-1)*psi
                d2psiBrho2 = (2*Ci*(rrho-1)**2-1)*2*Ci*psi
                ddeltaBrho = (rrho-1.0)*(Ai*theta*2.0/betai*((rrho-1)**2)**(1/(2*betai)-1) + 2.0*Bi*ai*((rrho-1)**2)**(ai-1))
                d2deltaBrho2 = 1.0/(rrho-1)*ddeltaBrho+(rrho-1)**2*(4*Bi*ai*(ai-1)*((rrho-1)**2)**(ai-2)+ 2*Ai**2*(1.0/betai)**2*(((rrho-1)**2)**(1.0/(2*betai)-1))**2+ Ai*theta*4.0/betai*(1.0/(2*betai)-1)*(((rrho-1)**2)**(1.0/(2*betai)-2)))
                ddeltabiBrho = bi*delta**(bi-1)*ddeltaBrho
                d2deltabiBrho2 = bi*(delta**(bi-1)*d2deltaBrho2+(bi-1)*delta**(bi-2)*(ddeltaBrho)**2)
                sum4 = sum4 + ni*(delta**bi*(2*dpsiBrho+rrho*d2psiBrho2) + 2*ddeltabiBrho*(psi+rrho*dpsiBrho)+ d2deltabiBrho2*rrho*psi )

        sum = sum1+sum2+sum3+sum4
        return(sum)

def getAtres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sum1 = 0  
	for i in range(0,7):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum1 = sum1 + ni*ti*rrho**di*tau**(ti-1) 

	sum2 = 0  
	for i in range(7,51):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum2 = sum2 + ni*ti*rrho**di*tau**(ti-1)*numpy.exp(-rrho**ci)

	sum3 = 0  
	#for i in range(51,54):
	for i in range(0,3):
		ci       = co2[i][1]
		di       = co2[i][2]
		ti       = co2[i][3]
		ni       = co2[i][4]  
		alphai   = co2[i][5]  
		betai    = co2[i][6]  
		gammai   = co2[i][7]  
		epsiloni = co2[i][8]  
		exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
		sum3 = sum3 + ni*rrho**di*tau**ti*numpy.exp(exponent)*(ti/tau-2*betai*(tau-gammai))


	sum4 = 0  
	#for i in range(51,54):
	for i in range(0,2):
		ai    = co3[i][1]
		bi    = co3[i][2]
		Bi    = co3[i][3]
		ni    = co3[i][4]  
		Ci    = co3[i][5]  
		Di    = co3[i][6]  
		Ai    = co3[i][7]  
		betai = co3[i][8]  
		exponent = -Ci*(rrho-1)**2-Di*(tau-1)**2
		psi = numpy.exp(exponent)
		theta = 1-tau + Ai*((rrho-1)**2)**(1/(2*betai))
		delta = theta**2+Bi*((rrho-1)**2)**ai
		dpsiBtau = -2*Di*(tau-1)*psi
		ddeltabiBtau = -2*theta*bi*delta**(bi-1)
		sum4 = sum4 + ni*rrho*(ddeltabiBtau*psi+delta**bi*dpsiBtau)
	
	sum = sum1+sum2+sum3+sum4
	return(sum)

def getAttres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sum1 = 0  
	for i in range(0,7):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum1 = sum1 + ni*ti*(ti-1)*rrho**di*tau**(ti-2) 

	sum2 = 0  
	for i in range(7,51):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum2 = sum2 + ni*ti*(ti-1)*rrho**di*tau**(ti-2)*numpy.exp(-rrho**ci)

	sum3 = 0  
	#for i in range(51,54):
	for i in range(0,3):
		ci       = co2[i][1]
		di       = co2[i][2]
		ti       = co2[i][3]
		ni       = co2[i][4]  
		alphai   = co2[i][5]  
		betai    = co2[i][6]  
		gammai   = co2[i][7]  
		epsiloni = co2[i][8]  
		exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
		sum3 = sum3 + ni*rrho**di*tau**ti*numpy.exp(exponent)*((ti/tau-2*betai*(tau-gammai))**2 - ti/tau**2 - 2*betai)


	sum4 = 0  
	#for i in range(51,54):
	for i in range(0,2):
		ai    = co3[i][1]
		bi    = co3[i][2]
		Bi    = co3[i][3]
		ni    = co3[i][4]  
		Ci    = co3[i][5]  
		Di    = co3[i][6]  
		Ai    = co3[i][7]  
		betai = co3[i][8]  
		exponent = -Ci*(rrho-1)**2-Di*(tau-1)**2
		psi = numpy.exp(exponent)
		theta = 1-tau + Ai*((rrho-1)**2)**(1/(2*betai))
		delta = theta**2+Bi*((rrho-1)**2)**ai
		dpsiBtau = -2*Di*(tau-1)*psi
		ddeltabiBtau = -2*theta*bi*delta**(bi-1)
		d2deltabiBtau2 = 2*bi*delta**(bi-1)+4*theta**2*bi*(bi-1)*delta**(bi-2)
		d2psiBtau2 = (2*Di*(tau-1)**2 -1 )*2*Di*psi
		
		sum4 = sum4 + ni*rrho*(d2deltabiBtau2*psi + 2*ddeltabiBtau*dpsiBtau + delta**bi*d2psiBtau2)
	
	sum = sum1+sum2+sum3+sum4
	return(sum)

def getAptres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sum1 = 0  
	for i in range(0,7):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum1 = sum1 + ni*di*ti*rrho**(di-1)*tau**(ti-1) 

	sum2 = 0  
	for i in range(7,51):
		ci = co1[i][1]
		di = co1[i][2]
		ti = co1[i][3]
		ni = co1[i][4]  
		#print co1[i][0],ci,di,ti,ni
		sum2 = sum2 + ni*ti*rrho**(di-1)*tau**(ti-1)*(di- ci*rrho**ci)*numpy.exp(-rrho**ci)

	sum3 = 0  
	#for i in range(51,54):
	for i in range(0,3):
		ci       = co2[i][1]
		di       = co2[i][2]
		ti       = co2[i][3]
		ni       = co2[i][4]  
		alphai   = co2[i][5]  
		betai    = co2[i][6]  
		gammai   = co2[i][7]  
		epsiloni = co2[i][8]  
		exponent = -alphai*(rrho-epsiloni)**2-betai*(tau-gammai)**2
		sum3 = sum3 + ni*rrho**di*tau**ti*numpy.exp(exponent)*(di/rrho-2*alphai*(rrho-epsiloni))*(ti/tau-2*betai*(tau-gammai))


	sum4 = 0  
	#for i in range(51,54):
	for i in range(0,2):
		ai    = co3[i][1]
		bi    = co3[i][2]
		Bi    = co3[i][3]
		ni    = co3[i][4]
		Ci    = co3[i][5]
		Di    = co3[i][6]
		Ai    = co3[i][7]
		betai = co3[i][8]
		exponent = -Ci*(rrho-1)**2-Di*(tau-1)**2
		psi = numpy.exp(exponent)
		theta = 1-tau + Ai*((rrho-1)**2)**(1/(2*betai))
		delta = theta**2+Bi*((rrho-1)**2)**ai
		dpsiBtau = -2*Di*(tau-1)*psi
		dpsiBrho = -2.0*Ci*(rrho-1)*psi
		ddeltabiBtau = -2*theta*bi*delta**(bi-1)
		d2psiBdrhodtau = 4*Ci*Di*(rrho-1)*(tau-1)*psi
		ddeltaBrho = (rrho-1.0)*(Ai*theta*2.0/betai*((rrho-1)**2)**(1/(2*betai)-1) + 2.0*Bi*ai*((rrho-1)**2)**(ai-1))
		d2deltabiBdrhodtau = -Ai*bi*2.0/betai*delta**(bi-1)*(rrho-1)*((rrho-1)**2)**(1/(2*betai)-1) - 2*theta*bi*(bi-1)*delta**(bi-2)*ddeltaBrho
		ddeltabiBrho = bi*delta**(bi-1)*ddeltaBrho
		ddeltabiBtau = -2*theta*bi*delta**(bi-1)
		sum4 = sum4 + ni*(delta**bi*(dpsiBtau + rrho*d2psiBdrhodtau) + rrho*ddeltabiBrho*dpsiBtau + ddeltabiBtau*(psi+rrho*dpsiBrho) + d2deltabiBdrhodtau*rrho*psi)
	
	sum = sum1+sum2+sum3+sum4
	return(sum)


#print "getAptres95"
#print getAptres95(838.025,500)
#print getAptres95(358,647.0)
def getsid95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sid = tau*(getAtid95(rho,T)) - getAid95(rho,T)
	return(sid) 
	

def getsres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	sres = tau*(getAtres95(rho,T)) - getAres95(rho,T)
	return(sres) 

def gets95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	s = tau*(getAtid95(rho,T)+getAtres95(rho,T)) - getAid95(rho,T) - getAres95(rho,T)
	return(s)

def getuid95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	uid = tau*(getAtid95(rho,T))
	return(uid) 

def getures95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	ures = tau*(getAtres95(rho,T))
	return(ures) 

def getu95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	u = tau*(getAtid95(rho,T) + getAtres95(rho,T))
	return(u)	

def gethid95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	hid = 1 + tau*(getAtid95(rho,T))
	return(hid) 

def gethres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	#pdb.set_trace()
	hres = tau*(getAtres95(rho,T)) + rrho*getApres95(rho,T)
	return(hres) 


def geth95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	h = gethid95(rho,T) + gethres95(rho,T)
	return(h) 

def getgid95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	gid = 1 + (getAid95(rho,T))
	return(gid) 

def getgres95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	gres = (getAres95(rho,T)) + rrho*getApres95(rho,T)
	return(gres) 

def getg95(rho,T):
	g = getgid95(rho,T) + getgres95(rho,T) 
	return(g)
	
def getcpid95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	cpid = -tau**2*getAttid95(rho,T) + 1
	return(cpid) 

def getcp95(rho,T):
	tau = Tc/T
	rrho = rho/rhoc
	cp = -tau**2*(getAttid95(rho,T) + getAttres95(rho,T)) + (1+rrho*getApres95(rho,T)-rrho*tau*getAptres95(rho,T))**2/(1 + 2*rrho*getApres95(rho,T) + rrho**2*getAppres95(rho,T))
	return(cp)


#get saturation properties

def maxwell(x,T):
	Psat = x[0]
	rholsat = x[1]
	rhogsat = x[2]
	# solve maxwell criterion for Psat, rholsat and rhogsat Eqns 6.9a-c
	a = numpy.array([Psat/(R*T*rholsat)-(1+rholsat/rhoc*getApres95(rholsat,T)),\
			 Psat/(R*T*rhogsat)-(1+rhogsat/rhoc*getApres95(rhogsat,T)),\
			 Psat/(R*T)*(1.0/rhogsat-1.0/rholsat)-numpy.log(rholsat/rhogsat)-(getAres95(rholsat,T)-getAres95(rhogsat,T))], numpy.float64)
	#print a
	return(a)

def getsatprops(T):
	x0 = numpy.array([100,1000.1,0.001], numpy.float64)
	sol = optimize.fsolve(maxwell,x0,args=(T),xtol=1E-10)
	Psat = sol[0] # Pa
	rholsat = sol[1] # kg/m3
	rhogsat = sol[2] # kg/m3
	hvap = kB*T*(geth95(rhogsat,T) - geth95(rholsat,T)) # kJ/mol
	return([Psat,rholsat,rhogsat,hvap])
	#return(hvap)


# get experimental properties as function of pressure and temperature



# Density kg/m3
def getrho(T,P):
	P = P*1E5
	#rho0 = 1000
	rho0 = 1000
	def func(x):
		return(P/(R*T*x)-(1+x/rhoc*getApres95(x,T))) 
	sol = optimize.fsolve(func,rho0,xtol=1E-10)
	return(sol)

# Residual gibss free energy kJ/mol 
def getgresTP(T,P):
	rho = getrho(T,P)
	tau = Tc/T
	rrho = rho/rhoc
	gres = kB*T*((getAres95(rho,T)) + rrho*getApres95(rho,T))
	return(gres) 

# gibss free energy kJ/mol 
def getgTP(T,P):
	rho = getrho(T,P)
	tau = Tc/T
	rrho = rho/rhoc
	g = kB*T*(getg95(rho,T))
	return(g) 

# Entropy kJ/mol/K
def getsresTP(T,P):
        rho = getrho(T,P)
        s = getsres95(rho,T)
        s = s*kB
        return(s)

def getsTP(T,P):
        rho = getrho(T,P)
        s = gets95(rho,T)
        s = s*kB
        return(s)

# Molar volume m3/mol
def getvresTP(T,P):
        rho = getrho(T,P)
	#pdb.set_trace()
        molvol = 1/rho*(0.0001*Molwt) - 1000*kB*T/(P*1E5)
        return(molvol)

# Residual enthalpy kJ/mol
def gethresTP(T,P):
	rho = getrho(T,P)
	hres = kB*T*gethres95(rho,T)
	return(hres)


# Internal energy kJ/mol 
def geturesTP(T,P):
        rho = getrho(T,P)
        u = getures95(rho,T)
        u = u*kB*T	
        return(u)

def getuTP(T,P):
        rho = getrho(T,P)
        u = getu95(rho,T)
        u = u*kB*T	
        return(u)
 
# enthalpy kJ/mol
def gethTP(T,P):
	rho = getrho(T,P)
	h = kB*T*geth95(rho,T)
	return(h)

#Heat of vaporization kJ/mol
def gethvap(T):
	a = getsatprops(T) 
	hvap = a[3]
	return(hvap)

# heat capacity J/mol/K
def getcpTP(T,P):
	rho = getrho(T,P)
	cp = 1000*kB*(getcp95(rho,T))
	return(cp)

# Residual heat capacity J/mol/K
def getcpresTP(T,P):
	rho = getrho(T,P)
	cpres = 1000*kB*(getcp95(rho,T) - getcpid95(rho,T))
	return(cpres)

def getz(Tinp,Pinp):
	z = Pinp*1E5*18.0154E-6/(getrho(Tinp,Pinp)*kB*Tinp)
	return(z)

dbeta = 1E-4
def getd2lnzBdbeta2(Tinp,Pinp):
	beta = 1.0/(kB*Tinp)
	beta1 = beta + dbeta
	beta2 = beta - dbeta
	T1 = 1.0/(kB*beta1)
	T2 = 1.0/(kB*beta2)
	delbeta = beta1-beta2
	d2lnzbdbeta2 = (-kB*(beta**2))*(( numpy.log(getz(T1,Pinp)) -2*numpy.log(getz(Tinp,Pinp)) + numpy.log(getz(T2,Pinp)) )/(0.5*delbeta)**2 )*1000
	return(d2lnzbdbeta2)

# Ideal heat capacity J/mol/K
def getcpidTV(T,P):
	rho = getrho(T,P) # real water density kg/m3
	#print rho
	cpid = 1000*kB*(getcpid95(rho,T))
	return(cpid)

# Ideal heat capacity J/mol/K
def getcpidTPig(T,P):
	rho = 18.0154E-3*P*1E5/(kB*1000*T)  # ideal gas density kg/m3
	#print rho
	cpid = 1000*kB*(getcpid95(rho,T))
	return(cpid)

Tn = 274.15
Pn = 1.01325

#print getcpidTV(Tn,Pn),getcpidTPig(Tn,Pn),getcpidTV(Tn,Pn)-getcpidTPig(Tn,Pn)

#print getcpresTP(Tn,Pn), getd2lnzBdbeta2(Tn,Pn), getcpresTP(Tn,Pn) - getd2lnzBdbeta2(Tn,Pn)  
# Thermal expansion coefficient  (1/v*dv/dt @ const p)  1/K
def getalphaTP(T,P):
        dt = 0.0001
        Tpdt = T + dt
        Tmdt = T - dt
        deltaT = Tpdt-Tmdt
        vatT= 1/getrho(T,P)
        vpdt= 1/getrho(Tpdt,P)
        vmdt= 1/getrho(Tmdt,P)
        ealpha = 1/vatT*(vpdt-vmdt)/deltaT
        return(ealpha)

# Isothermal compressibility (-1/v*dv/dp @ const t) 1/bar
def getkappaTP(T,P):
        dp = 0.0001
        Ppdp = P+ dp
        Pmdp = P- dp
        deltaP = Ppdp-Pmdp
        vatP=  1/getrho(T,P)
        vpdp=  1/getrho(T,Ppdp)
        vmdp=  1/getrho(T,Pmdp)
        ekappa = -1/vatP*(vpdp-vmdp)/deltaP
        return(ekappa)


#Ta = 274.15
#Pa = 1.01325
#print gethresTP(Ta,Pa) , getvresTP(Ta,Pa), Pa*1E5/1E3*getvresTP(Ta,Pa), gethresTP(Ta,Pa) - Pa*1E5/1E3*getvresTP(Ta,Pa) 
#pdb.set_trace()
#print geturesTP(372.150,1.2101325e+02)
#print gethresTP(372.150,1.2101325e+02)
#print getsresTP(372.150,1.2101325e+02)
#print getgresTP(372.150,1.2101325e+02)
#print 1.2101325e+02*getvresTP(372.150,1.2101325e+02)	
#print getuTP(298.150,1.01325)
#print gethTP(298.150,1.01325)
#print getsTP(298.150,1.01325)
#print getgTP(298.150,1.01325)
#print getuTP(372.150,1.01325)
#print gethTP(372.150,1.01325)
#print getsTP(372.150,1.01325)
#print getgTP(372.150,1.01325)
#print getsatprops(298.15)
#print getrho(298.15,1.10325)
#print gethresTP95(298.15,1.10325)
#print gethvap(298.15)
#print getcpresTP(298.15,1.10325)
#print getalphaTP(298.15,1.10325)
#print getkappaTP(298.15,1.10325)
#print kB*getcp95(getrho(1.01325,298.15),298.15)
#print kB*getcpid95(getrho(1.01325,298.15),298.15)
#for i in range(98):
#	print 274.15+i,getsatprops(274.15+i)

## corrections
## All corections in kJ/mol
M = 18.015268 # g/mol

# second viriral coeff
def expB(T):
        a = getsatprops(T)
        Psat  = a[0]
        rhogsat  = a[3]
        molv = 1.0/rhogsat*(M/1000) # multiplied by M/1000 to convert to m3/mol
        reB = (Psat*molv/R/T-1)*molv # units of B = m3/mol
        return(reB)

# Cni correction for the non ideality of the gas phase
def corrni(T):
        a = getsatprops(T)
        Psat  = a[0]
        dt = 0.0000001
        Tpdt = T + dt
        Tmdt = T - dt
        deltaT = Tpdt-Tmdt
        B = expB(T)
        Bpdt = expB(Tpdt)
        Bmdt = expB(Tmdt)
        reCni = Psat*(B-T*(Bpdt-Bmdt)/(deltaT))*1E-3 # multiply by 1E-3 to convert from J/mol to kJ/mol
        return(reCni)

muliq = 2.321 #input
mugas = 1.854989 # D
alphagas = 1.6633E-40 # Fm2 one Farad(F) = C2/J

Na = 6.0221415000000003E+023 # Avogadro's number
h = 3.9903127176E-10/1000 # molar Planck's constant kJ.s.mol-1 =Na*h
c = 299792458.0*100 # speed of light cm/s
hc = h*c # kJ/mol-cm
nuGintra = numpy.array([3755.7,3656.6,1594.6], numpy.float64) # 1/cm
nuLintra = numpy.array([3490,3280,1645], numpy.float64) # 1/cm
nuLinter = numpy.array([800,500,200,50], numpy.float64) # 1/cm
DtoSI = 1.0/c*1E-21 # C.m



def EvibQM(nui,T):
        sum1=0
        for i in range(len(nui)):
                sum1 = sum1 + hc*nui[i]/2.0 + (hc*nui[i])/(numpy.exp(hc*nui[i]/kB/T)-1)
        return(sum1)

def dEvibQMdT(nui,T):
        sum2=0
        for i in range(len(nui)):
                sum2 = sum2 + kB*((hc*nui[i]/(kB*T))**2 * numpy.exp(hc*nui[i]/(kB*T))/(numpy.exp(hc*nui[i]/(kB*T))-1)**2)
        return(sum2)

def EvibCM(numi,T):
        return(numi*kB*T)

def dEvibCMdT(numi,T):
        return(numi*kB)

def corrvib(T):
        corrHvibintra = EvibQM(nuGintra,T)- EvibQM(nuLintra,T)
        corrHvibinter = EvibQM(nuLinter,T)- EvibCM(len(nuLinter),T)
        recvib = corrHvibintra-corrHvibinter # fixed !! the second term has a negative sign. See JPM postma's thesis 
        return(recvib)

def corrcpvib(T):
        recorrCpvib = dEvibQMdT(nuLintra,T)+ dEvibQMdT(nuLinter,T)-dEvibCMdT(len(nuLinter),T)
        return(recorrCpvib)

def corrpol(muliqin):
        corrHpol = Na/2.0*((mugas-muliqin)*3.335641E-30)**2/(alphagas)/1000
	#print corrHpol
	#pdb.set_trace()
        return(corrHpol)


def allcorr(T,muliqin):
        tcorr = corrvib(T)+corrni(T)+corrpol(muliq)
        return(tcorr)

#print allcorr(298.0,muliq)
#print corrvib(298.0)/4.184
#print corrni(298.0)/4.184
