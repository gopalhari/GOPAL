#!/usr/bin/python 

import numpy
import pdb
# module to generate vector of parameters required by gentrr module to run initial simulations

# inputs to this module
# nlam = number of intermediate states between TIP4P and SPC/E
# T = vector of temperautres at which the simulations have to be run
# P = vector of pressures at which the simulations have to be run 

# Output of this module
# (FF,T,P)
# FF = force field parameters FF = [ROH(nm),RHH(nm),a,Qo,C6,C12]
# T = Temperature T = [T1,T2,T3...]
# P = Pressure P = [P1,P2,P3,..]


def genvec(fset,TP):

	finvec = fset
	[nTP,bt] = numpy.shape(TP) 
	nFF = (numpy.size(fset))/6
	ct = 0 
	#pdb.set_trace()
	outPvec = numpy.zeros([nFF*nTP,6+2], numpy.float64)
	for k in range(nFF):
		for i in range(nTP):
			if nFF ==1:
				outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],TP[i,0],TP[i,1]]
			else:
				outPvec[ct,:] = [finvec[k,0],finvec[k,1],finvec[k,2],finvec[k,3],finvec[k,4],finvec[k,5],TP[i,0],TP[i,1]]
			ct = ct+1
	return(outPvec)
				

def genvecTP(fset,T,P):
	finvec = fset
	pdb.set_trace

	nT = len(T)
	nP = len(P)
	nFF = 1
	ct = 0 
	outPvec = numpy.zeros([nFF*nT*nP,6+2], numpy.float64)
	if nP == 1:
		j = 0 
		if nT == 1:
			outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[0],P[j]]
			ct = ct+1
		else: 	
			for k in range(nT):
				outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[k],P[j]]
				ct = ct+1

	else:	
		for j in range(nP):
			if nT == 1:
				outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[0],P[j]]
				ct = ct+1
			else: 	
				for k in range(nT):
					outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[k],P[j]]
					ct = ct+1
	 	
	#if nT == 1:
	#	j = 0 
	#	if nP == 1:
	#		outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[j],P[0]]
	#		ct = ct+1
	#	else: 	
	#		for k in range(nP):
	#			outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[j],P[k]]
	#			ct = ct+1

	#else:	
	#	for j in range(nT):
	#		if nP == 1:
	#			outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[j],P[0]]
	#			ct = ct+1
	#		else: 	
	#			for k in range(nP):
	#				outPvec[ct,:] = [finvec[0],finvec[1],finvec[2],finvec[3],finvec[4],finvec[5],T[j],P[k]]
	#				ct = ct+1
	return(outPvec)
				
