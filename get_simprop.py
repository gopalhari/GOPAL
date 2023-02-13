#!/usr/bin/python


import numpy
import pdb
import os
import math
import getexp
##import old_mbar as pymbar
import old_mbar

	

def getres(scriptdir,jobdir,case,rjobs,resultsjn,nff,at,K,Nk,NTP,NTPUS,NumTP,verbose,pvar,ref_st,x):
	#pdb.set_trace()
	kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)
	m1=15.9994 
	m2=1.008 
	m3=1.008 
	M = m1+m2+m3
	delp = 1E-3 # bar
	dbeta = 1E-4

	def builduklt(k,l,m,Uij,Vij,Jij,u_klt,Pveccomb):
		Ti = Pveccomb[l+m,6] 
		Pi = Pveccomb[l+m,7]
		#print(Ti)
		betai = 1.0/(kB*Ti) 
		[u_klt[k,l+m,:], \
		u_klt[k,l+m+K,:], \
		u_klt[k,l+m+2*K,:], \
		u_klt[k,l+m+3*K,:], \
		u_klt[k,l+m+4*K,:]]= \
		[betai*(Uij+Pi*Vij)-Jij, \
		betai*(Uij+(Pi+delp)*Vij)-Jij, \
		betai*(Uij+(Pi-delp)*Vij)-Jij, \
		(betai+dbeta)*(Uij+Pi*Vij)-Jij, \
		(betai-dbeta)*(Uij+Pi*Vij)-Jij]	

		Pveccomb[l+m,6] = 1.0/(betai*kB)
		Pveccomb[l+m,7] = Pi
		Pveccomb[l+m+K,6] = 1.0/(betai*kB)
		Pveccomb[l+m+K,7] = Pi+delp
		Pveccomb[l+m+2*K,6] = 1.0/(betai*kB)
		Pveccomb[l+m+2*K,7] = Pi-delp
		Pveccomb[l+m+3*K,6] = 1.0/(kB*(betai+dbeta))
		Pveccomb[l+m+3*K,7] = Pi
		Pveccomb[l+m+4*K,6] = 1.0/(kB*(betai-dbeta))
		Pveccomb[l+m+4*K,7] = Pi
		return()
	ukltfile = jobdir+case+'_uklt_ext_cp_.npy'
	Pvecfile = jobdir+case+'_pveccomb_ext_cp_.npy'
	fkfile = jobdir+case+'_fk_ext_cp_.npy'


	if not (os.path.exists(ukltfile)):
		u_klt = numpy.zeros([at,K+4*K,Nk], numpy.float64)# +9K for +-dt,+-dp,+dp-dt,+dp+dt,-dp+dt,-dp-dt at each T,P point
		Pveccomb = numpy.load(Pvecfile)
	else:
		u_klt = numpy.load(ukltfile)
		Pveccomb = numpy.load(Pvecfile)

	for ki in range(0,nff):
		for k in range((ki)*NTP,(ki+1)*NTP):
			if not (os.path.exists(ukltfile)):
				lic=0  # lic count
				for l in range(0,nff*NTP,NTP):
				#pdb.set_trace()
				#Running through first NTP states and building uklt using E1
					for m in range(NTP):
						#print(k,l+m)
						if m == 0:
							#[Uij,Vij,Jij]=jobsn[rjobs]()
							[Uij,Vij,Jij]=resultsjn[rjobs]
							rjobs = rjobs+1
							#print(rjobs)
							builduklt(k,l,m,Uij,Vij,Jij,u_klt,Pveccomb)	
							# runnin through 2NTP to 2NTP+NTPUS-1 states
							for m2 in range(nff*NTP+(lic)*NTPUS,nff*NTP+(lic+1)*NTPUS):
								builduklt(k,l,m2-l,Uij,Vij,Jij,u_klt,Pveccomb)	
						else:	# rest states use the peij, vij and Jij from above
							builduklt(k,l,m,Uij,Vij,Jij,u_klt,Pveccomb)	
					lic = lic+1
			l=nff*(NTP+NTPUS)
			#Running through NTP+NTPUS states for f4 and building uklt using E1M4
			for m in range(NTP+NTPUS):
				if m == 0:
					#[Uij,Vij,Jij]=jobsn[rjobs]()
					[Uij,Vij,Jij]=resultsjn[rjobs]
					rjobs = rjobs+1
					builduklt(k,l,m,Uij,Vij,Jij,u_klt,Pveccomb)	
				else:	# rest states use the peij, vij and Jij from above
					builduklt(k,l,m,Uij,Vij,Jij,u_klt,Pveccomb)	


	if not (os.path.exists(ukltfile)):
		numpy.save(ukltfile,u_klt)
		numpy.save(Pvecfile,Pveccomb)


	N_k = numpy.zeros(5*K, int) # N_k[k] is the number of uncorrelated samples from state k
	for k in range(at):
		N_k[k] = Nk


	if not (os.path.exists(fkfile)):

		mbar = old_mbar.MBAR(u_klt, N_k, method = 'adaptive',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initialize = 'BAR')
		#mbar = pymbarsee.MBAR(u_klt_new, N_k, method = 'self-consistent-iteration',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initialize = 'BAR')
		if pvar == 0: 
			fij = mbar.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			#fij_id = mbarid.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			Deltaf_ij = fij[0]
			#Deltafid_ij = fij_id[0]
		else:
			(Deltaf_ij, dDeltaf_ij,theta) = mbar.getFreeEnergyDifferences()
			#(Deltafid_ij, dDeltafid_ij,thetaid) = mbarid.getFreeEnergyDifferences()

		numpy.save(fkfile,Deltaf_ij)
	else:
		fkg0 = numpy.load(fkfile)
		fkg = fkg0[0]
		mbar = old_mbar.MBAR(u_klt, N_k, maximum_iterations=50, method = 'adaptive',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initial_f_k=fkg)
		#mbar = pymbarsee.MBAR(u_klt_new, N_k, method = 'self-consistent-iteration',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initial_f_k=fkg)

		if pvar == 0: 
			fij = mbar.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			#fij_id = mbarid.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			Deltaf_ij = fij[0]
			#Deltafid_ij = fij_id[0]
		else:
			(Deltaf_ij, dDeltaf_ij,theta) = mbar.getFreeEnergyDifferences()
	#print "MBARclac time *******",t2-t1

	#pdb.set_trace()	
	N = 521.0
	Deltaf_ij = (Deltaf_ij)/N 
	UncDeltaf_ij = numpy.copy(Deltaf_ij) # Unchanged Delta_fij needed for pressure derivatives 
	if pvar == 1:
		dDeltaf_ij = dDeltaf_ij/N
		theta = theta/(N**2)
	else:
		pass

	for i in range(5*K):
		for j in range(5*K):
			Ti = Pveccomb[i,6]
			Pi = Pveccomb[i,7]
			Tj = Pveccomb[j,6]
			Pj = Pveccomb[j,7]
			#print "Deltaf_ij , Ti, Pi, Tj, Pj,ratio",Deltaf_ij[i,j],Ti,Pi,Tj,Pj, numpy.log(Pj*Ti/(Pi*Tj))
			Deltaf_ij[i,j] = Deltaf_ij[i,j] - (1.0+1.0/N)*numpy.log(Pj*Ti/(Pi*Tj)) # f_res(T,P,U0) = f(T,P,U0)_liquid - f(T,P,U0)_IG


	def getmuliq(st):
		rohlami = Pveccomb[st,0]
		rhhlami = Pveccomb[st,1]
		romlami = Pveccomb[st,2]
		sinthetaby2i = rhhlami/2/rohlami;
		anglami = math.asin(sinthetaby2i)*2;
		qolami = Pveccomb[st,3]
		qhlami = -qolami/2
		remu = 2.0*qhlami*(rohlami*numpy.cos(anglami/2)-romlami)/0.020822678	
		return(remu)

	def getz(Tinp,Pinp):
		z = Pinp*1E5*18.0154E-6/(getexp.getrho(Tinp,Pinp)*kB*Tinp)
		return(z)

	def getd2lnzBdbeta2(Tinp,Pinp):
		beta = 1.0/(kB*Tinp)
		beta1 = beta+dbeta
		beta2 = beta-dbeta
		T1 = 1.0/(kB*beta1)
		T2 = 1.0/(kB*beta2)
		delbeta = beta1-beta2
		d2lnzbdbeta2 = (-kB*(beta**2))*(( numpy.log(getz(T1,Pinp)) -2*numpy.log(getz(Tinp,Pinp)) + numpy.log(getz(T2,Pinp)) )/(0.5*delbeta)**2 )*1000
		return(d2lnzbdbeta2)
		
	#pdb.set_trace()
	# estimation of properties at reference state i.e. the first state for guessed model
	# this calculation has to be done only once 		
	a = ref_st
	#pa2kJpmol = (1e5)*(1e-09)**3 * 6.02214 * (1e23) / 1000 # convert Panm3 to kJ/mol
	Ta = Pveccomb[a,6]
	Pa = Pveccomb[a,7]
	betaa = 1.0/(kB*Ta)

	beta1a = betaa+dbeta
	beta2a = betaa-dbeta
	delbetaa = beta1a-beta2a
	Pra = Pa*1E5
	mvola = ((UncDeltaf_ij[2*K+a,K+a]/betaa/(2*delp*1E5))*1000)/((m1+m2+m3)/1000) # m3/kg
	mvolresa= (mvola*(m1+m2+m3)/1000 - 1000*kB*Ta/(Pa*1E5)) # residual volume m3/mol at standard state
	hresa = (Deltaf_ij[4*K+a,3*K+a])/(delbetaa) 
	uresa = hresa - (Pa*1E5*mvolresa)*1E-3 # kJ/mol
	#hresa = hresa - uresa

	expgresa = getexp.getgresTP(Ta,Pa) - kB*Ta*numpy.log(getz(Ta,Pa)) # Gresa kJ/mol
	expsresa = getexp.getsresTP(Ta,Pa) + kB*numpy.log(getz(Ta,Pa)) # Sresa kJ/mol/K + kB*log(z) to convert from P(T,rho) id.g.st to V(T,P) EOS
	#expuresa = getexp.gethresTP(Ta,Pa) - Pa*1E5/1E3*getexp.getvresTP(Ta,Pa)   # Uresa kJ/mol/K
	expuresa = getexp.geturesTP(Ta,Pa)   # Uresa kJ/mol/K

	exp_prop_file = jobdir+case+'_exp_prop_ext_cp.npy'

	def objfunc(st,nTP):
		starting_k = st
		sum = 0
		rescpres = 0
		resalpha = 0
		reskappa = 0
		res3 = 0
		res4 = 0
		res5 = 0
		resg = 0
		resvres = 0 
		reshres = 0 
		resures = 0 
		ressres = 0 
		resalpha = 0
		reskappa = 0
		resv = 0 
		sumdffpg = 0 
		sumdmvolres = 0

		resgdP = 0 
		resgdB = 0 
		varres0 = 0
		varres1 = 0
		varres2 = 0
		varres3 = 0
		varres4 = 0
		varres5 = 0
		varresg = 0 
		#BigT = numpy.zeros([nTP], numpy.float64)	
		#BigP = numpy.zeros([nTP], numpy.float64)	
		#Biggsim = numpy.zeros([nTP], numpy.float64)	

		if not (os.path.exists(exp_prop_file)):
			exp_prop_array =  numpy.zeros([nTP,2+5+2], numpy.float64)
			for j in range(nTP): # goes over all T and P
				l = starting_k + j #+ K*btp 
				T = Pveccomb[l,6]
				P = Pveccomb[l,7]
				exprho = getexp.getrho(T,P)
				expmvol = ( 1/exprho*(18.015268)/1000 ) # experimental resdiual volume m3/mol
				expmvolres = ( 1/exprho*(18.015268)/1000 - 1000*kB*T/(P*1E5) ) # experimental resdiual volume m3/mol

				# estimate ures
				expuresl = getexp.geturesTP(T,P) 
				# estimate sres
				expsresl = getexp.getsresTP(T,P) + kB*numpy.log(getz(T,P)) # Sres,i  in kJ/mol/K + kB*log(z) to convert from P(T,rho) id.g.st to V(T,P) EOS
				expsres_la = (expsresa - expsresl)*1000 # Sres,a - Sres,i  From IAPWS 95 J/mol/K
				# estimate gres
				expgresl = getexp.getgresTP(T,P) - kB*T*numpy.log(getz(T,P))    # Gresi(T,V) -lnZ in kJ/mol
				expg = (expgresa/(kB*Ta) - expgresl/(kB*T)) # Ga-Gl in reduced units
				expp_cpres = getexp.getcpresTP(T,P) 
				expp_kappa = getexp.getkappaTP(T,P) 

				exp_prop_array[j] =numpy.array([T,P,exprho,expmvolres,expuresl,expsres_la,expg,expp_cpres,expp_kappa], numpy.float64)  

			numpy.save(exp_prop_file,exp_prop_array)
		else:
			exp_prop_array = numpy.load(exp_prop_file)

		for j in range(nTP): # goes over all T and P
			
			l = starting_k + j #+ K*btp 
			T = Pveccomb[l,6]
			P = Pveccomb[l,7]
			# Retrieving experimental properties from array
			expmvolres = exp_prop_array[j,3]
			expuresl   = exp_prop_array[j,4]
			expsres_la = exp_prop_array[j,5]
			expg       = exp_prop_array[j,6]	
			exp_cpres  = exp_prop_array[j,7]	
			exp_kappa  = exp_prop_array[j,8]	
		
			# calculating vres
			beta = 1/(kB*T) 
			Pr = P*1E5
			mvol = ((UncDeltaf_ij[2*K+l,K+l]/beta/(2*delp*1E5))*1000) # m3/mol
			if pvar ==0:
				mvolres= (mvol - 1000*kB*T/(P*1E5)) # residual volume m3/mol
			else:
				mvolres= (mvol - 1000*kB*T/(P*1E5)) # residual volume m3/mol
				dmvolres = ((dDeltaf_ij[2*K+l,K+l]/beta/(2*delp*1E5))*1000) # error in estimates residual volume m3/mol
				


			resvres = resvres + ((mvolres - expmvolres)*1E9)**2 # Residual in Vres mm3/mol
			#resv = resv + ((mvol - expmvol)*1E9)**2 # Residual in V in (mm3/mol) 
			
			if verbose:
				print l,T,P,mvolres,resvres	
			if pvar == 1:
				#sumdmvolres = sumdmvolres+ (dmvolres*(2/expmvolres*(mvolres/expmvolres -1)))**2 used in previous iteration
				sumdmvolres = sumdmvolres+ (dmvolres*1E9*(1E9*2*(mvolres-expmvolres)))**2 # (mm3/mol)^2
			else:
				pass	
			
			# calculating hres refereced using uresa 
			beta1 = beta+dbeta
			beta2 = beta-dbeta
			delbeta = beta1-beta2
			if pvar ==0:
				hres = (Deltaf_ij[4*K+l,3*K+l])/(delbeta) 
				#return(hres)
			else:
				hres = (Deltaf_ij[4*K+l,3*K+l])/(delbeta) 
				dhres = dDeltaf_ij[4*K+l,3*K+l]/(delbeta)
				#return(hres,dhres)
		
			uresl = hres - (P*1E5*mvolres)*1E-3  # kJ/mol
			resures = resures + (uresl-expuresl)**2  # Residual in Ures

			sres_la = (kB*(-Deltaf_ij[l,a] + (betaa*hresa - beta*hres)))*1000 # Sres,a - Sres,i  From force field in J/mol/K
			if expsres_la == 0 :
				pass
			else:
				#ressres = ressres + (sres_la /expsres_la-1)**2 # Residual in Sres used in previous iteration
				ressres = ressres + (sres_la -expsres_la)**2 # Residual in Sres

			if pvar == 0:
				ffpg = Deltaf_ij[l,a]# simulated residual free energy
			else:
				ffpg = Deltaf_ij[l,a]
				dffpg = dDeltaf_ij[l,a]
			

			if expg == 0 :
				pass
			else:
				#resg = resg + (Deltaf_ij[l,a]/expg-1)**2
				resg = resg + (Deltaf_ij[l,a]-expg)**2
				if pvar ==0:
					pass
				else:
					sumdffpg = sumdffpg + (2*(Deltaf_ij[l,a]-expg)*dDeltaf_ij[l,a])**2 
		
			cpres = (-kB*(beta**2)*((Deltaf_ij[0+l,3*K+l]+Deltaf_ij[0+l,4*K+l])/(0.5*delbeta)**2))*1000  # in J/mol/K
			rescpres = rescpres + (cpres - exp_cpres)**2 # in (J/mol/K)^2
			kappa = ((UncDeltaf_ij[2*K+l,0+l] + UncDeltaf_ij[K+l,0+l])/(delp*1E5*1E-5)**2)/( UncDeltaf_ij[2*K+l,K+l]/(2*delp*1E5*1E-5) ) # in 1/bar
			reskappa = reskappa + ((kappa - exp_kappa)*1E5)**2 # in 1/pa^2

				
			#out6 = "%8.2e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n"%(T,P,x[0],x[1],x[2],x[3],x[4],x[5],(mvolres-expmvolres)*1E9,mvolres,expmvolres,uresl,expuresl,sres_la,expsres_la,Deltaf_ij[l,a],expg,cpres,exp_cpres,kappa,exp_kappa)
			#filename6 = scriptdir+"/loggopal16"+"_krsna_test2_debug_gopal2"
			#f6 = open(filename6, 'a')
			#f6.write(out6)
			#f6.close()


		resgres = resg /nTP
		resvres = resvres /nTP
		resures = resures /nTP
		ressres = ressres /nTP
		rescpres = rescpres /nTP#res0 =    res0/nTP
		reskappa = reskappa /nTP#res0 =    res0/nTP
		
		if pvar == 1:
			stdffpg = (sumdffpg)**0.5
			stdmvolres = (sumdmvolres)**0.5
			vardglim  = stdffpg - resgres
			varvollim = stdmvolres - resvres
			F = numpy.copy(fstiter) 
			mu = getmuliq(st) 
			out5 = "%8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %4.2f\n"%(F[0],F[1],F[2],F[3],F[4],F[5],resures,ressres,resgres,resvres,stdffpg,stdmvolres,t2-t0)
			filename = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/"+str(int(runnum))+"of"+str(int(2))+"/loggopal_ssp_gop22_sp2_a1_vres"
			f = open(filename, 'a')
			f.write(out5)
			f.close()
			return(resvres,resures,ressres,resgres,mu,vardglim,varvollim)
		else:
			#F = numpy.copy(fstiter) 
			#mu = getmuliq(st) 
			#t2 = time.time()
			#out5 = "%8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %4.2f\n"%(F[0],F[1],F[2],F[3],F[4],F[5],resures,ressres,resgres,resvres,rescpres,reskappa,t2-t0)
			#filename = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/"+str(int(runnum))+"of"+str(int(2))+"/loggopal_ext_cp"+runcase
			#f = open(filename, 'a')
			#f.write(out5)
			#f.close()
			return(resvres,resures,ressres,resgres,rescpres,reskappa)
		

	if pvar ==1:
		[resVres,resUres,resSres,resGres,mu,vardglim,varvollim] = objfunc(nff*(NTP+NTPUS),NTP+NTPUS) 

		#fval = resVres
		fval = resVres
		# International skeleton tables IST85 suggests allowable tolerances in properties of water v  = 0.01% h = 0.1%
		# MSE in v = (0.0001)**2 MSE in h = (0.001)**2 = 1E-6 no tolerances given for s and g so use the tolerance of h 
		c1 = resUres - 1E-1
		c2 = resSres - 1E-1
		c3 = resGres - 1E-1
		cmu1 = 1.9 - mu
		cmu2 = mu - 3.5
		print  "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(fval,c1,c2,c3,cmu1,cmu2,vardglim,varvollim)
	else:
		#[resVres,resUres,resSres,resGres,mu] = objfunc(nff*(NTP+NTPUS),NTP+NTPUS) 
		[resVres,resUres,resSres,resGres,resCpres,resKappa] = objfunc(nff*(NTP+NTPUS),NTP+NTPUS) 

		#sfval = sfval+resGres
		# International skeleton tables IST85 suggests allowable tolerances in properties of water v  = 0.01% h = 0.1%
		# MSE in v = (0.0001)**2 MSE in h = (0.001)**2 = 1E-6 no tolerances given for s and g so use the tolerance of h 
		#sc1 = sc1 + resUres/1000.0 #- 1E-5
		#sc2 = sc2 + resSres #- 1E-5
		#sc3 = sc3 + resVres #- 83#1E-1
		#sc4 = sc4 + resCpres #- 5#1E-1
		#sc5 = sc5 + resKappa #- 1E-5
		return(resGres,resUres,resSres,resVres,resCpres,resKappa)	

		#cmu1 = 1.9 - mu
		#cmu2 = mu - 3.5
		#print  "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(fval,c1,c2,c3,c4,c5)
		#return(fval,c1,c2,c3,c4,c5)
	#outf.close()


	##fval1,c11,c21,c31,c41,c51 = myfunc(x,1)
	##fval2,c12,c22,c32,c42,c52 = myfunc(x,5000)
	##
	##fval3 = fval1 + fval2
	##c13 = c11 + c12 -1E-5
	##c23 = c21 + c22 -1E-5
	##c33 = c31 + c32 - 83
	##c43 = c41 + c42 - 5
	##c53 = c51 + c52 -1E-5

#print "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(sfval,sc1,sc2,sc3,sc4,sc5)


