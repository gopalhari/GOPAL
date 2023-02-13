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
		##[u_klt[k,l+m,:], \
		##u_klt[k,l+m+K,:], \
		##u_klt[k,l+m+2*K,:], \
		##u_klt[k,l+m+3*K,:], \
		##u_klt[k,l+m+4*K,:]]= \
		##[betai*(Uij+Pi*Vij)-Jij, \
		##betai*(Uij+(Pi+delp)*Vij)-Jij, \
		##betai*(Uij+(Pi-delp)*Vij)-Jij, \
		##(betai+dbeta)*(Uij+Pi*Vij)-Jij, \
		##(betai-dbeta)*(Uij+Pi*Vij)-Jij]	

		##Pveccomb[l+m,6] = 1.0/(betai*kB)
		##Pveccomb[l+m,7] = Pi
		##Pveccomb[l+m+K,6] = 1.0/(betai*kB)
		##Pveccomb[l+m+K,7] = Pi+delp
		##Pveccomb[l+m+2*K,6] = 1.0/(betai*kB)
		##Pveccomb[l+m+2*K,7] = Pi-delp
		##Pveccomb[l+m+3*K,6] = 1.0/(kB*(betai+dbeta))
		##Pveccomb[l+m+3*K,7] = Pi
		##Pveccomb[l+m+4*K,6] = 1.0/(kB*(betai-dbeta))
		##Pveccomb[l+m+4*K,7] = Pi

		
		[u_klt[k,l+m,:], \
		u_klt[k,l+m+K,:], \
		u_klt[k,l+m+2*K,:], \
		u_klt[k,l+m+3*K,:], \
		u_klt[k,l+m+4*K,:], \
		u_klt[k,l+m+5*K,:], \
		u_klt[k,l+m+6*K,:], \
		u_klt[k,l+m+7*K,:], \
		u_klt[k,l+m+8*K,:]]= \
		[betai*(Uij+Pi*Vij)-Jij, \
		betai*(Uij+(Pi+delp)*Vij)-Jij, \
		betai*(Uij+(Pi-delp)*Vij)-Jij, \
		(betai+dbeta)*(Uij+Pi*Vij)-Jij, \
		(betai-dbeta)*(Uij+Pi*Vij)-Jij, \
		(betai+dbeta)*(Uij+(Pi+delp)*Vij)-Jij, \
		(betai-dbeta)*(Uij+(Pi+delp)*Vij)-Jij, \
		(betai+dbeta)*(Uij+(Pi-delp)*Vij)-Jij, \
		(betai-dbeta)*(Uij+(Pi-delp)*Vij)-Jij]	

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
		Pveccomb[l+m+5*K,6] = 1.0/(kB*(betai+dbeta))
		Pveccomb[l+m+5*K,7] = Pi+delp
		Pveccomb[l+m+6*K,6] = 1.0/(kB*(betai-dbeta))
		Pveccomb[l+m+6*K,7] = Pi+delp
		Pveccomb[l+m+7*K,6] = 1.0/(kB*(betai+dbeta))
		Pveccomb[l+m+7*K,7] = Pi-delp
		Pveccomb[l+m+8*K,6] = 1.0/(kB*(betai-dbeta))
		Pveccomb[l+m+8*K,7] = Pi-delp

		return()
	ukltfile = jobdir+case+'_uklt_ext_cp_hk.npy'
	Pvecfile = jobdir+case+'_pveccomb_ext_cp_hk.npy'
	fkfile = jobdir+case+'_fk_ext_cp_hk.npy'

	exp_prop_file = jobdir+case+'_exp_prop_ext_cp_hk.npy'

	if not (os.path.exists(ukltfile)):
		u_klt = numpy.zeros([at,K+4*K+4*K,Nk], numpy.float64)# +9K for +-dt,+-dp,+dp-dt,+dp+dt,-dp+dt,-dp-dt at each T,P point
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


	N_k = numpy.zeros(5*K+4*K, int) # N_k[k] is the number of uncorrelated samples from state k
	for k in range(at):
		N_k[k] = Nk


	if not (os.path.exists(fkfile)):

		mbar = old_mbar.MBAR(u_klt, N_k, method = 'adaptive',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initialize = 'BAR')
		#mbar = old_mbar.MBAR(u_klt, N_k, method = 'adaptive',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initialize = 'zeros')
		#mbar = old_mbar.MBAR(u_klt_new, N_k, method = 'self-consistent-iteration',relative_tolerance=1.0e-10,use_optimized=True,verbose=verbose,initialize = 'BAR')
		if pvar == 0: 
			fij = mbar.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			#fij_id = mbarid.getFreeEnergyDifferences(compute_uncertainty=False, return_theta = False)
			Deltaf_ij = fij[0]
			#Deltafid_ij = fij_id[0]
		else:
			(Deltaf_ij, dDeltaf_ij,theta) = mbar.getFreeEnergyDifferences(compute_uncertainty=True, return_theta = True)
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
			(Deltaf_ij, dDeltaf_ij,theta) = mbar.getFreeEnergyDifferences(compute_uncertainty=True, return_theta = True)
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

	for i in range(5*K+4*K):
		for j in range(5*K+4*K):
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

	print 1.0/mvola,UncDeltaf_ij[2*K+a,K+a], Deltaf_ij[2*K+a,K+a]

	mvolresa= (mvola*(m1+m2+m3)/1000 - 1000*kB*Ta/(Pa*1E5)) # residual volume m3/mol at standard state
	hresa = (Deltaf_ij[4*K+a,3*K+a])/(delbetaa) 
	uresa = hresa - (Pa*1E5*mvolresa)*1E-3 # kJ/mol
	#hresa = hresa - uresa
	h_t_a = (UncDeltaf_ij[4*K+a,3*K+a])/(delbetaa) 
	u_t_a = h_t_a - (Pa*1E5*mvola)*1E-3 # kJ/mol

	expgresa = getexp.getgresTP(Ta,Pa) - kB*Ta*numpy.log(getz(Ta,Pa)) # Gresa kJ/mol
	expsresa = getexp.getsresTP(Ta,Pa) + kB*numpy.log(getz(Ta,Pa)) # Sresa kJ/mol/K + kB*log(z) to convert from P(T,rho) id.g.st to V(T,P) EOS
	#expuresa = getexp.gethresTP(Ta,Pa) - Pa*1E5/1E3*getexp.getvresTP(Ta,Pa)   # Uresa kJ/mol/K
	expuresa = getexp.geturesTP(Ta,Pa)   # Uresa kJ/mol/K


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
		resden = 0 
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
			#exp_prop_array =  numpy.zeros([nTP,2+5+2], numpy.float64)
			exp_prop_array =  numpy.zeros([nTP,2+5+2+2], numpy.float64)
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

				#additional
				exp_hres = getexp.gethresTP(T,P)
				exp_alpha = getexp.getalphaTP(T,P)		

				exp_prop_array[j] =numpy.array([T,P,exprho,expmvolres,expuresl,expsres_la,expg,expp_cpres,expp_kappa,exp_alpha,exp_hres], numpy.float64)  

			numpy.save(exp_prop_file,exp_prop_array)
		else:
			exp_prop_array = numpy.load(exp_prop_file)

		for j in range(nTP): # goes over all T and P
			
			l = starting_k + j #+ K*btp 
			T = Pveccomb[l,6]
			P = Pveccomb[l,7]
		
			# Retrieving experimental properties from array
			exprho     = exp_prop_array[j,2]
			expmvolres = exp_prop_array[j,3]
			expuresl   = exp_prop_array[j,4]
			expsres_la = exp_prop_array[j,5]
			expg       = exp_prop_array[j,6]	
			exp_cpres  = exp_prop_array[j,7]	
			exp_kappa  = exp_prop_array[j,8]	
			exp_alpha  = exp_prop_array[j,9]	
			exp_hres   = exp_prop_array[j,10]	
		
			# calculating vres
			beta = 1.0/(kB*T) 
			Pr = P*1E5
			mvol = ((UncDeltaf_ij[2*K+l,K+l]/beta/(2*delp*1E5))*1000) # m3/mol
			mvolres= (mvol - 1000*kB*T/(P*1E5)) # residual volume m3/mol
			resvres = resvres + ((mvolres - expmvolres)*1E9)**2 # sum errors in Vres (mm3/mol)^2
			
			# calculating density
			dens = 0.0180154/mvol # density in kg/m3
			resden = resden +  (dens - exprho)**2
		
			
				
			# calculating hres  
			beta1 = beta+dbeta
			beta2 = beta-dbeta
			delbeta = beta1-beta2
			hres = (Deltaf_ij[4*K+l,3*K+l])/(delbeta) 
			reshres = reshres + (hres-exp_hres)**2


			#calculating h_t total enthalpy
			h_t = (UncDeltaf_ij[4*K+l,3*K+l])/(delbeta)

			# calculating ures = hres - PVres
			uresl = hres - (P*1E5*mvolres)*1E-3  # kJ/mol
			resures = resures + (uresl-expuresl)**2  # Sum errors in Ures (kJ/mol)^2

			#calculating u_t total internal energy
			u_t = h_t - (P*1E5*mvol)*1E-3  # kJ/mol

			# calculating Sres_la
			sres_la = (kB*(-Deltaf_ij[l,a] + (betaa*hresa - beta*hres)))*1000 # Sres,a - Sres,i  From force field in J/mol/K
			ressres = ressres + (sres_la -expsres_la)**2 # Sum errors in Sres (J/mol/K)^2

			#calculating St_la
			s_t_la = (kB*(-UncDeltaf_ij[l,a] + (betaa*h_t_a - beta*h_t)))*1000 # Sres,a - Sres,i  From force field in J/mol/K
			

			# calculating Gres_la
			ffpg = Deltaf_ij[l,a]# simulated residual free energy
			resg = resg + (Deltaf_ij[l,a]-expg)**2 # sum errors in Gres

			# calculating G_t_la
			ffpg_t = UncDeltaf_ij[l,a]# simulated residual free energy
 
		
			# calculating Cpres
			cpres = (-kB*(beta**2)*((Deltaf_ij[0+l,3*K+l]+Deltaf_ij[0+l,4*K+l])/(0.5*delbeta)**2))*1000  # in J/mol/K
			rescpres = rescpres + (cpres - exp_cpres)**2 # in (J/mol/K)^2

			# calculating Cp
			cp_t = (-kB*(beta**2)*((UncDeltaf_ij[0+l,3*K+l]+UncDeltaf_ij[0+l,4*K+l])/(0.5*delbeta)**2))*1000  # in J/mol/K

			# calculating kappa
			kappa = ((UncDeltaf_ij[2*K+l,0+l] + UncDeltaf_ij[K+l,0+l])/(delp*1E5*1E-5)**2)/( UncDeltaf_ij[2*K+l,K+l]/(2*delp*1E5*1E-5) ) # in 1/bar
			reskappa = reskappa + ((kappa - exp_kappa)*1E5)**2 # in 1/pa^2
	


			# calculating alpha Isobaric thermal expansion coefficient 1/K
			alp = (-kB*(-UncDeltaf_ij[2*K+l,K+l]/(2*delp*1E5)*beta+beta*beta*((UncDeltaf_ij[6*K+l,5*K+l]-UncDeltaf_ij[8*K+l,7*K+l])/(2*delp*1E5)/(delbeta))))/((UncDeltaf_ij[2*K+l,K+l]/(2*delp*1E5)))       
			resalpha = resalpha + (alp-exp_alpha)**2 # in 1/K^2


			if pvar ==1:
				
				#error in Gres_la kJ/mol
				dffpg = dDeltaf_ij[l,a]

				#error in density
				#ddens = (m1+m2+m3)/1000/((UncDeltaf_ij[2*K+st,K+st]/beta/(2*delp*1E5))*1000)/((UncDeltaf_ij[2*K+st,K+st]/beta/(2*delp*1E5))*1000)*((dDeltaf_ij[2*K+st,K+st]/beta/(2*delp*1E5))*1000)
				# error in mvolres	
				dmvolres = ((dDeltaf_ij[2*K+l,K+l]/beta/(2*delp*1E5))*1000) # error in estimates residual volume m3/mol
				
				#error in density
				ddens = 0.0180154*dmvolres/mvol**2 # kg/m3

				#error in ures # kJ/mol
				mu = 1.0/delbeta
				nu = P/beta/(2*delp)
				
				covu =   mu**2*theta[3*K+l,3*K+l] - mu**2*theta[3*K+l,4*K+l] - mu*nu*theta[3*K+l,K+l] +mu*nu*theta[3*K+l,2*K+l] \
					-mu**2*theta[4*K+l,3*K+l] + mu**2*theta[4*K+l,4*K+l] + mu*nu*theta[4*K+l,K+l] -mu*nu*theta[4*K+l,2*K+l] \
					-mu*nu*theta[K+l,3*K+l]   + mu*nu*theta[K+l,4*K+l]   + nu**2*theta[K+l,K+l]   -nu**2*theta[K+l,2*K+l] \
					+mu*nu*theta[2*K+l,3*K+l] - mu*nu*theta[2*K+l,4*K+l] - nu**2*theta[2*K+l,K+l] +nu**2*theta[2*K+l,2*K+l]
				
				err_ures = numpy.sqrt(numpy.abs(covu))	# error in uresl  From force field in kJ/mol

				# error in hres kJ/mol
				dhres = dDeltaf_ij[4*K+l,3*K+l]/(delbeta)

				# Error in entropy Sres_la
				ms = betaa/delbeta
				ns = -beta/delbeta
				
				covs = ms**2*theta[3*K+a,3*K+a] - ms**2*theta[3*K+a,4*K+a] + ms*ns*theta[3*K+a,3*K+l] - ms*ns*theta[3*K+a,4*K+l] - ms*theta[3*K+a,a] + ms*theta[3*K+a,l] \
				      -ms**2*theta[4*K+a,3*K+a] + ms**2*theta[4*K+a,4*K+a] - ms*ns*theta[4*K+a,3*K+l] + ms*ns*theta[4*K+a,4*K+l] + ms*theta[4*K+a,a] - ms*theta[4*K+a,l] \
				      +ms*ns*theta[3*K+l,3*K+a] - ms*ns*theta[3*K+l,4*K+a] + ns**2*theta[3*K+l,3*K+l] - ns**2*theta[3*K+l,4*K+l] - ns*theta[3*K+l,a] + ns*theta[3*K+l,l] \
				      -ms*ns*theta[4*K+l,3*K+a] + ms*ns*theta[4*K+l,4*K+a] - ns**2*theta[4*K+l,3*K+l] + ns**2*theta[4*K+l,4*K+l] + ns*theta[4*K+l,a] - ns*theta[4*K+l,l] \
				             -ms*theta[a,3*K+a] +        ms*theta[a,4*K+a] -           theta[a,3*K+l] +        ns*theta[a,4*K+l] +        theta[a,a] -        theta[a,l] \
					     +ms*theta[l,3*K+a] -	 ms*theta[l,4*K+a] +        ns*theta[l,3*K+l] -        ns*theta[l,4*K+l] -        theta[l,a] +        theta[l,l]	 	

				err_sla = 1000*kB*numpy.sqrt(numpy.abs(covs))	# error in Sres,a - Sres,i  From force field in J/mol/K

				# Error in kappa
				KX = (UncDeltaf_ij[2*K+l,l]+UncDeltaf_ij[K+l,l])
				KY = (UncDeltaf_ij[2*K+l,K+l])
				varKX =  theta[l,l]+theta[2*K+l,2*K+l]-2*theta[l,2*K+l] \
					+ theta[l,l]+theta[K+l,K+l]-2*theta[l,K+l] \
					+2*(theta[l,l]-theta[l,K+l]-theta[2*K+l,l]+theta[2*K+l,K+l])
				varKY = theta[2*K+l,K+l]
				covKXY = 2*theta[l,K+l]-2*theta[l,2*K+l]-theta[2*K+l,K+l]+theta[2*K+l,2*K+l]-theta[K+l,K+l]+theta[K+l,2*K+l]
				Kerr = 2/(delp*1E5)*numpy.sqrt(numpy.abs(varKX*1.0/KY**2+(varKY)*(KX/KY**2)**2+(2*covKXY)*(1.0/KY)*(-KX/KY*3)))


				#Error in Cpres Heat capacity
				covcp = theta[3*K+l,3*K+l]-2*theta[3*K+l,l]+theta[3*K+l,4*K+l] \
					+ 2*(-theta[l,3*K+l]+2*theta[l,l]-theta[l,4*K+l]) \
					+theta[4*K+l,3*K+l]-2*theta[4*K+l,l]+theta[4*K+l,4*K+l]
				cpreserr = kB*(beta**2)/((0.5*delbeta)**2)*1000*(numpy.sqrt(covcp))
				
				# Error in alpha Isobaric thermal expansion coefficient             
				alpX = UncDeltaf_ij[6*K+l,5*K+l]-UncDeltaf_ij[8*K+l,7*K+l]
				alpY = UncDeltaf_ij[2*K+l,K+l]  
				varalpX = theta[5*K+l,5*K+l]-theta[5*K+l,6*K+l]-theta[5*K+l,7*K+l]+theta[5*K+l,8*K+l] \
					 -theta[6*K+l,5*K+l]+theta[6*K+l,6*K+l]+theta[6*K+l,7*K+l]-theta[6*K+l,8*K+l] \
					  -theta[7*K+l,5*K+l]+theta[7*K+l,6*K+l]+theta[7*K+l,7*K+l]-theta[7*K+l,8*K+l] \
					 +theta[8*K+l,5*K+l]-theta[8*K+l,6*K+l]-theta[8*K+l,7*K+l]+theta[8*K+l,8*K+l]  
				varalpY = theta[K+l,K+l] + theta[2*K+l,2*K+l] -2*theta[K+l,2*K+l]
				covalpXY = theta[K+l,5*K+l]-theta[K+l,6*K+l]-theta[K+l,7*K+l]+theta[K+l,8*K+l] \
						    -theta[2*K+l,5*K+l]+theta[2*K+l,6*K+l]+theta[2*K+l,7*K+l]-theta[2*K+l,8*K+l]
				covaralp  = varalpX*(1.0/alpY)**2 + varalpY*(alpX/alpY**2)**2 + 2*covalpXY*(1.0/alpY)*(-alpX/alpY**2)
				alperr = kB*beta**2/delbeta*numpy.sqrt(numpy.abs(covaralp))

	


			out6 = "%8.2e %8.5e \
				%8.8e %8.8e %8.8e %8.8e %8.8e %8.8e\
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e \
				%8.8e %8.8e %8.8e %8.8e %8.8e\
				\n"%(T,P,\
				x[0],x[1],x[2],x[3],x[4],x[5],\
				dens,ddens,exprho,dens-exprho,\
				mvolres*1E9,dmvolres*1E9,expmvolres*1E9,(mvolres-expmvolres)*1E9,\
				uresl,err_ures,expuresl,uresl-expuresl,\
				sres_la,err_sla,expsres_la,sres_la-expsres_la,\
				ffpg,dffpg,expg,ffpg-expg,\
				cpres,cpreserr,exp_cpres,cpres-exp_cpres,\
				kappa,Kerr,exp_kappa,kappa-exp_kappa,\
				alp,alperr,exp_alpha,alp-exp_alpha,\
				hres,dhres,exp_hres,hres-exp_hres,\
				h_t,u_t,s_t_la,ffpg_t,cp_t)
			filename6 = scriptdir+"/loggopal22"+"_krsna_detailed"
			f6 = open(filename6, 'a')
			f6.write(out6)
			f6.close()


		resgres = resg /nTP
		resvres = resvres /nTP
		resures = resures /nTP
		ressres = ressres /nTP
		rescpres = rescpres /nTP#res0 =    res0/nTP
		reskappa = reskappa /nTP#res0 =    res0/nTP
		resalpha = resalpha /nTP#res0 =    res0/nTP
		resden   = resden /nTP#res0 =    res0/nTP
		reshres   = reshres /nTP#res0 =    res0/nTP
		
		return(resvres,resures,ressres,resgres,rescpres,reskappa,resalpha,resden,reshres)
		

	[resVres,resUres,resSres,resGres,resCpres,resKappa,resAlpha,resDen,resHres] = objfunc(nff*(NTP+NTPUS),NTP+NTPUS) 

	#sfval = sfval+resGres
	# International skeleton tables IST85 suggests allowable tolerances in properties of water v  = 0.01% h = 0.1%
	# MSE in v = (0.0001)**2 MSE in h = (0.001)**2 = 1E-6 no tolerances given for s and g so use the tolerance of h 
	#sc1 = sc1 + resUres/1000.0 #- 1E-5
	#sc2 = sc2 + resSres #- 1E-5
	#sc3 = sc3 + resVres #- 83#1E-1
	#sc4 = sc4 + resCpres #- 5#1E-1
	#sc5 = sc5 + resKappa #- 1E-5
	return(resGres,resUres,resSres,resVres,resCpres,resKappa,resAlpha,resDen,resHres)	

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


