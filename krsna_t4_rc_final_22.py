#!/usr/bin/python


import numpy
import pdb
import os
import sys
import math
import time
###from scipy import interpolate
###import scipy.optimize as spop
## Import modules 
#import geneinitPvec 
import genPvec_v2 
import gentrr
import getenev3bnew_pp_a3
import getexp
import old_mbar 
import pp
import get_simprop2

#pdb.set_trace()

t0 = time.time()
x = numpy.genfromtxt(sys.argv[1], dtype=float)
#pinput = numpy.genfromtxt(sys.argv[2], dtype=float)
scriptdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/0of2/"
##pythonexec = "/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python"
##nodefile = scriptdir+"unode"
##port = ':'+str(35000)
##nodelist = open(nodefile,'r').read().splitlines()
##nodewithport = [nstr + port for nstr in nodelist]
##nodes=tuple(nodewithport)
##
##plist = [1,100,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]# 15
#plist = [1,100,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000]# 16
#plist = [1,100]# 16
plist = [1]# 15

sfval = 0 
sc1 = 0 
sc2 = 0 
sc3 = 0 
sc4 = 0 
sc5 = 0
sc6 = 0
sc7 = 0
sc8 = 0
jobsn = []
Njobs=0
ppservers=()
ncpus =int(40)
#job_server = pp.Server(ncpus=0,ppservers=nodes,secret="krsna")
job_server = pp.Server(ncpus=ncpus,ppservers=ppservers)

scroh = 1.0
scrhh = 1.0
scrom = 1.0
scq = 1.0
scc6 = 1.0
scc12 = 1.0	

smat = numpy.array([scroh,scrhh,scrom,scq,scc6,scc12], numpy.float64)

#########################################################################
#gopal15 # spce + opc + tip4p-ew + gopal3 + tip3p
#########################################################################
##fstguess = numpy.array([[0.1,0.163300,0.0,-0.8476,0.0026188975686492,2.6354497260417803e-06],[0.08724,0.137116310721,0.01594,-1.3582,0.003590505607210,3.62005237431366E-06],[0.09572,0.151390,0.0125,-1.04844,0.002735499717669941,2.7465336819642542e-06],[8.95656663e-02,1.41094825e-01,1.51767981e-02,-1.27044796e+00,3.37613292e-03,3.39225806e-06],[0.09789,0.15985249656554,0,-0.8952,2.79757931E-03,2.86184161E-06]], numpy.float64)#0.08956567 1.81400774 0.01517680 -1.27044796 0.31647900 0.84002111
#########################################################################
#########################################################################
#gopal17 # spce + opc + tip4p-ew + gopal4 + opc3
#########################################################################
#fstguess = numpy.array([[0.1,0.163300,0.0,-0.8476,0.0026188975686492,2.6354497260417803e-06],[0.08724,0.137116310721,0.01594,-1.3582,0.003590505607210,3.62005237431366E-06],[0.09572,0.151390,0.0125,-1.04844,0.002735499717669941,2.7465336819642542e-06],[8.93580228e-02,1.41166797e-01,1.51538350e-02,-1.27055475e+00,3.31129529e-03,3.31824835e-06],[0.09789,0.15985249656554,0,-0.8952,2.79757931E-03,2.86184161E-06]], numpy.float64)#0.08935802 1.82127344 0.01515383 -1.27055475 0.31633834 0.82608920

#gopal19 # gopal_f + opc 
#########################################################################
#fstguess = numpy.array([[8.93580228e-02,1.41166797e-01,1.51538350e-02,-1.27055475e+00,3.31129529e-03,3.31824835e-06],[0.08724,0.137116310721,0.01594,-1.3582,0.003590505607210,3.62005237431366E-06]], numpy.float64)#0.08935802 1.82127344 0.01515383 -1.27055475 0.31633834 0.82608920
# gopal 22 #spce + opc + tip4p-ew + gopal4 + opc3 
fstguess = numpy.array([[0.1,0.163300,0.0,-0.8476,0.0026188975686492,2.6354497260417803e-06],[0.08724,0.137116310721,0.01594,-1.3582,0.003590505607210,3.62005237431366E-06],[0.09572,0.151390,0.0125,-1.04844,0.002735499717669941,2.7465336819642542e-06],[8.93580228e-02,1.41166797e-01,1.51538350e-02,-1.27055475e+00,3.31129529e-03,3.31824835e-06],[0.09789,0.15985249656554,0,-0.8952,2.79757931E-03,2.86184161E-06]], numpy.float64)#0.08935802 1.82127344 0.01515383 -1.27055475 0.31633834 0.82608920

theta = x[1] # in radians
rhh = 2*x[0]*math.sin(theta/2)
x[1] =rhh

sigLJ = x[4] # in nm
epsLJ = x[5] # in kJ/mol


x[4] = 4*epsLJ*sigLJ**6
x[5] = 4*epsLJ*sigLJ**12

rohlami = x[0]
rhhlami = x[1]
sinthetaby2i = rhhlami/2/rohlami;
#fstiter = numpy.array([(x[0])*1e-02,(x[1])*1e-01,(x[2])*1e-02,x[3],(x[4])*1e-03,(x[5])*1e-06], numpy.float64)
fstiter = numpy.divide(x,smat)

gamma = 1E-5 # how much we want to descend in the direction of the gradient 
kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)
#trrfac = 6  #17,19
trrfac = 6*4 #22
nsteps = trrfac*100000 # total number of steps in production run
nstxout = 600    # frequency of sampling 
pvar = 1 # 0 variance in property not calculated if != 0 variance calculated
wts = numpy.array([1.0,1.0,1.0,1.0,1.0,0.0], numpy.float64) 
numproc = 1
cluster = True # use cluster or not  
verbose = True
optimization = True # If True optimization carried out if false only property prediction at given T and P carried out using the firstguess water model

N = 521.0 # number of water molecules in the box
nbs = 1 # number of bootstrap 
etol = 1E-7#1E-8 # tolerance for the objective funtion
err = 10 # initialize error in the obejectve function 
tk = 1 # step length control paramter
count = 0 #initialize the count of the number of iteration.
maxcount = 1000 # maximum number of times iterations could be run
C2K = 273.15# celcius to Kelvin + 0.01 K because the NIST liquid phase data starts from 273.16 K

m1=15.9994 
m2=1.008 
m3=1.008 
M = m1+m2+m3
delp = 1E-3 # bar
dbeta = 1E-4

#### Checking if the paramter combination has been calculated before 

ptol = 1e-8
skip_prev_calc=0
docalc = 1

if skip_prev_calc == 1:
	#resfile = scriptdir+"/loggopal15_krsna_vns"
	resfile = scriptdir+"/log19_1"
	premat = numpy.genfromtxt(resfile, dtype=float)
	#pdb.set_trace()
	for i in range(len(premat)):
		if (docalc == 1 ):
			if (abs(x[0]-premat[i,0]) <= ptol) and (abs(x[1]-premat[i,1]) <= ptol) and (abs(x[2]-premat[i,2]) <= ptol) and (abs(x[3]-premat[i,3]) <= ptol) and (abs((x[4]-premat[i,4]))*1e3 <= ptol) and (abs((x[5]-premat[i,5]))*1e6 <= ptol):
				t0 = time.time()
				sfval= premat[i,6]
				sc1= premat[i,7]
				sc2= premat[i,8]
				sc3= premat[i,9]
				sc4= premat[i,10]
				sc5= premat[i,11]
				t2 = time.time()
				#theta = premat[i,1] # in radians
				#rhh = 2*premat[i,0]*math.sin(theta/2)
				np = 1 # since we are just copying the results
				out52 = "%8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %4.2f\n"%(x[0],x[1],x[2],x[3],x[4],x[5],sfval/np,sc1/np,sc2/np,sc3/np,sc4/np,sc5/np,t2-t0)
				filename = scriptdir+"/loggopal22"+"_krsna_test4_vns"
				f = open(filename, 'a')
				f.write(out52)
				f.close()
				docalc = 0

				print "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(sfval/np,sc1/np,sc2/np,sc3/np-90,sc4/np,sc5/np)
			else:
				pass
		else:
			pass

# Create fresh uklt, Pvec, fk
freshmats = 0

if docalc == 1:
	for pressure in plist:

		case = "gopal22"+"_"+str(pressure)
		runcase = "gopal22"+"_wd_"+str(pressure) # for directory and logfile
		#case = "t4ew"
		dir ="/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/"+case+"/"
		runnum = 0
		#jobdir = "/home/hp4t/optscan/jobdir/"+str(int(runnum))+"of"+str(int(2))+"/"
		#jobdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/"+str(int(runnum))+"of"+str(int(2))+"/mads2_gop22_2big_a1/"
		jobdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/"+str(int(runnum))+"of"+str(int(2))+"/mads2_t_"+runcase+"/"

		if (os.path.exists(jobdir)):
			pass
		else:
			makedir = "mkdir"+" "+jobdir
			os.system(makedir)

		if freshmats==1:
			ukltfile = jobdir+case+'_uklt_ext_cp_hk.npy'
			Pvecfile = jobdir+case+'_pveccomb_ext_cp_hk.npy'
			fkfile = jobdir+case+'_fk_ext_cp_hk.npy'
			expfile = jobdir+case+'_exp_prop_ext_cp_hk.npy'
			# Deleting old files
			rmuklt = "rm"+" "+ ukltfile
			rmPvec = "rm"+" "+ Pvecfile
			rmfk = "rm"+" "+ fkfile
			rmexp = "rm"+" "+ expfile
			os.system(rmuklt)
			os.system(rmPvec)
			os.system(rmfk)
			os.system(rmexp)
			##resetting freshmats
			freshmat=0

			
		Nk = int(numpy.float(nsteps)/numpy.float(nstxout))

		#Ts = C2K+numpy.array([1,4,10,20,30,40,50,60,70,80,90,99], numpy.float64) # sampled temperatures case 15
		#Ts = C2K+numpy.array([1,4,10,20,30,40,50,60,70,80,90,99], numpy.float64) # sampled temperatures case 16
		#Ts = C2K+numpy.array([1,4,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99], numpy.float64) # sampled temperatures case 17
		#Ts = C2K+numpy.arange(1,100,1) # sampled temperatures case 19
		Ts = C2K+numpy.array([1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99], numpy.float64) # sampled temperatures case 22

		Ps = 0.01325+numpy.array([pressure], numpy.float64) # sampled pressures case 13,14
		#pdb.set_trace()
		TsPs = numpy.zeros([len(Ts)*len(Ps),2], numpy.float64)
		c =0
		for i in range(len(Ps)):
			for j in range(len(Ts)):
				TsPs[c,0] = Ts[j]
				TsPs[c,1] = Ps[i]
				c = c+1 

		NTP = len(TsPs)
			
				
		Tus = C2K+numpy.array([6], numpy.float64) # unsampled temperatures original case 13
		#Tus = C2K+numpy.arange(1,100,1) # unsampled temperatures original case 13
		Pus = 0.01325+numpy.array([pressure], numpy.float64) # unsampled pressures

		####     This section makes sure that the Temperatures and pressures combinations are not repeated while generating unsampled T and P combinations.
		def presentq(ele):
			ans = False
			for i in range(len(Ts)*len(Ps)):
				if ele[0] == TsPs[i,0] and ele[1] == TsPs[i,1]:
					ans = True
					return(ans)
			return(ans) 

		c =0 	
		for i in range(len(Pus)):
			for j in range(len(Tus)):
				ele = numpy.array([Tus[j],Pus[i]], numpy.float64)
				if not presentq(ele):
					if c == 0 :
						TusPus = numpy.concatenate([[ele]])
					else: 
						TusPus = numpy.concatenate((TusPus,[ele]), axis=0)
					c = c+1	

		NTPUS = len(TusPus)

		#pdb.set_trace()
		# Generate the vector of parameters for guess water model, T and P (FF,T,P)
		Pvecs = genPvec_v2.genvec(fstguess,TsPs) #  sampled
		Pvecus = genPvec_v2.genvec(fstguess,TusPus) #  unsampled

		[at,bt] = numpy.shape(Pvecs)
		[aus,bus] = numpy.shape(Pvecus)
		nff = int(numpy.size(fstguess)/6) # number of water models used i.e. number of water models used

		#NumTP = 3*(NTP+NTPUS) 
		NumTP = nff*(NTP+NTPUS) 
		K = (nff+1)*(NTP+NTPUS) # total number of states = sampled force fields + iterated force field
		unsf1 = genPvec_v2.genvec(fstiter,TsPs)
		unsf2 = genPvec_v2.genvec(fstiter,TusPus)

		# combine all vectors in a single vector
		Pvec0 = numpy.concatenate((Pvecs,Pvecus),axis=0) 
		Pvec01 = numpy.concatenate((Pvec0,unsf1),axis=0) 
		Pvec02 = numpy.concatenate((Pvec01,unsf2),axis=0) 

		Pveccomb = Pvec02

		ref_st = NumTP # starting state for the iterated force field used as the reference state where G is set 0
		unsf3 = numpy.zeros([4*K+4*K,6+2], numpy.float64)	# [+dt, -dt, +dp, -dp, +dp-dt, +dp+dt, -dp+dt, -dp-dt; six water model parameters + 2 T and P] 
		Pvec03 = numpy.concatenate((Pvec02,unsf3),axis=0) 
		Pveccomb = Pvec03

		ukltfile = jobdir+case+'_uklt_ext_cp_hk.npy'
		Pvecfile = jobdir+case+'_pveccomb_ext_cp_hk.npy'
		fkfile = jobdir+case+'_fk_ext_cp_hk.npy'
		exp_prop_file = jobdir+case+'_exp_prop_ext_cp_hk.npy'
		
		Nk = int(numpy.float(nsteps)/numpy.float(nstxout))
		#pdb.set_trace()

		# If ukltfile is not present initialise a new array else load old ukltfile
		##if not (os.path.exists(ukltfile)):
		##	u_klt = numpy.zeros([at,K+4*K+4*K,Nk], numpy.float64)# +9K for +-dt,+-dp,+dp-dt,+dp+dt,-dp+dt,-dp-dt at each T,P point
		##else:
		##	u_klt = numpy.load(ukltfile)

		# If Pvecfile is not present save Pvecfile else load old Pvecfile
		if not (os.path.exists(Pvecfile)):
			numpy.save(Pvecfile,Pveccomb)
		else:
			pass
		numproc = 1
		for ki in range(0,nff,1): 
			for k in range(ki*(NTP),(ki+1)*NTP,1):
				#print (k)
				if not (os.path.exists(ukltfile)):
					for l in range(0,nff*NTP,NTP):
						if (l==k):
							etype='simp' # just extract energy
						else:
							etype='MRE' # just extract energy
						#print(k,l,Njobs)
						#if pressure ==5000:
							#print dir
							#pdb.set_trace()
						jobsn.append(job_server.submit(getenev3bnew_pp_a3.makexvg, (Pvecs[k],k,Pveccomb[l],l,nstxout,NumTP,dir,delp,dbeta,numproc,runnum,jobdir,etype), (), ("numpy","subprocess","pdb","math","os","time","getexp",)))
						Njobs = Njobs+1
					
				#if pressure ==5000:
				#	print dir
				#	pdb.set_trace()

				l = nff*(NTP+NTPUS) # job results will contain MRE1M4
				if verbose:	
					print(pressure,k,l)
				#print dir,jobdir
				#time.sleep(10)	
				etype = 'MRE' # Map rerun and extract energy
				jobsn.append(job_server.submit(getenev3bnew_pp_a3.makexvg, (Pvecs[k],k,Pveccomb[l],l,nstxout,NumTP,dir,delp,dbeta,numproc,runnum,jobdir,etype), (), ("numpy","subprocess","pdb","math","os","time","getexp",)))
				Njobs = Njobs+1


	#pdb.set_trace()
	cjobs=0
	Njobs = len(jobsn)
	while cjobs < Njobs:
		#time.sleep(5)
		cjobs=0
		for i in range(Njobs):
			if jobsn[i].finished:
				cjobs=cjobs+1
		if verbose:
			print cjobs,Njobs
			#pdb.set_trace()
			#time.sleep(10)	
		if cjobs < Njobs:
			compjobs=0

		
	#pdb.set_trace()
	resultsjn = []
	for jn in range(len(jobsn)):
		resultsjn.append(jobsn[jn]())

	##job_server.destroy()

	##job_server2 = pp.Server(ncpus=0,ppservers=nodes,secret="krsna")
	job_server2 = pp.Server(ncpus=len(plist),ppservers=ppservers)
		
	rjobs=0 #This has to stay out the for loop because all jobs are in one array. After rjobs for p1 ends rjobs for p2 pressure starts 

	Nk = int(numpy.float(nsteps)/numpy.float(nstxout))
	#job_server.destroy()

	#job_server = pp.Server(ncpus=0,ppservers=nodes,secret="krsna")
	jobsn2 = []
	tot_jobs = len(jobsn)
	iter = 0 
 
	jobperpress = int(tot_jobs/len(plist))
	for pressure in plist:
		case = "gopal22"+"_"+str(pressure)
		runcase = "gopal22"+"_wd_"+str(pressure) # for directory and logfile
		dir ="/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/"+case+"/"
		runnum = 0
		jobdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/optscan/jobdir2/"+str(int(runnum))+"of"+str(int(2))+"/mads2_t_"+runcase+"/"

		rjobs = iter*jobperpress # strating rjob has to be supplied to each task
		#pdb.set_trace()
		[resGres,resUres,resSres,resVres,resCpres,resKappa,resAlpha,resDen,resHres] = get_simprop2.getres(scriptdir,jobdir,case,rjobs,resultsjn,nff,at,K,Nk,NTP,NTPUS,NumTP,verbose,pvar,ref_st,x)
	###	jobsn2.append(job_server2.submit(get_simprop.getres, (scriptdir,jobdir,case,rjobs,resultsjn,nff,at,K,Nk,NTP,NTPUS,NumTP,verbose,pvar,ref_st,x), (), ("numpy","pdb","os","math","getexp","old_mbar",)))
		iter = iter+1
		sfval = sfval+resGres
		sc1 = sc1 + resUres/1000.0 #- 1E-5
		sc2 = sc2 + resSres #- 1E-5
		sc3 = sc3 + resVres #- 83#1E-1
		sc4 = sc4 + resCpres #- 5#1E-1
		sc5 = sc5 + resKappa #- 1E-5
		sc6 = sc6 + resAlpha #- 1E-5
		sc7 = sc7 + resDen #- 1E-5
		sc8 = sc8 + resHres #- 1E-5


	#pdb.set_trace()
	np = len(plist)
	t2 = time.time()
	out5 = "%8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e %4.2f\n"%(x[0],x[1],x[2],x[3],x[4],x[5],sfval/np,sc1/np,sc2/np,sc3/np,sc4/np,sc5/np,sc6/np,sc7/np,sc8/np,t2-t0)
	filename = scriptdir+"/loggopal22"+"_krsna_summary"
	f = open(filename, 'a')
	f.write(out5)
	f.close()

	#print "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(sfval/np,sc1/np,sc2/np,sc3/np-100,sc4/np-10,sc5/np)
	#print "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(sfval/np,sc4/np-10,sc5/np-0.05)
	print "%4.8e %4.8e %4.8e %4.8e %4.8e %4.8e\n" %(sfval/np,sc1/np,sc2/np,sc3/np-90,sc4/np,sc5/np)


