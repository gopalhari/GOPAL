#!/usr/bin/python

import os
import numpy
import subprocess
import pdb
import math
import time
import getexp
#import expprop


	


def makexvg(XPvec,Xj,UPvec,Ui,nstxout,numTP,dirin,delp,dbeta,numproc,runnum,jobdirin,etype):
	#dir = "/home/hp4t/water_p_sweep/opt/"
	kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)
	pa2kJpmol = (1e5)*(1e-09)**3 * 6.02214 * (1e23) / 1000 # convert Panm3 to kJ/mol
	#delp = 0.9
	#dbeta = 1E-4
	 
	precision='single'
	gloc = "/home/hpaliwal/soft/gromacs462files/"

	os.environ['GMXLIB'] = "/home/hpaliwal/soft/gromacs462files/"
	os.environ['OMP_NUM_THREADS'] = "1"
	fresh = 'yes'
	#fresh = 'no'
	#dir = dirin
	scriptdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/"
	trr2xtc =scriptdir+"trr2xtc"

	pi = 3.14159265359
	pe = 5
	te = 7
	pv = 17
	vol = 15

	enedir = dirin+"/ene1/"
	trrdir = dirin+"/trr1/"
	mdpdir = dirin+"/mdp1/"
	topdir = dirin+"/top1/"
	tprdir = dirin+"/tpr1/"
	trrdir = dirin+"/trr1/"
	grodir = dirin+"/gro1/"
	cptdir = dirin+"/cpt1/"
	logdir = dirin+"/log1/"
	dgdldir = dirin+"/dgdl1/"
	mdp = scriptdir+"run.mdp"
	inputgro = scriptdir+"prod0.gro"
	templatetop = scriptdir+"/Lsolv_template.top"
	baseshell = scriptdir+"/shell_base.sh"

	m1=15.9994;
	m2=1.008;
	m3=1.008;
	M = m1+m2+m3;
	Nmol = 521.0

	if (precision == 'double'):
		grompp = "grompp_d"
		mdrun="mdrun_d -pin off"
		g_energy = "g_energy_d"
	else:
		grompp = "grompp"
		mdrun="mdrun -pin off"
		g_energy = "g_energy"
	def getmapJ(OHi,OHj,HHi,HHj):
		Kij = OHi/OHj*(numpy.sqrt((1.0-(HHi/2.0/OHi)*(HHi/2.0/OHi))/(1.0-(HHj/2.0/OHj)*(HHj/2.0/OHj))))
		bc = -m1/M

		l11ij = m1/M+(1.0+bc)*Kij
		l12ij = m2/M-(1.0+bc)*Kij/2.0
		l13ij = m3/M-(1.0+bc)*Kij/2.0
		l21ij = m1/M+bc*Kij
		l22ij = m2/M-bc*Kij/2.0+0.5*HHi/HHj
		l23ij = m3/M-bc*Kij/2.0-0.5*HHi/HHj
		l31ij = m1/M+bc*Kij
		l32ij = m2/M-bc*Kij/2.0-0.5*HHi/HHj
		l33ij = m3/M-bc*Kij/2.0+0.5*HHi/HHj

		J2 = l11ij*l22ij*l33ij-l11ij*l23ij*l32ij-l12ij*l21ij*l33ij+l12ij*l31ij*l23ij+l13ij*l21ij*l32ij-l13ij*l31ij*l22ij
		return(3*Nmol*numpy.log(numpy.fabs(J2)))


	def dotrr2xtc(in1_inputtrr,in2_reruntrr,in3_romlami,in4_rohlami,in5_rohlamj,in6_anglami,in7_anglamj):
		# Mapping trajecotry from jth state  to ith state
		if (os.path.exists(in2_reruntrr)):
			remtrr = "rm "+in2_reruntrr
			os.system(remtrr)	 
		outtranscmd = trr2xtc+" -i "+in1_inputtrr+" -o "+in2_reruntrr+" -romi "+str(in3_romlami)+" -rohi "+str(in4_rohlami)+" -rohj "+str(in5_rohlamj)+" -angi "+str(in6_anglami)+" -angj "+str(in7_anglamj)
		return(outtranscmd)

	def dorerun(in1_tprrerunfile,in2_reruntrr,in3_outtrr,in4_outgro,in5_outlog,in6_outedr):	 
		# Reevaluating the jth state samples with ith state potential definitions
		outrerunmd =  gloc+mdrun+" -nt "+str(numproc)+" -s "+in1_tprrerunfile+" -rerun "+in2_reruntrr+" -o "+in3_outtrr+" -c "+in4_outgro+" -g "+in5_outlog+" -e "+in6_outedr
		#print outrerunmd
		return(outrerunmd)

	def doextractPEpv(in1_i,in2_j,in3_outedr,in4_tprfile):
		# extract pe and pv from the outedr, rerun energy file
		energyfile = jobdirin+"RerunU"+str(int(in1_i))+"X"+str(int(in2_j))+"R"+str(runnum)+".xvg"
		outgetene = "echo "+str(pe)+" "+str(vol)+" | "+gloc+g_energy+" -f "+in3_outedr+" -s "+in4_tprfile+" -o "+energyfile
		#outgetene = "echo "+str(te)+" "+str(vol)+" | "+gloc+g_energy+" -f "+in3_outedr+" -s "+in4_tprfile+" -o "+energyfile
		return(outgetene)
	#print "XPvec",XPvec	
	#print "UPvec",UPvec
	#print Ui,Xj
	# Read parameters from UPvec
	rohlami = UPvec[0]
	rhhlami = UPvec[1]
	romlami = UPvec[2]
	sinthetaby2i = rhhlami/2/rohlami;
	anglami = math.asin(sinthetaby2i)*180/numpy.pi*2;
	costhetaby2i = math.sqrt(1-sinthetaby2i*sinthetaby2i);

	alami = romlami/(2*rohlami*costhetaby2i);
	blami = alami

	qolami = UPvec[3]
	qhlami = -qolami/2
	C6olami =UPvec[4]
	C12olami=UPvec[5]
	sigmaoi = (C12olami/C6olami)**(1.0/6.0)
	#C6olami = 4*UPvec[5]*(UPvec[4])**6.0
        #C12olami = 4*UPvec[5]*(UPvec[4])**12.0

	Tlami =  UPvec[6] # Temperature at which the simulation was run
	Plami =  UPvec[7] # Pressure at which the simulation was run
	i = Ui # sequence number of the water like model potential definition

	# Read parameters from XPvec
	rohlamj = XPvec[0]
	rhhlamj = XPvec[1]
	romlamj = XPvec[2]
	sinthetaby2j = rhhlamj/2/rohlamj;
	anglamj = math.asin(sinthetaby2j)*180/numpy.pi*2;
	costhetaby2j = math.sqrt(1-sinthetaby2j*sinthetaby2j);

	alamj = romlamj/(2*rohlamj*costhetaby2j);
	blamj = alamj

	qolamj = XPvec[3]
	qhlamj = -qolamj/2
	C6olamj =XPvec[4]
	C12olamj=XPvec[5]
	sigmaoj = (C12olamj/C6olamj)**(1.0/6.0)

	#C6olamj = 4*XPvec[5]*(XPvec[4])**6.0
        #C12olamj = 4*XPvec[5]*(XPvec[4])**12.0


	Tlamj =  XPvec[6] # Temperature at which the simulation was run
	Plamj =  XPvec[7] # Pressure at which the simulation was run
	j = Xj # sequence number of the state with samples

	nstenergy = nstxout #since frequency energy = frequency of samples 
	
	# Make a new top file with the ith state topology parameters
	newtop = jobdirin+"LsolvU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".top"      
	copytop = "cp "+templatetop+" "+newtop
	os.system(copytop)
	
	# Make a new mdp file with ith state Temperature and Pressure
	newenemdp = jobdirin+"rerunU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".mdp"      
	copyenemdp = "cp "+mdp+" "+newenemdp
	os.system(copyenemdp)
	#pdb.set_trace()
	
	topcmd1 = "sed -i -e 's/fillrohlam/"+str(rohlami)+"/' "+newtop
	topcmd2 = "sed -i -e 's/fillrhhlam/"+str(rhhlami)+"/' "+newtop
	topcmd3 = "sed -i -e 's/fillromlam/"+str(romlami)+"/' "+newtop
	topcmd4 = "sed -i -e 's/fillalam/"+str(alami)+"/' "+newtop
	topcmd5 = "sed -i -e 's/fillblam/"+str(blami)+"/' "+newtop
	topcmd6 = "sed -i -e 's/fillqolam/"+str(qolami)+"/' "+newtop
	topcmd7 = "sed -i -e 's/fillqhlam/"+str(qhlami)+"/' "+newtop
	topcmd8 = "sed -i -e 's/fillC6olam/"+str(C6olami)+"/' "+newtop
	topcmd9 = "sed -i -e 's/fillC12olam/"+str(C12olami)+"/' "+newtop

	kappa = getexp.getkappaTP(Tlami,Plami)[0] # isothermal compressibility at that temperature and pressure
	mdpcmd13 = "sed -i -e 's/filltemperature/"+str(Tlami)+"/' "+newenemdp
	mdpcmd14 = "sed -i -e 's/fillpressure/"+str(Plami)+"/' "+newenemdp
	mdpcmd15 = "sed -i -e 's/fillnsteps/"+str(0)+"/' "+newenemdp
	mdpcmd16 = "sed -i -e 's/fillnstxout/"+str(0)+"/' "+newenemdp
	mdpcmd17 = "sed -i -e 's/fillnstvout/"+str(0)+"/' "+newenemdp
	mdpcmd18 = "sed -i -e 's/fillnstenergy/"+str(nstenergy)+"/' "+newenemdp
	mdpcmd19 = "sed -i -e 's/fillkappa/"+str(kappa)+"/' "+newenemdp
	
	getnewtop = topcmd1+";"+topcmd2+";"+topcmd3+";"+topcmd4+";"+topcmd5+";"+topcmd6+";"+topcmd7+";"+topcmd8+";"+topcmd9
	getnewmdp = mdpcmd13+";"+mdpcmd14+";"+mdpcmd15+";"+mdpcmd16+";"+mdpcmd17+";"+mdpcmd18+";"+mdpcmd19
	
	tprrerunfileI = jobdirin+"prodrerunU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".tpr"
	
	prodgro   = grodir+"prodU"+str(int(j))+"X"+str(int(j))+".gro" # since some i states are not sampled. 
	
	genreruntpr =  gloc+grompp+" -f "+newenemdp+" -po "+jobdirin+"mdrerunout"+"R"+str(runnum)+".mdp"+" -c "+prodgro+" -o "+tprrerunfileI+" -p "+newtop+" -maxwarn 4"

        gettprI = getnewtop+";"+getnewmdp+";"+genreruntpr
	#os.system(gettprI)
	
	tprrerunfileJ = tprdir+"prodrerun"+str(int(j))+".tpr" # Potential definition from the jth state should be present already
	prodedrJ   = enedir+"prodU"+str(int(j))+"X"+str(int(j))+".edr" # energy file from the jth simulation

	inputtrr = trrdir+"prodU"+str(int(j))+"X"+str(int(j))+".trr" # samples from the jth simulation

	reruntrr = jobdirin+"transU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".trr" # mapped trajecory Tij(xj) 
	outtrr   = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".trr"
	outgro   = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".gro"
	outedr   = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".edr"
	outlog   = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".log"
	outdgdl  = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".xvg"
	outcpt   = jobdirin+"outU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".cpt"
         
	energyfile = jobdirin+"RerunU"+str(int(i))+"X"+str(int(j))+"R"+str(runnum)+".xvg"

	# Mapping ,rerun and energy extraction done according to the need to reduce unecessary time consumption  
	
	#pdb.set_trace()
	if etype == 'simp':
		if fresh == 'yes':
			if (os.path.exists(prodedrJ)):
				combine3 = doextractPEpv(i,j,prodedrJ,tprrerunfileJ)
			else:
				rerunmd =  dorerun(tprrerunfileI,inputtrr,outtrr,outgro,outlog,prodedrJ)
				getene = doextractPEpv(i,j,prodedrJ,tprrerunfileI)
				combine3 =gettprI+";"+rerunmd+";"+getene
			

		else:
			#print "no file is being generated only energy being read from the ",j," state"
			combine3 = "echo "+"Reading..."+energyfile

	else: # mapping rerun and extraction required
		transcmd = dotrr2xtc(inputtrr,reruntrr,romlami,rohlami,rohlamj,anglami,anglamj) # mapping
		rerunmd =  dorerun(tprrerunfileI,reruntrr,outtrr,outgro,outlog,outedr) # rerun
		getene = doextractPEpv(i,j,outedr,tprrerunfileI)
		if fresh == 'yes' :
			#if (rohlami == rohlamj) and (rhhlami==rhhlamj) and (romlami==romlamj): # no mapping required
			#	combine3 =rerunmd+";"+getene
			#else:
			combine3 =transcmd+";"+gettprI+";"+rerunmd+";"+getene
		else:
			combine3 = "echo "+"Reading..."+energyfile
			##else: # use the first of numTP states edr to get energies and reweight with different T and P
			##	if fresh == 'yes':
			##		getene = doextractPEpv(i,j,onlyTPedr,tprrerunfileI)
			##		combine3 =getene
			##	else:
			##		combine3 = "echo "+"Reading..."+energyfile	

		#print combine3	
		#pdb.set_trace()
	
		
	#os.system(combine3)

	# for running NOMAD 
	#proc = subprocess.Popen([combine3], stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	#(comout,comerr) = proc.communicate()
	#proc.stdout.close()
	#proc.stdin.close()
	#proc.stderr.close()
	#pdb.set_trace()
	proc = subprocess.Popen([combine3], stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	#print proc.stdout.read()
	#print proc.stderr.read()
	(comout,comerr) = proc.communicate()
	#print comerr
	proc.stdout.close()
	proc.stdin.close()
	proc.stderr.close()

	#kvij = (sigmaoi/sigmaoj)**3.0
	kvij = 1.0
        while not (os.path.exists(energyfile)):
		#print(i,j,'wating...\n')
		time.sleep(2)

	uandv = numpy.genfromtxt(energyfile, skip_header=21)
	Uij = uandv[:,1]
	Vij = kvij*uandv[:,2]
	#pv = uandv[:,3]
	#print pv[0],Plami,Plami*pa2kJpmol*v[0]
	#pdb.set_trace()
	#betai = 1/kB/Tlami
        #redu = betai*(Uij + Plami*Vij*pa2kJpmol  ) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj)  
        #redu1 = betai*(Uij + (1+delp)*Plami*Vij*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu2 = betai*(Uij + (1-delp)*Plami*Vij*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu3 = (1+dbeta)*betai*(Uij + Plami*Vij*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu4 = (1-dbeta)*betai*(Uij + Plami*Vij*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu5 = (1+dbeta)*betai*(pe + (1+delp)*Plami*v*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu6 = (1-dbeta)*betai*(pe + (1+delp)*Plami*v*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu7 = (1+dbeta)*betai*(pe + (1-delp)*Plami*v*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
        #redu8 = (1-dbeta)*betai*(pe + (1-delp)*Plami*v*pa2kJpmol) - getmapJ(rohlami,rohlamj,rhhlami,rhhlamj) 
	
	#print getmapJ(rohlami,rohlamj,rhhlami,rhhlamj)
	
	#removefiles ="rm "+" "+outtrr+" "+outgro+" "+outlog+" "+outcpt+" "+mdpdir+"/*#* "+enedir+"/*#* "+trrdir+"/*#* "+tprdir+"/*#* "+jobdir+"/*#* "+jobdir+"/*.trr* "+jobdir+"/*.top* "+jobdir+"/*tpr* "+jobdir+"/*.mdp* "+jobdir+"/*.xvg* "+jobdir+"/*.cpt* "+jobdir+"/*.gro* "+jobdir+"/*.edr* "+jobdir+"/*#* "
	removefiles ="rm "+jobdirin+"/*#* "
	
	#Turn off the file deletion for making all the energy files with fresh = yes and then turn fresh = no to just read the energy files
	# Turn on the file deletion once the debugging is done.
	#os.system(removefiles)
	proc = subprocess.Popen([removefiles], stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	(comout,comerr) = proc.communicate()
	proc.stdout.close()
	proc.stdin.close()
	proc.stderr.close()

	Jij = kvij*getmapJ(rohlami,rohlamj,rhhlami,rhhlamj)	
	#return(redu,redu1,redu2,redu3,redu4,redu5,redu6,redu7,redu8,pe,v*pa2kJpmol,Jij)
	#return(redu,redu1,redu2,redu3,redu4,Uij,Vij*pa2kJpmol,Jij)
	return(Uij,Vij*pa2kJpmol,Jij)
