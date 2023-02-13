#!/usr/bin/python

import os
import numpy
import subprocess
import pdb
import math
import getexp

def maketrr(inputPvec,seqn,nstepsinp,nstxoutinp,cluster,dirin,case):
	
	numproc=1
	gloc = "/home/hpaliwal/soft/gromacs462files/"
	#precision = 'double'
	precision = 'single'
	os.environ['GMXLIB'] = "/home/hpaliwal/soft/gromacs462files/"
	os.environ['OMP_NUM_THREADS'] = "1"
	
	#dir = "/home/hp4t/water_p_sweep/opt/"
	dir = dirin
	scriptdir = "/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/water_p_sweep/scripts/"
	#trr2xtc = dir+"trr2xtc"
	trr2xtc = scriptdir+"trr2xtc"
	
	pi = 3.14159265359
	pe = 5
	te = 7
	pv = 17
	
	enedir = dir+"/ene1/"
	trrdir = dir+"/trr1/"
	mdpdir = dir+"/mdp1/"
	topdir = dir+"/top1/"
	tprdir = dir+"/tpr1/"
	trrdir = dir+"/trr1/"
	grodir = dir+"/gro1/"
	cptdir = dir+"/cpt1/"
	logdir = dir+"/log1/"
	dgdldir = dir+"/dgdl1/"
	minmdp = scriptdir+"min.mdp"
	eqmdp = scriptdir+"eq.mdp"
	mdp = scriptdir+"run.mdp"
	#inputgro = scriptdir+"prod0.gro"
	inputgro = scriptdir+"prod_gopal1.gro"
	templatetop = scriptdir+"/Lsolv_template.top"
	baseshell = scriptdir+"/shell_base.sh"
	
	def makedir(dirname):
		if (os.path.exists(dirname)):
			pass
		else:
			makedir = "mkdir"+" "+dirname
			os.system(makedir)
	
	# make directories if they don't exist
	makedir(dir)
	makedir(enedir)
	makedir(trrdir)
	makedir(mdpdir)
	makedir(topdir)
	makedir(tprdir)
	makedir(trrdir)
	makedir(grodir)
	makedir(cptdir)
	makedir(logdir)
	makedir(dgdldir)
	
	
	
	# Read parameters from inputPvec 
	rohlam = inputPvec[0]
	rhhlam = inputPvec[1]
	romlam = inputPvec[2]
        sinthetaby2 = rhhlam/2/rohlam;
        anglam = math.asin(sinthetaby2)*180/numpy.pi*2;
        costhetaby2 = math.sqrt(1-sinthetaby2*sinthetaby2);

        alam = romlam/(2*rohlam*costhetaby2);
        blam = alam

	qolam = inputPvec[3]
	qhlam = -qolam/2
	C6olam = inputPvec[4]
	C12olam = inputPvec[5]
	Tlam = inputPvec[6] # Temperature at which the simulation has to be run
	Plam = inputPvec[7] # Pressure at which the simulation has to be run
	i = seqn # sequence number of the water like model
	kappa = getexp.getkappaTP(Tlam,Plam)[0] # isothermal compressibility at that temperature and pressure
	
	nstepseq = 100000#50000
	nstxouteq =0
	nstvouteq =0
	nstenergyeq =0
	
	nsteps = nstepsinp
	nstxout = nstxoutinp
	nstvout = nstxout
	nstenergy = nstxout
	
	
	# Make a new top file with the input topology parameters
	newtop = topdir+"Lsolv"+str(int(i))+".top"      
	copytop = "cp "+templatetop+" "+newtop
	os.system(copytop)
	
	# Make a new minimization mdp file with the new Temperature and Pressure
	newminmdp = topdir+"runmin"+str(int(i))+".mdp"      
	copyminmdp = "cp "+minmdp+" "+newminmdp
	os.system(copyminmdp)
	
	# Make a new mdp file with the new Temperature and Pressure
	neweqmdp = topdir+"runeq"+str(int(i))+".mdp"      
	copyeqmdp = "cp "+eqmdp+" "+neweqmdp
	os.system(copyeqmdp)
	
	newmdp = topdir+"run"+str(int(i))+".mdp"      
	copymdp = "cp "+mdp+" "+newmdp
	os.system(copymdp)
	
	newenemdp = topdir+"rerun"+str(int(i))+".mdp"      
	copyenemdp = "cp "+mdp+" "+newenemdp
	os.system(copyenemdp)
	#pdb.set_trace()
	
	topcmd1 = "sed -i -e 's/fillrohlam/"+str(rohlam)+"/' "+newtop
	topcmd2 = "sed -i -e 's/fillrhhlam/"+str(rhhlam)+"/' "+newtop
	topcmd3 = "sed -i -e 's/fillromlam/"+str(romlam)+"/' "+newtop
	topcmd4 = "sed -i -e 's/fillalam/"+str(alam)+"/' "+newtop
	topcmd5 = "sed -i -e 's/fillblam/"+str(blam)+"/' "+newtop
	topcmd6 = "sed -i -e 's/fillqolam/"+str(qolam)+"/' "+newtop
	topcmd7 = "sed -i -e 's/fillqhlam/"+str(qhlam)+"/' "+newtop
	topcmd8 = "sed -i -e 's/fillC6olam/"+str(C6olam)+"/' "+newtop
	topcmd9 = "sed -i -e 's/fillC12olam/"+str(C12olam)+"/' "+newtop
	
	mdpcmd01 = "sed -i -e 's/filltemperature/"+str(Tlam)+"/' "+newminmdp
	mdpcmd02 = "sed -i -e 's/fillpressure/"+str(Plam)+"/' "+newminmdp
	mdpcmd03 = "sed -i -e 's/fillnsteps/"+str(nstepseq)+"/' "+newminmdp
	mdpcmd04 = "sed -i -e 's/fillnstxout/"+str(nstxouteq)+"/' "+newminmdp
	mdpcmd05 = "sed -i -e 's/fillnstvout/"+str(nstvouteq)+"/' "+newminmdp
	mdpcmd06 = "sed -i -e 's/fillnstenergy/"+str(nstenergyeq)+"/' "+newminmdp
	mdpcmd19 = "sed -i -e 's/fillkappa/"+str(kappa)+"/' "+newminmdp
	
	mdpcmd1 = "sed -i -e 's/filltemperature/"+str(Tlam)+"/' "+neweqmdp
	mdpcmd2 = "sed -i -e 's/fillpressure/"+str(Plam)+"/' "+neweqmdp
	mdpcmd3 = "sed -i -e 's/fillnsteps/"+str(nstepseq)+"/' "+neweqmdp
	mdpcmd4 = "sed -i -e 's/fillnstxout/"+str(nstxouteq)+"/' "+neweqmdp
	mdpcmd5 = "sed -i -e 's/fillnstvout/"+str(nstvouteq)+"/' "+neweqmdp
	mdpcmd6 = "sed -i -e 's/fillnstenergy/"+str(nstenergyeq)+"/' "+neweqmdp
	mdpcmd20 = "sed -i -e 's/fillkappa/"+str(kappa)+"/' "+neweqmdp
	
	mdpcmd7 = "sed -i -e 's/filltemperature/"+str(Tlam)+"/' "+newmdp
	mdpcmd8 = "sed -i -e 's/fillpressure/"+str(Plam)+"/' "+newmdp
	mdpcmd9 = "sed -i -e 's/fillnsteps/"+str(nsteps)+"/' "+newmdp
	mdpcmd10 = "sed -i -e 's/fillnstxout/"+str(nstxout)+"/' "+newmdp
	mdpcmd11 = "sed -i -e 's/fillnstvout/"+str(nstvout)+"/' "+newmdp
	mdpcmd12 = "sed -i -e 's/fillnstenergy/"+str(nstenergy)+"/' "+newmdp
	mdpcmd21 = "sed -i -e 's/fillkappa/"+str(kappa)+"/' "+newmdp

	mdpcmd13 = "sed -i -e 's/filltemperature/"+str(Tlam)+"/' "+newenemdp
	mdpcmd14 = "sed -i -e 's/fillpressure/"+str(Plam)+"/' "+newenemdp
	mdpcmd15 = "sed -i -e 's/fillnsteps/"+str(0)+"/' "+newenemdp
	mdpcmd16 = "sed -i -e 's/fillnstxout/"+str(0)+"/' "+newenemdp
	mdpcmd17 = "sed -i -e 's/fillnstvout/"+str(0)+"/' "+newenemdp
	mdpcmd18 = "sed -i -e 's/fillnstenergy/"+str(nstenergy)+"/' "+newenemdp
	mdpcmd22 = "sed -i -e 's/fillkappa/"+str(kappa)+"/' "+newenemdp
	
	getnewtop = topcmd1+";"+topcmd2+";"+topcmd3+";"+topcmd4+";"+topcmd5+";"+topcmd6+";"+topcmd7+";"+topcmd8+";"+topcmd9
	getnewmdp =mdpcmd01+";"+mdpcmd02+";"+mdpcmd03+";"+mdpcmd04+";"+mdpcmd05+";"+mdpcmd06+";"+mdpcmd1+";"+mdpcmd2+";"+mdpcmd3+";"+mdpcmd4+";"+mdpcmd5+";"+mdpcmd6+";"+mdpcmd7+";"+mdpcmd8+";"+mdpcmd9+";"+mdpcmd10+";"+mdpcmd11+";"+mdpcmd12+";"+mdpcmd13+";"+mdpcmd14+";"+mdpcmd15+";"+mdpcmd16+";"+mdpcmd17+";"+mdpcmd18+";"+mdpcmd19+";"+mdpcmd20+";"+mdpcmd21+";"+mdpcmd22
	
	#os.system(topcmd1)
	os.system(getnewtop)
	os.system(getnewmdp)
	
	
	tprminfile = tprdir+"prodmin"+str(int(i))+".tpr"
	tpreqfile = tprdir+"prodeq"+str(int(i))+".tpr"
	tprfile = tprdir+"prod"+str(int(i))+".tpr"
	tprrerunfile = tprdir+"prodrerun"+str(int(i))+".tpr"
	
	prodmintrr   = trrdir+"prodminU"+str(int(i))+"X"+str(int(i))+".trr"
	prodmingro   = grodir+"prodminU"+str(int(i))+"X"+str(int(i))+".gro"
	prodminedr   = enedir+"prodminU"+str(int(i))+"X"+str(int(i))+".edr"
	prodminlog   = logdir+"prodminU"+str(int(i))+"X"+str(int(i))+".log"
	prodmindgdl  = dgdldir+"prodminU"+str(int(i))+"X"+str(int(i))+".xvg"
	prodmincpt   = cptdir+"prodminU"+str(int(i))+"X"+str(int(i))+".cpt"
	
	prodeqtrr   = trrdir+"prodeqU"+str(int(i))+"X"+str(int(i))+".trr"
	prodeqgro   = grodir+"prodeqU"+str(int(i))+"X"+str(int(i))+".gro"
	prodeqedr   = enedir+"prodeqU"+str(int(i))+"X"+str(int(i))+".edr"
	prodeqlog   = logdir+"prodeqU"+str(int(i))+"X"+str(int(i))+".log"
	prodeqdgdl  = dgdldir+"prodeqU"+str(int(i))+"X"+str(int(i))+".xvg"
	prodeqcpt   = cptdir+"prodeqU"+str(int(i))+"X"+str(int(i))+".cpt"
	
	prodtrr   = trrdir+"prodU"+str(int(i))+"X"+str(int(i))+".trr"
	prodgro   = grodir+"prodU"+str(int(i))+"X"+str(int(i))+".gro"
	prodedr   = enedir+"prodU"+str(int(i))+"X"+str(int(i))+".edr"
	prodlog   = logdir+"prodU"+str(int(i))+"X"+str(int(i))+".log"
	proddgdl  = dgdldir+"prodU"+str(int(i))+"X"+str(int(i))+".xvg"
	prodcpt   = cptdir+"prodU"+str(int(i))+"X"+str(int(i))+".cpt"
	
	if (precision == 'double'):
		# Equilibration
		geneqtpr =  gloc+"grompp_d -f "+neweqmdp+" -po "+mdpdir+"mdeqout.mdp"+" -c "+inputgro+" -o "+tpreqfile+" -p "+newtop+" -maxwarn 4"
		runeqmd =  gloc+"mdrun_d -pin off"+" -nt "+str(numproc)+" -s "+tpreqfile+" -o "+prodeqtrr+" -c "+prodeqgro+" -g "+prodeqlog+" -e "+prodeqedr+" -cpo "+prodeqcpt
		
		# Production run
		gentpr =  gloc+"grompp_d -f "+newmdp+" -po "+mdpdir+"mdout.mdp"+" -c "+prodeqgro+" -o "+tprfile+" -p "+newtop+" -maxwarn 4"
		runmd =  gloc+"mdrun_d -pin off"+" -nt "+str(numproc)+" -s "+tprfile+" -o "+prodtrr+" -c "+prodgro+" -g "+prodlog+" -e "+prodedr+" -cpo "+prodcpt
		
		# Tpr for rerun
		genreruntpr =  gloc+"grompp_d -f "+newenemdp+" -po "+mdpdir+"mdrerunout.mdp"+" -c "+prodgro+" -o "+tprrerunfile+" -p "+newtop+" -maxwarn 4"
	else:
		# Minimization
		genmintpr =  gloc+"grompp -f "+newminmdp+" -po "+mdpdir+"mdminout.mdp"+" -c "+inputgro+" -o "+tprminfile+" -p "+newtop+" -maxwarn 4"
		runminmd =  gloc+"mdrun -pin off"+" -nt "+str(numproc)+" -s "+tprminfile+" -o "+prodmintrr+" -c "+prodmingro+" -g "+prodminlog+" -e "+prodminedr+" -cpo "+prodmincpt
	
		# Equilibration
		#geneqtpr =  gloc+"grompp -f "+neweqmdp+" -po "+mdpdir+"mdeqout.mdp"+" -c "+inputgro+" -o "+tpreqfile+" -p "+newtop+" -maxwarn 4"
		geneqtpr =  gloc+"grompp -f "+neweqmdp+" -po "+mdpdir+"mdeqout.mdp"+" -c "+prodmingro+" -o "+tpreqfile+" -p "+newtop+" -maxwarn 4"
		runeqmd =  gloc+"mdrun -pin off"+" -nt "+str(numproc)+" -s "+tpreqfile+" -o "+prodeqtrr+" -c "+prodeqgro+" -g "+prodeqlog+" -e "+prodeqedr+" -cpo "+prodeqcpt
		
		# Production run
		gentpr =  gloc+"grompp -f "+newmdp+" -po "+mdpdir+"mdout.mdp"+" -c "+prodeqgro+" -o "+tprfile+" -p "+newtop+" -maxwarn 4"
		runmd =  gloc+"mdrun -pin off "+" -nt "+str(numproc)+" -s "+tprfile+" -o "+prodtrr+" -c "+prodgro+" -g "+prodlog+" -e "+prodedr+" -cpo "+prodcpt
		
		# Tpr for rerun
		genreruntpr =  gloc+"grompp -f "+newenemdp+" -po "+mdpdir+"mdrerunout.mdp"+" -c "+prodgro+" -o "+tprrerunfile+" -p "+newtop+" -maxwarn 4"
	
	combinegmacs = genmintpr+";"+runminmd+";"+geneqtpr+";"+runeqmd+";"+gentpr+";"+runmd+";"+genreruntpr 
	
	#combinegmacs = geneqtpr
        
	if (cluster):
		#pdb.set_trace()
	        newjob = scriptdir+case+"U"+str(int(i))+"X"+str(int(i))+".sh"
		copyshell = "cp "+baseshell+" "+newjob
		os.system(copyshell)
		open(newjob,"a+").write(combinegmacs)
		submitshell = "qsub "+newjob
		os.system(submitshell)
	else:	
		os.system(combinegmacs)	

        #rmcommand = " rm "+trrdir+"/*#*"+"; rm "+grodir+"/*#*"+"; rm "+enedir+"/*#*"+"; rm "+logdir+"/*#*"+";  rm "+dgdldir+"/*#*"+"; rm "+cptdir+"/*#*"+"; rm "+mdpdir+"/*#*"
		
	#os.system(rmcommand)	
	#pdb.set_trace()
	
