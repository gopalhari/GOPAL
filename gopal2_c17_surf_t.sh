#!/bin/sh
####PBS -l nodes=2:ppn=20
#PBS -l nodes=1:ppn=20
#PBS -l walltime=240:00:00
#PBS -e "$PBS_JOBID".err
#PBS -o "$PBS_JOBID".out

####PBS -W group_list=shirtsgroup
####PBS -q shirtscluster
####PBS -l walltime=120:00:00
####PBS -l nodes=1:ppn=4
###PBS -l mem=8gb
####PBS -m abe

set -x
export main=/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/
export GMXLIB=/home/hpaliwal/soft/gromacs462files/
export OMP_NUM_THREADS=1
export PYTHONPATH=$main/water_p_sweep/scripts/:$PYTHONPATH
export PATH=$main/water_p_sweep/scripts/:$PATH

cd $main/water_p_sweep/scripts/optscan/jobdir2/0of2/

echo "PBS job id is $PBS_JOBID"
echo "PBS nodefile is at $PBS_NODEFILE"
NPROCS=$(wc -l < "$PBS_NODEFILE")
echo "NPROCS is $NPROCS"
###cat "$PBS_NODEFILE" | uniq > unode
###
###mynodes=`cat $PBS_NODEFILE | uniq`
###
###
###
#### Kill Parallel Python
###for node in $mynodes
###do
###    ssh $node 'killall python'
###    ssh $node 'pkill python'
###    sleep 5	
###done
###
###
###
###for node in $mynodes
###do
###    ssh -f $node '/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python /home/hpaliwal/soft/python/Python-2.7.6/build/bin/ppserver.py -p 35000 -s "krsna" -r &'
###    ##ssh -f $node '/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python /home/hpaliwal/soft/python/Python-2.7.6/build/bin/ppserver.py -p 35000 -w 60 -s "krsna" -r &'
###    sleep 5
###done
#### execute the program
#/home/hpaliwal/PhD_backup/working_PhD_water_exercise_scripts/nomad.3.8.1/bin/nomad gop2_c16_surf_t

#/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final2.py gopal_f
/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc.py gopal_f
#/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py gopal_f
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py opc3
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py opc 
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip4p-ew
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip3p
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py spce
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip4p_2005
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip4p
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py spc
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip4pst
##/home/hpaliwal/soft/python/Python-2.7.6/build/bin/python krsna_t4_rc_final.py tip3pst
##
#### Kill Parallel Python
###for node in $mynodes
###do
###    ssh $node 'killall python'
###done
