#!/bin/bash
# Rerun.sh
# By Gaurav, 2023/6/12
hmDir=$PWD
for N in 1 20 50 100 150 200 250 300
	do
		mkdir N_${N}
		cp lammps_sim.in N_${N}
		cp submit_N.lsf N_${N} 
		cp prod.restart1 N_${N}
		cd N_${N}
		sed -i "s/NUMBER/${N}/g" "lammps_sim.in"
		sed -i "s/NUMBER/${N}/g" "submit_N.lsf"
		bsub < submit_N.lsf
		cd ${hmDir}
		echo $PWD
	done



