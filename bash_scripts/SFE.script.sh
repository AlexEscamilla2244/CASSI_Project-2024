#!/bin/bash                                                                                                                                                                         
#SBATCH -J YSO_on_snaps -p development -N 1 --ntasks-per-node 56 -t 2:00:00 -A AST21002 


DIR1="/work2/10071/alexescamilla2244/frontera/CASSI_Project-2024/output/YSOobjects"
DIR2="/scratch3/03532/mgrudic/STARFORGE_RT/STARFORGE_v1.2/M2e4_R10_Z1_S0_A2_B0.1_I1_Res271_n2_sol0.5_42/output/dustemission/"

PYTHON_SCRIPT="/work2/10071/alexescamilla2244/frontera/CASSI_Project-2024/scripts/e_ff.py"


FILES1=($(ls -v -1 $DIR1/*.hdf5)) 
FILES2=($(ls -v -1 $DIR2/*.hdf5 | sed -n '64,$p'))

for i in ${!FILES1[@]}; do

   yso=${FILES1[$i]}
   dust=${FILES2[$i]}
   python $PYTHON_SCRIPT $yso $dust 
done
