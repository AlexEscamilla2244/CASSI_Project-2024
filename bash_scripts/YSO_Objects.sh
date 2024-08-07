#!/bin/bash                                                                     
#SBATCH -J YSO_on_snaps -p development -N 1 --ntasks-per-node 56 -t 2:00:00 -A AST21002

DIRECTORY="/scratch3/03532/mgrudic/STARFORGE_RT/STARFORGE_v1.2/M2e4_R10_Z1_S0_A2_B0.1_I1_Res271_n2_sol0.5_42/output/"

YSO="./CASSI_Research/starforge_tools/mockobs/yso.py"

exit="/work2/10071/alexescamilla2244/frontera/output"

for file in "$DIRECTORY"/snapshot*.hdf5; do
    if [ -f "$file" ]; then
        echo "Checking $file for field"
        
        # Check if the file contains the specific field
        if h5ls -r "$file" | grep -q "/PartType5/Coordinates"; then
            echo "Processing $file"
            python $YSO "$file"
        else
            echo "Field not found in $file, skipping."
        fi
    fi
done
