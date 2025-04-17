#!/bin/bash
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --time=1:00:00
#SBATCH --nodes=8
#SBATCH --job-name=Chabaoum
#SBATCH --output=/network/iss/cenir/analyse/meeg/CARACAS/Test_Max/logs/log-%a
#SBATCH --error=/network/iss/cenir/analyse/meeg/CARACAS/Test_Max/logs/err/err-%a
#SBATCH --array=1-375

echo started on $HOSTNAME
date
module load matlab
tic="`date +%s`"
cmd=$( printf "script_02_ds003690_precomp(%d)" ${SLURM_ARRAY_TASK_ID})
matlab -nodesktop -nodisplay -nosoftwareopengl -r "cd /network/iss/cenir/analyse/meeg/CARACAS/Test_Max/code; $cmd;exit;"
let toc="`date +%s`"
let sec="$toc - $tic"
let heu="$sec / 3600"
let min="($sec - $heu * 3600) / 60"
echo Elapsed time: $heu H : $min min
