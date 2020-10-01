#!/bin/bash -l

#SBATCH -J SpaFHy_Peat
#SBATCH --account=project_2003225
#SBATCH -o my_output_%j
#SBATCH -e my_output_err_%j
#SBATCH -t 00:15:00
#SBATCH --cpus-per-task=10
#SBATCH -p serial 
#SBATCH --partition=test
#SBATCH --mem-per-cpu=4G
#SBATCH -N 2
#SBATCH -n 2
#SBATCH --mail-type=END
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kyle.eyvindson@luke.fi

module load r-env
module load geoconda/3.7
#pip3 install contextily==1.0rc2


for i in {0..400..200}
	do
	srun -J "$i" -N 1 -n 1 python yasso_run_multiprocessing.py --f /scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/ --o /scratch/project_2003225/GIT/Yasso/Output_files_Yasso/ --s $i &
	echo $i
	
	done


wait