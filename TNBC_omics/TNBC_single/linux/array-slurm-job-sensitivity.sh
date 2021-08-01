#!/bin/bash -l

#SBATCH
#SBATCH --job-name=chang_sbml_marcc_test
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
# number of cpus (threads) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=szhan121@jhu.edu

#### load and unload modules you may need
# module unload openmpi/intel
# module load mvapich2/gcc/64/2.0b
# module list

#### process job list data
filename=exp_roster_sensitivity_param.dat
declare -a myArray
mapfile -t myArray < $filename
nrJobs=${#myArray[@]}

#IS_MARCC_JOB=false
IS_MARCC_JOB=true

BIN="./TNBC_s_sim"

TSLICE="-t 800

BRIEF="-B"
#BRIEF=""

STAT="-S --stats-interval 1"

GRID=""
#GRID="-G 3 --grid-interval 40"

#### execute code 
time $BIN $TSLICE ${myArray[$((SLURM_ARRAY_TASK_ID-1))]} $BRIEF $STAT $GRID
echo "Finished with job $SLURM_ARRAY_TASK_ID"
#### mpiexec by default launches number of tasks requested