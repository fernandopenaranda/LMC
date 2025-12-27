#!/bin/bash
#SBATCH --job-name=array-job     # create a short name for your job
#SBATCH --output=slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across allnodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 ifmulti-threaded tasks)
#SBATCH --mem-per-cpu=8G         # memory per cpu-core (4G is default)
#SBATCH --time=23:00:00          # total run time limit (HH:MM:SS)
#SBATCH --array=0-99             # job array with index values 0, 1, 2, 3, 4
#SBATCH --error=slurm-%A.%a.err
#SBATCH --mail-user=fernando.penaranda@dipc.org

pathtofile=$1
Ezmin=$2
Ezmax=$3
evals=$4
N=$5
eta=$6

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Array length: $SLURM_ARRAY_TASK_MAX"

julia --compiled-modules=no $pathtofile $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_ARRAY_JOB_ID $Ezmin $Ezmax $evals $N $eta