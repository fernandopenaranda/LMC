#!/bin/bash
#SBATCH --job-name=array-job     # create a short name for your job
#SBATCH --output=slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across allnodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 ifmulti-threaded tasks)
#SBATCH --mem-per-cpu=8G         # memory per cpu-core (4G is default)
#SBATCH --time=23:00:00          # total run time limit (HH:MM:SS)
#SBATCH --array=0-499             # job array with index values 0, 1, 2, 3, 4
#SBATCH --error=slurm-%A.%a.err
#SBATCH --mail-user=fernando.penaranda@dipc.org

pathtofile=$1
pdPID=$2
evals=$3
T=$4
tau=$5
which_observable=$6

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Array length: $SLURM_ARRAY_TASK_MAX"
echo "Which observable: $which_observable"

CMD="/scratch/ferpe/julia-1.9.4/bin/julia --compiled-modules=no $pathtofile $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_ARRAY_JOB_ID $pdPID $evals $T $tau $which_observable"
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
    printf "JOBID=%s CMD=%s\n" "$SLURM_ARRAY_JOB_ID" "$CMD" >> "julia_command_${which_observable}.txt"
fi
exec $CMD