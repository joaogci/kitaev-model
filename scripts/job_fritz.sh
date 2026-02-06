#!/bin/bash -l
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# Most variables are self-explanatory, one exceptions might be Nthreads, which is referring to the number of OpenMP threads per MPI task.
# In general we found that ALF does not profit from hyper-threading such that we suggest to only use physical cores.
#
#SBATCH --job-name=kitaev-hyper
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joao.carvalho-inacio@uni-wuerzburg.de
# Wall clock limit (HH:MM:SS):
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE

#available partitions: singlenode, multinode
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10

unset SLURM_EXPORT_ENV
module purge
module load gcc/14.2.0
module load gsl

# the follwing environment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be adapted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyper-threading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following environment variables are Intel specific.
#export KMP_AFFINITY=verbose,granularity=fine,compact
export KMP_AFFINITY=granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache
export I_MPI_PIN_ORDER=scatter

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1


# bash ./out_to_in.sh
srun /home/atuin/b133ae/b133ae21/kitaev_hyperbolic/kitaev-model/src/main /home/atuin/b133ae/b133ae21/kitaev_hyperbolic/kitaev-model

