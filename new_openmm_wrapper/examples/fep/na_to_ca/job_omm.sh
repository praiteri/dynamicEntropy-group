#!/bin/bash --login
#SBATCH --job-name=MYJOB
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0185-gpu
##SBATCH --qos=high

export MPICH_GPU_SUPPORT_ENABLED=1 #This allows for GPU-aware MPI communication among GPUs
export OMP_NUM_THREADS=1           #This controls the real CPU-cores per task for the executable
export ROCR_VISIBLE_DEVICES=0

#Loading needed modules (adapt this for your own purposes):
module load craype-accel-amd-gfx90a

module use /software/setonix/unsupported
module load rocm/6.3.1
module load python/3.11.6
module load py-numpy/1.26.1

source /software/projects/pawsey0185/paolo/python_venv/work/bin/activate

# module load rocm/5.7.3
# module load singularity/3.11.4
# cntr=/software/projects/pawsey0184/paolo/setonix/containers/openmm8.sif
# cmd="srun singularity run -B ${PWD}:${HOME} ${cntr}"

n=`ls fep.*.yaml | wc -l`
for ((i=0;i<$n;i++));do
    runOpenMM_lite fep.${i}.yaml > md.${i}.log
done
