#!/bin/bash --login
#SBATCH --job-name=MYJOB
#SBATCH --partition=gpu-dev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=1
#SBATCH --time=5:00
#SBATCH --account=pawsey0184-gpu

export MPICH_GPU_SUPPORT_ENABLED=1 #This allows for GPU-aware MPI communication among GPUs
export ROCR_VISIBLE_DEVICES=0

echo "Node list: " $SLURM_JOB_NODELIST

#Loading needed modules (adapt this for your own purposes):
module load craype-accel-amd-gfx90a

module use /software/setonix/unsupported
module load rocm/6.3.1
module load python/3.11.6
module load py-numpy/1.26.4

source /software/projects/pawsey0185/paolo/python_venv/work/bin/activate

cmd="srun -N 1 -n 1 -c 1"

${cmd} runOpenMM md.yaml --log md.0.log

python ffmd.py > md.1.log
