#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH -p all
#SBATCH -o logs/%x_%j.out 
#SBATCH -e logs/%x_%j.err 
#SBATCH --job-name=arcos_run1

echo "Job name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"

WELL="B2"
SITE=4
TIFFPATH="/mnt/imaging.data/pgagliardi/MCF10A_TimeLapse_Geminin-Drugs/2019-05-03_Exp3_Azd53-Bez-Azd85/TIFFs/"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate arcos_tutorial

WORK_DIR=/home/slin/ARCOS-geminin/

cd $WORK_DIR

python ce-arcos.py --well $WELL --site $SITE --tiffpath $TIFFPATH

echo "Job ended at: $(date)"