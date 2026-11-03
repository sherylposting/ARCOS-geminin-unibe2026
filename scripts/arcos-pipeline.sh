#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=00:10:00
#SBATCH -p all
#SBATCH --array=0-3
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --job-name=arcos-pipeline_run3

# arcos-pipeline.sh: one array task per well/site; steps run sequentially:
# segment -> track -> ce-arcos -> plots

set -e  # stop this task if any step fails

# Global variables - remember to edit
WORK_DIR=/home/slin/ARCOS-geminin/
TIFFPATH="/mnt/imaging.data/pgagliardi/MCF10A_TimeLapse_Geminin-Drugs/2019-05-03_Exp3_Azd53-Bez-Azd85/TIFFs/"
# remember to edit the SBATCH --array appropriately to the total wells * sites you specify
# ex. --array=0-2 for three wells * 1 site 
WELLS=(A3 B2 B3 B4)
SITES=(1)

cd $WORK_DIR

# Map array index to well/site
WELL=${WELLS[$((SLURM_ARRAY_TASK_ID / ${#SITES[@]}))]}
SITE=${SITES[$((SLURM_ARRAY_TASK_ID % ${#SITES[@]}))]}
ARGS="--well $WELL --site $SITE --tiffpath $TIFFPATH"

echo "Filepath: $TIFFPATH"
echo "Task $SLURM_ARRAY_TASK_ID: well=$WELL site=$SITE"

# Adjust to your conda install
source ~/miniconda3/etc/profile.d/conda.sh
conda activate arcos_tutorial

#echo "=== segment === started at: $(date)";  python segment.py  $ARGS
#echo "=== tracking === started at: $(date)";    python tracking.py    $ARGS
#echo "=== ce-arcos === started at: $(date)"; python ce-arcos.py $ARGS
echo "=== napari-save === started at: $(date)";    python napari-save.py    $ARGS
#echo "=== plots === started at: $(date)";    python plots.py    $ARGS

echo "arcos-pipeline finished at: $(date) for well=$WELL site=$SITE"