#!/bin/bash
#SBATCH -J psbulk_deseq2
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH -p kellis 
#SBATCH --time=3:00:00 
#SBATCH --mem=30G 
#SBATCH --output=/home/cboix/data/DEVTRAJ/db/out/psbulk_deseq2_%j.log 
#SBATCH --export=ALL 
#SBATCH --requeue 
#SBATCH --array=1-668

# Run as:
# sbatch --array=1-668 $SBINDIR/differential_expression/calculate_degs/run_differential_04_DESeq2.sh

start=`date +%s`

# Cluster type:
if [[ "$SLURM_JOB_ID" == "" ]]; then
    TASK=$SGE_TASK_ID 
    JOBTYPE=SGE
else 
    TASK=$SLURM_ARRAY_TASK_ID
    JOBTYPE=SLURM
fi

cmd="pconf -t DEVTRAJ -s multiRegion -g b37"
bash -c "$cmd"

# Create this file concatenating runs from with setup_differential_runs.R
# INFOFILE=$SDBDIR/nebula_wRUV_joint_runlist.tsv
INFOFILE=$SDBDIR/DEG_multiRegion_SI_ACE_runlist.tsv
echo "[STATUS] Running task: $TASK"

# Ok with r_env, (r4-base for nebula):
MCENV=r_env
conda activate $MCENV

RCMD=$SFTDIR/miniconda3/envs/$MCENV/bin/R
SCRIPTDIR=$SBINDIR/differential_expression/calculate_degs/

# Run DESeq2 on pseudobulk:
cmd="$RCMD --slave -f ${SCRIPTDIR}/run_differential_04_DESeq2.R --args $INFOFILE $TASK"
echo "$cmd"
bash -c "$cmd"

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
