#!/bin/bash
#SBATCH -J nb_res
#SBATCH -N 1 
#SBATCH -n 32
#SBATCH -p kellis 
#SBATCH --time=24:00:00 
#SBATCH --mem=90G 
#SBATCH --output=/home/cboix/data/DEVTRAJ/db/out/nebula_res_%j.log 
#SBATCH --export=ALL 
#SBATCH --requeue 
#SBATCH --array=1-668

# Run as:
# sbatch --time=72:00:00 --array=2-139 $SBINDIR/differential_expression/calculate_degs/run_differential_01_NB_wRUV_selected_genes.sh

start=`date +%s`

# Cluster type:
if [[ "$SLURM_JOB_ID" == "" ]]; then
    TASK=$SGE_TASK_ID 
    JOBTYPE=SGE
else 
    TASK=$SLURM_ARRAY_TASK_ID
    JOBTYPE=SLURM
fi

pconf -t DEVTRAJ -s multiRegion -g b37

# Create this file concatenating runs from with setup_differential_runs_selected_genes.R
INFOFILE=$DBDIR/multiRegion/nebula_wRUV_selgenes_runlist.tsv
echo "[STATUS] Running task: $TASK"

# Needs r4-base running:
conda activate r4-base 

MTXDIR=$DBDIR/multiRegion/matrices/mtx
cd $MTXDIR
RCMD=$SFTDIR/miniconda3/envs/r4-base/bin/R
SCRIPTDIR=$BINDIR/multiRegion/differential_expression/calculate_degs/

# Run nebula with RUV-seq:
cmd="$RCMD --slave -f ${SCRIPTDIR}/run_differential_01_NB_wRUV_selected_genes.R --args $INFOFILE $TASK"

echo "$cmd"
bash -c "$cmd"

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
