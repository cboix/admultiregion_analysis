#!/bin/bash
#SBATCH -J mast_res_all
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH -p kellis 
#SBATCH --time=48:00:00 
#SBATCH --mem=20G 
#SBATCH --output=/home/cboix/data/DEVTRAJ/db/out/mast_res_all_%j.log 
#SBATCH --export=ALL 
#SBATCH --requeue 
# -----------------------------------------------------
# Runs of MAST + RE for DGE in all regions in a subtype
# Updated: 01/29/2021
# -----------------------------------------------------

# OLD version of analysis

# Arguments:
file=$1
ind=$2
chunksize=$3
path=$4

start=`date +%s`

base=${file##*majorcelltype.}
pref=${base%%.rda}
cts=${pref%.*}
celltype=${cts%.*}
subtype=${cts#*.}
region=${pref##*.}
echo "[STATUS] Queueing: $celltype $subtype allregions $ind"

# Needs r4-base running:
conda activate r4-base 

source $HOME/data/DEVTRAJ/bin/mosaicism/config_mosaic.sh ALZ

MTXDIR=$DBDIR/multiRegion/matrices/mtx
cd $MTXDIR
RCMD=$SFTDIR/miniconda3/envs/r4-base/bin/R

# Run MAST DE with rand effects across ALL regions:
cmd="$RCMD --slave -f $BINDIR/multiRegion/run_differential_allregions_batched.R --args $celltype $subtype $chunksize $ind $path"

echo "$cmd"
bash -c "$cmd"

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
