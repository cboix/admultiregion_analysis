#!/bin/bash
#SBATCH -J mast_res
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH -p kellis 
#SBATCH --time=12:00:00 
#SBATCH --mem=20G 
#SBATCH --output=/home/cboix/data/DEVTRAJ/db/out/mast_res_%j.log 
#SBATCH --export=ALL 
#SBATCH --requeue 

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
echo "[STATUS] Queueing: $celltype $subtype $region $ind"

# Needs r4-base running:
conda activate r4-base 

source $HOME/data/DEVTRAJ/bin/mosaicism/config_mosaic.sh ALZ


MTXDIR=$DBDIR/multiRegion/matrices/mtx
cd $MTXDIR
RCMD=$SFTDIR/miniconda3/envs/r4-base/bin/R

# Run MAST DE with rand effects:
cmd="$RCMD --slave -f $BINDIR/multiRegion/run_differential_ctbyregion_batched.R --args $celltype $subtype $region $chunksize $ind $path"

echo "$cmd"
bash -c "$cmd"

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
