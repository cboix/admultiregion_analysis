#!/bin/bash
# ---------------------------------------------------------
# Run differential celltype analysis for multi-region data:
# Updated: 10/14/2020
# OLD version of analysis
# ---------------------------------------------------------
source $HOME/data/DEVTRAJ/bin/mosaicism/config_mosaic.sh

conda activate mv_env

# Use the correct R version:
RCMD="/broad/compbio/cboix/software/miniconda3/envs/mv_env/bin/R"

ctlist=(Ast Mic_Immune Oli Opc Vasc_Epithelia Exc Inh)
ctlist=(Exc Inh Ast Opc Oli Vasc_Epithelia Mic_Immune)
for celltype in ${ctlist[@]}; do
    echo ${celltype} 
    cmd="$RCMD --slave -f $BINDIR/multiRegion/plot_differential_04_filter_low_lmm.R --args $celltype"
    echo "$cmd"
    if [[ "$celltype" == "Exc" ]]; then
        for ind in `seq 1 8`; do
            qsub -cwd -l h_vmem=80G -l h_rt=12:00:00 -N deg_multireg_$celltype -j y -b y -V -r y -o $DBDIR/out/ "$cmd $ind"
        done
    elif [[ "$celltype" == "Inh" ]]; then
        for ind in `seq 1 6`; do
            qsub -cwd -l h_vmem=50G -l h_rt=12:00:00 -N deg_multireg_$celltype -j y -b y -V -r y -o $DBDIR/out/ "$cmd $ind"
        done
    else
        qsub -cwd -l h_vmem=50G -l h_rt=12:00:00 -N deg_multireg_$celltype -j y -b y -V -r y -o $DBDIR/out/ "$cmd"
    fi
done


# Run the regression on fractions analysis:
RCMD="/broad/compbio/cboix/software/miniconda3/envs/mv_env/bin/R"
cmd="$RCMD --slave -f $BINDIR/multiRegion/calculate_metadata_01b_fractions_logreg_glmer.R"
qsub -cwd -t 913-1140 -l h_vmem=5G -l h_rt=0:45:00 -N frac_multireg_glmm -j y -b y -V -r y -o $DBDIR/out/ "$cmd"
# TODO: Add the nft_regional, etc.



# ------------------------------------------
# Run MAST RE models from the mtx directory:
# ------------------------------------------
conda activate r4-base

MTXDIR=$DBDIR/multiRegion/matrices/mtx
cd $MTXDIR

# for file in `ls all_brain_regions_filt_preprocessed_scanpy*.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Vasc_Epithelia*.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Mic_Immune*Mic*.rda`; do 
for file in `ls all_brain_regions_filt_preprocessed_scanpy*Oli*.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Exc*.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*.rda`; do 
    echo $file
    for ind in `seq 1 15`; do 
        # Run on SLURM:
        sbatch --time=24:00:00 $BINDIR/multiRegion/run_differential_ctbyregion_control.sh $file $ind 1250 $path
    done
done

MTXDIR=$DBDIR/multiRegion/matrices/mtx
cd $MTXDIR

# for file in `ls all_brain_regions_filt_preprocessed_scanpy*.EC.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Mic*.EC.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Ast*.EC.rda`; do 
#     echo $file
#     for ind in `seq 1 30`; do 
#         # Run on SLURM:
#         sbatch $BINDIR/multiRegion/run_differential_allregions_control.sh $file $ind 633 $path
#     done
# done


# Further parallelize....
# pathlist=(nft plaq_d plaq_n cogdxad nrad)
pathlist=(nft plaq_d plaq_n nrad)
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Ast*Ast*.EC.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Mic*.EC.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Opc*.EC.rda`; do 
# for file in `ls all_brain_regions_filt_preprocessed_scanpy*Vasc_Epithelia*.EC.rda`; do 
for file in `ls all_brain_regions_filt_preprocessed_scanpy*Oli*.EC.rda`; do 
    echo $file
    for pth in ${pathlist[@]}; do 
        echo $pth
        for ind in `seq 1 120`; do 
            # Run on SLURM:
            sbatch -n 16 $BINDIR/multiRegion/run_differential_allregions_control.sh $file $ind 160 $pth
        done
    done
done



