#!/usr/bin/env bash
# Get modules for major cell types using scdemon:
# TODO: At some point, put up as an array job instead
# ---------------------------------------------------
pconf -t DEVTRAJ -s multiRegion -g b37
MBDIR=${SBINDIR}/modules

conda activate scanpy_env


# Overall runs, with defaults (5% cut, z=4.5)
# -------------------------------------------
# NOTE: Resolution 3 might be too much for some cell types
ctlist=(Ast Mic_Immune Vasc_Epithelia Oli Opc Inh Exc All Glial)
# ctlist=(Oli Opc Inh Exc All Glial)
for ct in ${ctlist[@]}; do
    # python $MBDIR/modules_01_params_bootstraps.py --ct ${ct}
    python $MBDIR/modules_02_basic_graphs.py --ct ${ct} --res 3
done

# Run neuronal subsets as well:
subsets=(ECneurons HCneurons THneurons CTXneurons)
for st in ${subsets[@]}; do
    # python $MBDIR/modules_01_params_bootstraps.py --ct Exc --st ${st}
    python $MBDIR/modules_02_basic_graphs.py --ct Exc --st ${st} --res 3
    # python $MBDIR/modules_01_params_bootstraps.py --ct Exc --st ${st} --filt 0.01 --z 4.75
    # python $MBDIR/modules_02_basic_graphs.py --ct Exc --st ${st} --filt 0.01 --z 4.75
done


# Inh and Exc with more genes, stricter cutoff:
# ---------------------------------------------
neulist=(Inh Exc)
for ct in ${neulist[@]}; do
    python $MBDIR/modules_01_params_bootstraps.py --ct ${ct} --filt 0.01 --z 4.75
    python $MBDIR/modules_02_basic_graphs.py --ct ${ct} --filt 0.01 --z 4.75
done


# Tarball all of the results:
# ---------------------------
MODDIR=${SDBDIR}/modules
cd $MODDIR
tar -cvzf stats.tar.gz enrstats_* allgene_assign* module_df_* allpath*

# Pack representation data:
MODDIR=${SDBDIR}/modules
cd $MODDIR
tar -cvzf repr.tar.gz graph_*tsv module_cell_scores*.tsv.gz

