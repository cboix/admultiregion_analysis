#!/usr/bin/R
# ---------------------------------------------------------
# Run all of the modules resource creation steps, in order: 
# Updated 11/29/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')

# Scripts prefix and location:
modbindir = paste0(sbindir, 'modules/')
fnprefix = 'modules_degenr_'

# The set of scripts to run (in order):
fnlist = c('01_plot_overall.R',
    '02_plot_enrichments.R',
    '03_plot_pseudobulk_modules.R',
    '04_metadata_enrichments.R',
    '05_plot_scores_with_annotation.R',
    '06_plot_scoreumaps_with_enr.R',
    '07_replot_representations.R')


# List to make representative figure panels:
replist = c('01_plot_overall.R',
    '02_plot_enrichments.R',
    '03_plot_pseudobulk_modules.R',
    '04_metadata_enrichments.R',
    '13_representative_panel.R')

# The set of scripts to run (in order):
fnlist = c('01_plot_overall.R',
    # '02_plot_enrichments.R',
    # '03_plot_pseudobulk_modules.R',
    '04_metadata_enrichments.R',
    '05_plot_scores_with_annotation.R',
    '06_plot_scoreumaps_with_enr.R',
    '07_replot_representations.R')



# Runsets that we can run/have pre-processed results:
# ---------------------------------------------------
graph_id = 'boot'
runlist = c('Ast','Oli','Inh','Opc','Mic_Immune','Vasc_Epithelia',
    'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')
# runlist = c('ECneurons', 'THneurons', 'HCneurons', 'CTXneurons')
runlist = c('Ast', 'Mic_Immune', 'Vasc_Epithelia', 'Oli','Inh','Opc',
    'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons', 'Glial', 'All')


# Run through all of the analyses, one by one:
# --------------------------------------------
# runset = 'Vasc_Epithelia'
for (runset in runlist[1:length(runlist)]){
    for (fnfile in fnlist[1:length(fnlist)]){
    # for (fnfile in replist){
        print(paste("Will run:", fnfile, "with", runset, graph_id))

        # Run the script:
        commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id)}
        source(paste0(modbindir, fnprefix, fnfile))

        print(paste("Finished running:", fnfile))
    }
}

