# Auxiliary plotting settings:
library(ComplexHeatmap)
library(circlize)


# Global parameters for heatmaps:
# -------------------------------
ht_opt$heatmap_row_names_gp = gpar(fontsize = 5)
ht_opt$heatmap_column_names_gp = gpar(fontsize = 5)
ht_opt$heatmap_row_title_gp = gpar(fontsize = 5.5)
ht_opt$heatmap_column_title_gp = gpar(fontsize = 5.5)
ht_opt$legend_title_gp = gpar(fontsize = 5.5, font=2)
ht_opt$legend_labels_gp = gpar(fontsize = 5)
ht_opt$legend_grid_height = unit(2.5, 'mm')
ht_opt$legend_grid_width = unit(2.5, 'mm')
gridtxt.fs = 4.75

htsc = 2.5
ht_opt$DENDROGRAM_PADDING = unit(.5 / htsc, 'mm')
ht_opt$DIMNAME_PADDING = unit(1 / htsc, 'mm')
ht_opt$COLUMN_ANNO_PADDING = unit(1 / htsc, 'mm')
ht_opt$ROW_ANNO_PADDING = unit(1 / htsc, 'mm')
ht_opt$HEATMAP_LEGEND_PADDING = unit(2 / htsc, 'mm')
ht_opt$ANNOTATION_LEGEND_PADDING = unit(2 / htsc, 'mm')


# Functions:
# ----------
saveHeatmap = function(ht, pltprefix, w, h){
    require(ComplexHeatmap)
    # Save to PDF and to PNG:
    pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
    draw(ht, ht_gap=unit(0.5, 'mm'))
    dev.off()
    png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
    draw(ht, ht_gap=unit(0.5, 'mm'))
    dev.off()
    print(pltprefix)
}


saveGGplot = function(gp, pltprefix, w, h){
    require(ggplot2)
    ggsave(paste0(pltprefix, '.pdf'), gp, units='in', dpi=450, width=w, height=h)
    ggsave(paste0(pltprefix, '.png'), gp, units='in', dpi=450, width=w, height=h)
    print(pltprefix)
}

