#! /usr/bin/env Rscript

### Part 0: Load Packages
library(monocle3)
library(Matrix)
library(gplots)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dplyr)
options(stringsAsFactors=F)

cds = readRDS("all_cells_labeled.rds")
stele_cds = readRDS("stele_cells_labeled.rds")
pericycle_cds = readRDS("pericycle_cells_labeled.rds")
ce_cds = readRDS("ce_cells_labeled.rds")

# Figure 1 
plot_cells(cds, color_cells_by="experiment", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_manual(values= c("00 hours"="cyan", "08 hours"="mediumpurple", "20 hours"="blue4"))
plot_cells(pericycle_cds, color_cells_by="experiment", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_manual(values= c("00 hours"="cyan", "08 hours"="mediumpurple", "20 hours"="blue4"))
plot_cells(pericycle_cds, genes="AT1G74560", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))
plot_cells(pericycle_cds, genes="AT1G22530", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))

# Figure 2 - Figue 4
norm_mat <- t(t(counts(pericycle_cds))/size_factors(pericycle_cds))
valid_genes <- read.table("validation_genes.tsv"); valid_genes <- valid_genes[order(valid_genes$V1), ]; valid_genes <- valid_genes[!(valid_genes$V2 %in% c("NAP", "MED21")), ]
valid_genes$V3 <- factor(valid_genes$V3, levels=c("Control", "Stemness", "Chromatin", "Cell_Cycle"))
valid_genes <- valid_genes[ match( c("ARF7", "ARF19", "AIL6", "TMO6", "BRXL1", "BRXL4", "OFP8", "LRP1", "NRP1", "NRP2", "RACK1A", "RACK1B", "RACK1C", "SMR6", "SMR11", "HDT4", "HD2B", "HDA3", "ORTH1"), valid_genes$V2 ), ]
validation_norm_mat <- norm_mat[ rownames(norm_mat) %in% valid_genes$V1, ]
xpp <- rowMeans(validation_norm_mat[ , colData(pericycle_cds)$celltype == "Xylem Pole Pericycle" ])
lrp <- rowMeans(validation_norm_mat[ , colData(pericycle_cds)$celltype == "Lateral Root Primordia" ])
mature <- rowMeans(validation_norm_mat[ , colData(pericycle_cds)$celltype == "Mature Pericycle" ])
ppp <- rowMeans(validation_norm_mat[ , colData(pericycle_cds)$celltype == "Phloem Pole Pericycle" ])
tmp <- cbind(xpp, lrp, mature, ppp); tmp <- tmp[match(valid_genes$V1, rownames(tmp)), ]; rownames(tmp) <- valid_genes$V2; colnames(tmp) <- c("Xylem Pole Pericycle", "Lateral Root Primordia", "Mature Pericycle", "Phloem Pole Pericycle")
tmp2 <- tmp[ match( valid_genes$V2, rownames(tmp) ), ];
pheatmap(tmp2[ rownames(tmp2) %in% valid_genes$V2[valid_genes$V3 %in% c("Control", "Chromatin")], ], angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F, gaps_row=c(2))
pheatmap(tmp2[ rownames(tmp2) %in% valid_genes$V2[valid_genes$V3 %in% c("Control", "Cell_Cycle")], ], angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F, gaps_row=c(2))
pheatmap(tmp2[ rownames(tmp2) %in% valid_genes$V2[valid_genes$V3 %in% c("Control", "Stemness")], ], angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F, gaps_row=c(2))

# Figure 5
fig_genes <- c("AT3G17185", "AT5G62000", "AT2G33860", "AT5G60450", "AT4G23980", "AT3G23030", "AT4G25490", "AT3G48100", "AT1G67710")
fig_genes_names <- c("TAS3", "ARF2", "ARF3", "ARF4", "ARF9", "IAA2", "CBF1", "ARR5", "ARR11")
sub_norm_mat <- norm_mat[ rownames(norm_mat) %in% fig_genes, ]
xpp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Xylem Pole Pericycle" ])
lrp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Lateral Root Primordia" ])
mature <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Mature Pericycle" ])
ppp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Phloem Pole Pericycle" ])
tmp <- cbind(xpp, lrp, mature, ppp); tmp <- tmp[match(fig_genes, rownames(tmp)), ]; rownames(tmp) <- fig_genes_names; colnames(tmp) <- c("Xylem Pole Pericycle", "Lateral Root Primordia", "Mature Pericycle", "Phloem Pole Pericycle")
pheatmap(tmp, angle_col=45, cellwidth=46.06, col=viridis(13), scale='row', cluster_cols=F, cluster_rows=F); pheatmap(tmp, angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F)
plot_cells(ce_cds, color_cells_by="experiment", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_manual(values= c("00 hours"="cyan", "08 hours"="mediumpurple", "20 hours"="blue4"))
plot_cells(ce_cds, genes="AT5G39610", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))
plot_cells(ce_cds, genes="AT4G21980", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))


### Supplemental 1
celltypes = unique(colData(cds)$Round1_Partition_labels)
col_int <- c("avg_hair_cell", "avg_cortex", "avg_non_hair", "avg_pericycle", "avg_xylem_pole_pericycle", "avg_columella", "avg_lateral_root", "avg_phloem_pole_pericycle", "avg_maturing_xylem", "avg_endodermis", "avg_phloem_cc")
marker_avg <- colData(cds)[ , colnames(colData(cds)) %in% col_int ]
mat <- matrix(0, length(celltypes), length(col_int)); colnames(mat) <- col_int; rownames(mat) <- celltypes
for (i in 1:length(celltypes)) {
    mat[ i, ] <- colMeans( as.matrix(marker_avg[ colData(cds)$Round1_Partition_labels == celltypes[i], ]) )
}
pheatmap(mat, angle_col=45, cellwidth=46.06, scale='column', col=viridis(13)); pheatmap(mat, angle_col=45, cellwidth=46.06, scale='column');

celltypes = unique(colData(stele_cds)$celltype)
col_int <- c("avg_pericycle", "avg_xylem_pole_pericycle", "avg_lateral_root", "avg_phloem_pole_pericycle", "avg_maturing_xylem", "avg_phloem_cc")
marker_avg <- colData(stele_cds)[ colData(stele_cds)$celltype %in% celltypes, colnames(colData(stele_cds)) %in% col_int ]
mat <- matrix(0, length(celltypes), length(col_int)); colnames(mat) <- col_int; rownames(mat) <- celltypes
for (i in 1:length(celltypes)) {
    mat[ i, ] <- colMeans( as.matrix( marker_avg[ rownames(marker_avg) %in% colnames(stele_cds)[ colData(stele_cds)$celltype == celltypes[i] ] , ] ) )
}
pheatmap(mat[c(1, 2, 3, 6, 7, 5, 4), c(5, 6, 2, 3, 1, 4)], angle_col=45, cellwidth=46.06, scale='column', cluster_rows=F, cluster_cols=F); 
plot_cells(stele_cds, color_cells_by="experiment", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_manual(values= c("00 hours"="cyan", "08 hours"="mediumpurple", "20 hours"="blue4"))

### Supplemental 2 (See GLM and MWW script for pseudotime)

## heatmap
norm_mat <- t(t(counts(pericycle_cds))/size_factors(pericycle_cds))
fig_genes <- c("AT2G42430", "AT3G58190", "AT3G11260", "AT4G23750", "AT5G18560", "AT1G75640", "AT5G56580", "AT1G34110")
fig_gene_names <- c("LBD16", "LBD29", "WOX5", "CRF2", "PUCHI", "MUS", "MKK6", "RGI5")
sub_norm_mat <- norm_mat[ rownames(norm_mat) %in% fig_genes, ]
xpp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Xylem Pole Pericycle" ])
lrp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Lateral Root Primordia" ])
mature <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Mature Pericycle" ])
ppp <- rowMeans(sub_norm_mat[ , colData(pericycle_cds)$celltype == "Phloem Pole Pericycle" ])
tmp <- cbind(xpp, lrp, mature, ppp); tmp <- tmp[match(fig_genes, rownames(tmp)), ]; rownames(tmp) <- fig_genes_names; colnames(tmp) <- c("Xylem Pole Pericycle", "Lateral Root Primordia", "Mature Pericycle", "Phloem Pole Pericycle")
pheatmap(tmp, angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F); pheatmap(tmp, angle_col=45, cellwidth=46.06, scale='row', cluster_cols=F, cluster_rows=F)

### Supplemental Figure 11
celltypes = unique(colData(ce_cds)$celltype)
col_int <- c("avg_cortex", "avg_endodermis")
marker_avg <- colData(ce_cds)[ , colnames(colData(ce_cds)) %in% col_int ]
mat <- matrix(0, length(celltypes), length(col_int)); colnames(mat) <- col_int; rownames(mat) <- celltypes
for (i in 1:length(celltypes)) {
    mat[ i, ] <- colMeans( as.matrix(marker_avg[ colData(ce_cds)$celltype == celltypes[i], ]) )
}
rownames(mat) <- c("Lateral Root Endodermis", "Endodermis", "Cortex")
pheatmap(mat, angle_col=45, cellwidth=46.06, scale='column', col=viridis(13)); pheatmap(mat, angle_col=45, cellwidth=46.06, scale='column');


### Supplement 13 (See GLM and MWW Script for pseudotime)
plot_cells(ce_cds, genes="AT1G72490", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))
plot_cells(ce_cds, genes="AT5G13080", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))
plot_cells(ce_cds, genes="AT2G17500", show_trajectory_graph=F, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F, alpha=0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))

