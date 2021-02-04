#!/usr/bin/env Rscript

### Part 0: Load Packages
library(monocle3)
library(Matrix)
library(gplots)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dplyr)
options(stringsAsFactors=F)

### Part 1: Load Files
## data files available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158761
## Note: unzip the files before running script!!

# path to matrix.mtx, barcodes.tsv, and genes.tsv
EXPR_MTX_PATH = "GSE158761_matrix.mtx"
CELL_ANN_PATH = "GSE158761_barcodes.tsv"
GENE_ANN_PATH = "GSE158761_genes.tsv"

# path to tsv files used for gene list used for stele and cortex & endodermis preprocessing
STELE_PREPROCESSING_PATH = "stele_markers_used.tsv"
CE_PREPROCESSING_PATH = "ce_markers_used.tsv" 

# load data
EXPR_MTX = readMM(EXPR_MTX_PATH)
# label columns and rows of the barcodes table
CELL_ANN = read.table(CELL_ANN_PATH, sep="\t");  colnames(CELL_ANN) = c("barcodes", "experiment", "celltype"); rownames(CELL_ANN) = CELL_ANN$barcodes
# extract sample number from barcodes
CELL_ANN$sample_number = as.numeric(unlist(strsplit(CELL_ANN$barcodes, "-"))[c(F, T)])
# label columns and rows of the genes table
GENE_ANN = read.table(GENE_ANN_PATH); colnames(GENE_ANN) = c("id", "gene_short_name"); rownames(GENE_ANN) = GENE_ANN$id

# load gene list
stele_preprocessing_genes = read.table(STELE_PREPROCESSING_PATH); stele_preprocessing_genes = stele_preprocessing_genes[, 1]
ce_preprocessing_genes = read.table(CE_PREPROCESSING_PATH); ce_preprocessing_genes = ce_preprocessing_genes[ , 1]

### Part 2: All cells UMAP
# make new cell data set object
# Note: the input data has already been filtered to remove genes that did have at least 10 UMIs or appear in at least 10 cells
# as well as genes upregulated by protoplasting
cds = new_cell_data_set(EXPR_MTX, cell_metadata=CELL_ANN, gene_metadata=GENE_ANN)
# preprocess and run PCA 
cds = preprocess_cds(cds, num_dim=100, method="PCA", norm_method="log", scaling=T, residual_model_formula_str="~ sample_number", verbose=T)
# reduce dimension to 2 dimensions via UMAP
cds = reduce_dimension(cds, reduction_method="UMAP", preprocess_method="PCA", umap.metric="cosine", umap.min_dist=0.1, umap.n_neighbors=15L, umap.nn_method="annoy", umap.fast_sgd=F, verbose=T)
# call cell clusters
cds = cluster_cells(cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=c(10^seq(-6,0)), verbose=T)
# draw trajectory graphs for each major partition
cds = learn_graph(cds, use_partition=T, close_loop=F, learn_graph_control=list(prune_graph=T), verbose=T)
# transfer inital round of clustering, and partition information to colData dataframe
colData(cds)$Round1_clusters = cds@clusters[["UMAP"]][["clusters"]]
colData(cds)$Round1_partition = cds@clusters[["UMAP"]][["partitions"]]

## calculate average expression of marker genes
# parse/ define marker genes for each marker line
markers_qc = read.table('markerLinesEnrichedGenes/AGL42_QUIESCENT_CENTER.txt', header=T)
markers_phloem_and_cc = read.table('markerLinesEnrichedGenes/APL_PHLOEM+CC.txt', header=T)
markers_hair_cell = read.table('markerLinesEnrichedGenes/COBL9_HAIR_CELL.txt', header=T)
markers_cortex = read.table('markerLinesEnrichedGenes/CORTEX.txt', header=T)
markers_non_hair = read.table('markerLinesEnrichedGenes/GL2_NON_HAIR_CELL.txt', header=T)
markers_pericycle = read.table('markerLinesEnrichedGenes/J2661_PERICYCLE.txt', header=T)
markers_xylem_pole_pericycle = read.table('markerLinesEnrichedGenes/JO121_XYLEM_POLE_PERICYCLE.txt', header=T)
markers_lateral_root_cap = read.table('markerLinesEnrichedGenes/LRC_LATERAL_ROOT_CAP.txt', header=T)
markers_columella = read.table('markerLinesEnrichedGenes/PET111_COLUMELLA.txt', header=T)
markers_lateral_root = read.table('markerLinesEnrichedGenes/RM1000_LATERAL_ROOT.txt', header=T)
markers_phloem_pole_pericycle = read.table('markerLinesEnrichedGenes/S17_PHLOEM_POLE_PERICYCLE.txt', header=T)
markers_maturing_xylem = read.table('markerLinesEnrichedGenes/S18_MATURING_XYLEM.txt', header=T)
markers_protophloem = read.table('markerLinesEnrichedGenes/S32_PROTOPHLOEM.txt', header=T)
markers_developing_xylem = read.table('markerLinesEnrichedGenes/S4_DEVELOPING_XYLEM.txt', header=T)
markers_endodermis = read.table('markerLinesEnrichedGenes/SCR5_ENDODERMIS.txt', header=T)
markers_phloem_cc = read.table('markerLinesEnrichedGenes/SUC2_PHLOEMCC.txt', header=T)
# split genes seperated by ;
markers_qc =  unlist(lapply(as.character(markers_qc$locus), strsplit, ";"))
markers_phloem_and_cc = unlist(lapply(as.character(markers_phloem_and_cc$LOCUS), strsplit, ";"))
markers_hair_cell = unlist(lapply(as.character(markers_hair_cell$locus), strsplit, ";"))
markers_cortex = unlist(lapply(as.character(markers_cortex$LOCUS), strsplit, ";"))
markers_non_hair = unlist(lapply(as.character(markers_non_hair$locus), strsplit, ";"))
markers_pericycle = unlist(lapply(as.character(markers_pericycle$LOCUS), strsplit, ";"))
markers_xylem_pole_pericycle = unlist(lapply(as.character(markers_xylem_pole_pericycle$locus), strsplit, ";"))
markers_lateral_root_cap = unlist(lapply(as.character(markers_lateral_root_cap$locus), strsplit, ";"))
markers_columella = unlist(lapply(as.character(markers_columella$locus), strsplit, ";"))
markers_lateral_root = unlist(lapply(as.character(markers_lateral_root$locus), strsplit, ";"))
markers_phloem_pole_pericycle = unlist(lapply(as.character(markers_phloem_pole_pericycle$locus), strsplit, ";"))
markers_maturing_xylem = unlist(lapply(as.character(markers_maturing_xylem$locus), strsplit, ";"))
markers_protophloem = unlist(lapply(as.character(markers_protophloem$locus), strsplit, ";"))
markers_developing_xylem = unlist(lapply(as.character(markers_developing_xylem$locus), strsplit, ";"))
markers_endodermis = unlist(lapply(as.character(markers_endodermis$locus), strsplit, ";"))
markers_phloem_cc = unlist(lapply(as.character(markers_phloem_cc$locus), strsplit, ";"))
# parse Parizot Markers
Parizot_markers = read.table("markerLinesEnrichedGenes/Parizot_S2_Pericycle_Markers.txt", sep="\t", stringsAsFactors=F);
Parizot_markers = Parizot_markers[, 1:2]; colnames(Parizot_markers) = c("gene", "celltype")
# calculate average expression of marker genes
cds_mat = t(t(counts(cds))/size_factors(cds)) # Divide by Size Factor
colData(cds)$avg_qc = apply(cds_mat[rownames(cds_mat) %in% markers_qc, ], 2, mean)
colData(cds)$avg_phloem_and_cc = apply(cds_mat[rownames(cds_mat) %in% markers_phloem_and_cc, ], 2, mean)
colData(cds)$avg_hair_cell = apply(cds_mat[rownames(cds_mat) %in% markers_hair_cell, ], 2, mean)
colData(cds)$avg_cortex = apply(cds_mat[rownames(cds_mat) %in% markers_cortex, ], 2, mean)
colData(cds)$avg_non_hair = apply(cds_mat[rownames(cds_mat) %in% markers_non_hair, ], 2, mean)
colData(cds)$avg_pericycle = apply(cds_mat[rownames(cds_mat) %in% markers_pericycle, ], 2, mean)
colData(cds)$avg_xylem_pole_pericycle = apply(cds_mat[rownames(cds_mat) %in% markers_xylem_pole_pericycle, ], 2, mean)
colData(cds)$avg_lateral_root_cap = apply(cds_mat[rownames(cds_mat) %in% markers_lateral_root_cap, ], 2, mean)
colData(cds)$avg_columella = apply(cds_mat[rownames(cds_mat) %in% markers_columella, ], 2, mean)
colData(cds)$avg_lateral_root = apply(cds_mat[rownames(cds_mat) %in% markers_lateral_root, ], 2, mean)
colData(cds)$avg_phloem_pole_pericycle = apply(cds_mat[rownames(cds_mat) %in% markers_phloem_pole_pericycle, ], 2, mean)
colData(cds)$avg_maturing_xylem = apply(cds_mat[rownames(cds_mat) %in% markers_maturing_xylem, ], 2, mean)
colData(cds)$avg_protophloem = apply(cds_mat[rownames(cds_mat) %in% markers_protophloem, ], 2, mean)
colData(cds)$avg_developing_xylem = apply(cds_mat[rownames(cds_mat) %in% markers_developing_xylem, ], 2, mean)
colData(cds)$avg_endodermis = apply(cds_mat[rownames(cds_mat) %in% markers_endodermis, ], 2, mean)
colData(cds)$avg_phloem_cc = apply(cds_mat[rownames(cds_mat) %in% markers_phloem_cc, ], 2, mean)
colData(cds)$avg_Parizot_LRP = apply(cds_mat[rownames(cds_mat) %in% Parizot_markers$gene[Parizot_markers$celltype == "LRP (P)"], ], 2, mean)
colData(cds)$avg_Parizot_pericycle = apply(cds_mat[rownames(cds_mat) %in% Parizot_markers$gene[Parizot_markers$celltype == "Pericycle (P)"], ], 2, mean)
colData(cds)$avg_Parizot_ppp = apply(cds_mat[rownames(cds_mat) %in% Parizot_markers$gene[Parizot_markers$celltype == "Phloem Pole Pericycle (P)"], ], 2, mean)
colData(cds)$avg_Parizot_xpp = apply(cds_mat[rownames(cds_mat) %in% Parizot_markers$gene[Parizot_markers$celltype == "Xylem Pole Pericycle (P)"], ], 2, mean)

### Part 3: Stele Cells UMAP
# Note: since partition numbering and cluster number can vary, cells will be subsetted using celltype lalels 
non_stele_cells = colData(cds)$celltype == "Epidermis" | colData(cds)$celltype == "Columella/ Root Cap" | colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Endodermis" | colData(cds)$celltype == "Lateral Root Endodermis"
# subset data for stele cells
stele_cds = cds[, !non_stele_cells]
# redo preprocessing and PCA
stele_cds <- preprocess_cds(stele_cds, num_dim=100, method="PCA", norm_method="log", scaling=T, use_genes=stele_preprocessing_genes, residual_model_formula_str="~ sample_number", verbose=T)
# redo UMAP
stele_cds <- reduce_dimension(stele_cds, reduction_method="UMAP", preprocess_method="PCA", umap.metric="cosine", umap.min_dist=0.1, umap.n_neighbors=15L, umap.nn_method="annoy", umap.fast_sgd=F, verbose=T)
# identify new clusters in UMAP
stele_cds <- cluster_cells(stele_cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=0.001, verbose=T)
# draw trajectories ontop of UMAP
stele_cds <- learn_graph(stele_cds, use_partition=T, close_loop=F, learn_graph_control=list(prune_graph=T), verbose=T)
# transfer second round of clustering and partition information to colData dataframe
colData(stele_cds)$Round2_clusters <- stele_cds@clusters[["UMAP"]][["clusters"]]
colData(stele_cds)$Round2_partition <- stele_cds@clusters[["UMAP"]][["partitions"]]


### Part 4: XPP, LRP, MP, and PPP UMAP
non_pericycle_cells = colData(stele_cds)$celltype == "Phloem" | colData(stele_cds)$celltype == "Xylem" | colData(stele_cds)$celltype == "Ambiguous Stele Cells"
# subset data for pericycle cells
pericycle_cds = stele_cds[, !non_pericycle_cells]
# redo preprocessing w/ pericycle marker gene sets
pericycle_cds <- preprocess_cds(pericycle_cds, num_dim=100, method="PCA", norm_method="log", use_genes=stele_preprocessing_genes, scaling=T, residual_model_formula_str="~ sample_number", verbose=T)
# redo UMAP
pericycle_cds <- reduce_dimension(pericycle_cds, reduction_method="UMAP", preprocess_method="PCA", umap.metric="cosine", umap.min_dist=0.01, umap.n_neighbors=15L, umap.nn_method="annoy", umap.fast_sgd=F, verbose=T)
# recluster and redraw trajectories
pericycle_cds <- cluster_cells(pericycle_cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=0.0005, verbose=T)
pericycle_cds <- learn_graph(pericycle_cds, use_partition=T, close_loop=F, learn_graph_control=list(prune_graph=T), verbose=T)
pericycle_cds <- cluster_cells(pericycle_cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=0.001, verbose=T)
# transfer second round of clustering, and partition information to colData dataframe
colData(pericycle_cds)$Round3_clusters <- pericycle_cds@clusters[["UMAP"]][["clusters"]]
colData(pericycle_cds)$Round3_partition <- pericycle_cds@clusters[["UMAP"]][["partitions"]]

### Part 5: Cortex and Endodermis UMAP
cortex_n_endodermis_cells = colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Endodermis" | colData(cds)$celltype == "Lateral Root Endodermis"
# subset data for cortex and endodermis cells
ce_cds = cds[, cortex_n_endodermis_cells]
ce_cds <- preprocess_cds(ce_cds, num_dim=100, method="PCA", norm_method="log", scaling=T, residual_model_formula_str="~ sample_number", use_genes=ce_preprocessing_genes, verbose=T)
ce_cds <- reduce_dimension(ce_cds, reduction_method="UMAP", preprocess_method="PCA", umap.metric="cosine", umap.min_dist=0.1, umap.n_neighbors=15L, umap.nn_method="annoy", umap.fast_sgd=F, verbose=T)
ce_cds <- cluster_cells(ce_cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=0.001, verbose=T)
ce_cds <- learn_graph(ce_cds, use_partition=T, close_loop=F, learn_graph_control=list(prune_graph=T), verbose=T)
ce_cds <- cluster_cells(ce_cds, reduction="UMAP", k=20, louvain_iter=1, partition_qval=0.05, weight=T, resolution=c(10^seq(-6,0)), verbose=T)
colData(ce_cds)$Round2_clusters <- ce_cds@clusters[["UMAP"]][["clusters"]]
colData(ce_cds)$Round2_partition <- ce_cds@clusters[["UMAP"]][["partitions"]]

