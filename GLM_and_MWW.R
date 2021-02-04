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


### GLM 
# load pericycle data
cds = readRDS("pericycle_cells_labeled.rds")

## XPP vs. Lateral Root Primordia
lm = fit_models(cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Lateral Root Primordia" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "stele_XPPvLRP_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv  
write.table(ot, "stele_XPPvLRP_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP vs. Mature Pericycle 
lm = fit_models(cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Mature Pericycle" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "stele_XPPvMAT_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_XPPvMATUREPERICYCLE_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP and Mature Pericycle pseudotime
cds_xpp_mat = cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Mature Pericycle" ]
cds_xpp_mat <- order_cells(cds_xpp_mat) # this is interactive!!!! Select XPP cells as root
colData(cds_xpp_mat)$pseudotime <- cds_xpp_mat@principal_graph_aux$UMAP$pseudotime # save pseudotime info
lm <- fit_models(cds_xpp_mat,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "pseudotime_lm_xpp_mat_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_XPPtoMAT_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP and Lateral Root pseudotime
cds_xpp_lrp = cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Lateral Root Primordia" ]
cds_xpp_lrp <- order_cells(cds_xpp_lrp) # this is interactive!!!! Select XPP cells as root
colData(cds_xpp_lrp)$pseudotime <- cds_xpp_lrp@principal_graph_aux$UMAP$pseudotime # save pseudotime info
lm <- fit_models(cds_xpp_lrp,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "pseudotime_lm_xpp_lrp_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_XPPtoLRP_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# plotting
xpp_tab = read.table("xpp_lrp_mat_pseudotime.txt", sep="\t", header=T) # slim version of supplemental table 1
plot_cells(cds, genes = xpp_tab, show_trajectory_graph=T, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F) + xlim(-4.5, 2) + ylim(-6, 0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))


## Ambiguous Stele Cells vs. all other stele cells
# load stele data
cds = readRDS("stele_cells_labeled.rds")
# create feature column to run linear model
colData(cds)$is_ambiguous = colData(cds)$celltype == "Ambiguous Stele Cells"
lm = fit_models(cds,  model_formula_str = "~is_ambiguous", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "stele_ASCvall_lm.rds") # save model
ct <- coefficient_table(cds_lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_ASCvall_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## Cortex/Enodermis/Lateral Root Endodermis
# load cortex / endodermis data
cds = readRDS("ce_cells_labeled.rds")
# Endodermis vs. Cotex
lm = fit_models(cds[, colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Endodermis" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "stele_ENDODERMISvCORTEX_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_ENDODERMISvCORTEX_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# Lateral Root Endodermis vs. Cortex
lm = fit_models(cds[, colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Lateral Root Endodermis" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "stele_LREvCORTEX_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "stele_LREvCORTEX_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# Endodermis Pseudotime
cds_branch1 <- cds[, !( colData(cds)$Round2_clusters %in% c(1, 2, 5, 7, 8, 10, 13, 14) ) ]
cds_branch2 <- cds[, !( colData(cds)$Round2_clusters %in% c(2, 3, 10, 15) ) ]

cds_branch1 <- order_cells(cds) # this is interactive!!!! Select Endodermis cells as root
colData(cds_branch1)$pseudotime <- cds_branch1@principal_graph_aux$UMAP$pseudotime # saveRDS(cei_r2_branch1, 'cei_branch1.rds')
cds_branch2 <- order_cells(cds_branch2) # this is interactive!!!! Select Endodermis cells as root
colData(cds_branch2)$pseudotime <- cds_branch2@principal_graph_aux$UMAP$pseudotime # saveRDS(cei_r2_branch2, 'cei_branch2.rds')

lm <- fit_models(cds_branch1,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "lm_endodermis_mainbranch_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "endodermis_mainbranch_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

lm <- fit_models(cds_branch2,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "lm_endodermis_lrebranch_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "endodermis_lrebranch_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# plotting
endo_tab = read.table("endodermis_pseudotime.txt", sep="\t", header=T) # slim version of supplemental table 3
plot_cells(cds, genes=endo_tab, cell_size=1.5, label_branch_points=F, label_leaves=F, label_cell_groups=F)+xlim(-5, 6)+ylim(-2.5, 3.5) + coord_equal()  + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))


### MWW
# load pericycle data
cds = readRDS("pericycle_cells_labeled.rds")

# size factor normalized count matrix
norm_mat <- t(t(counts(cds))/size_factors(cds))

## XPP vs. LRP MWW
wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Xylem Pole Pericycle"], x[colData(cds)$celltype == "Lateral Root Primordia"])[['p.value']]))
wilcox_p_values_adj <- p.adjust(wilcox_p_values, method='BH')

# calculate LRP stats
lrp_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Lateral Root Primordia"])
lrp_var <- apply(norm_mat[, colData(cds)$celltype == "Lateral Root Primordia"], 1, var)
lrp_n <- sum(colData(cds)$celltype == "Lateral Root Primordia")

# calculate XPP stats
xpp_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Xylem Pole Pericycle"])
xpp_var <- apply(norm_mat[, colData(cds)$celltype == "Xylem Pole Pericycle"], 1, var)
xpp_n <- sum(colData(cds)$celltype == "Xylem Pole Pericycle")

# create output table
xpp_lrp_log2FoldChange <- log2(xpp_mean/lrp_mean)
out <- cbind(rownames(cds), xpp_mean, sqrt(xpp_var/xpp_n), lrp_mean, sqrt(lrp_var/lrp_n), xpp_lrp_log2FoldChange, wilcox_p_values_adj)
colnames(out) <- c("id", "xpp_mean", "xpp_stand_err", "lrp_mean", "lrp_stand_err", "xpp_lrp_log2FoldChange", "wilcox_p_values_adj")
write.table(out, "XPPvLRP_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

## XPP vs. Mature Pericycle MWW
wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Xylem Pole Pericycle"], x[colData(cds)$celltype == "Mature Pericycle"])[['p.value']]))
wilcox_p_values_adj <- p.adjust(wilcox_p_values, method='BH')

# calculate Mature Pericycle stats
mat_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Mature Pericycle"])
mat_var <- apply(norm_mat[, colData(cds)$celltype == "Mature Pericycle"], 1, var)
mat_n <- sum(colData(cds)$celltype == "Mature Pericycle")

# create output table
xpp_mat_log2FoldChange <- log2(xpp_mean/mat_mean)
out <- cbind(rownames(cds), xpp_mean, sqrt(xpp_var/xpp_n), mat_mean, sqrt(mat_var/mat_n), xpp_mat_log2FoldChange, wilcox_p_values_adj)
colnames(out) <- c("id", "xpp_mean", "xpp_stand_err", "mat_mean", "mat_stand_err", "xpp_mat_log2FoldChange", "wilcox_p_values_adj")
write.table(out, "XPPvMATUREPERICYCLE_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

# load cortex / endodermis data
cds = readRDS("ce_cells_labeled.rds")

# size factor normalized count matrix
norm_mat <- t(t(counts(cds))/size_factors(cds))

## cortex vs. endodermis
endodermis_wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Cortex"], x[colData(cds)$celltype == "Endodermis"])[['p.value']]))
endodermis_wilcox_p_values_adj <- p.adjust(endodermis_wilcox_p_values, method='BH')

## cortex vs. lre
lre_wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Cortex"], x[colData(cds)$celltype == "Lateral Root Endodermis"])[['p.value']]))
lre_wilcox_p_values_adj <- p.adjust(lre_wilcox_p_values, method='BH')

# calculate cortex stats
cortex_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Cortex"])
cortex_var <- apply(norm_mat[, colData(cds)$celltype == "Cortex"], 1, var)
cortex_n <- sum(colData(cds)$celltype == "Cortex")

# calculate endodermis stats
endodermis_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Endodermis"])
endodermis_var <- apply(norm_mat[, colData(cds)$celltype == "Endodermis"], 1, var)
endodermis_n <- sum(colData(cds)$celltype == "Endodermis")
write.table(cbind(rownames(cds), endodermis_mean, cortex_mean, endodermis_wilcox_p_values_adj), "ENDODERMISvCORTEX_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

# calculate lre stats
lre_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Lateral Root Endodermis"])
lre_var <- apply(norm_mat[, colData(cds)$celltype == "Lateral Root Endodermis"], 1, var)
lre_n <- sum(colData(cds)$celltype == "Lateral Root Endodermis")
write.table(cbind(rownames(cds), lre_mean, cortex_mean, lre_wilcox_p_values_adj), "LREvCORTEX_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)
