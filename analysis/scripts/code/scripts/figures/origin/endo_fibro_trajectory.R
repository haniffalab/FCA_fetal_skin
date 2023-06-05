#!/usr/bin/env Rscript

# # convert anndata to cds
# conda activate sceasy
#
# loadPackage(sceasy)
# 
# cds <- sceasy:::anndata2cds(
#     '../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.raw_autoqc.downsampled.h5ad',
#     outFile = '../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.raw_autoqc.downsampled.counts.cds.rds',
#     main_layer = 'X'
# )
# saveRDS(cds, '../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.raw_autoqc.downsampled.counts.cds.rds')

# # now run monocle
# conda activate monocle3

loadPackage(monocle3)
loadPackage(reticulate)
loadPackage(plyr)

reticulate::use_condaenv('monocle3')
anndata <- reticulate::import('anndata')

get_earliest_principal_node <- function(cds, group_by, root_group){
  cell_ids <- which(colData(cds)[, group_by] == root_group)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

run_monocle <- function(cds, batch=NA, add_pca=NA, add_umap=NA, res=1e-2, use_partition=F, learn_graph_control=NULL) {
    cds@colData <- droplevels(cds@colData)
    cds1 <- preprocess_cds(cds, num_dim=50)
    if (!is.na(batch)) {
        cds1 <- align_cds(cds1, num_dim=50, alignment_group=batch)
    }
    if (!is.na(add_pca) || !is.na(add_umap)) {
        embeds <- SimpleList()
        if (!is.na(add_pca)) embeds$PCA <- add_pca
        if (!is.na(add_umap)) embeds$UMAP <- add_umap
        reducedDims(cds1) <- embeds
    }
    if (is.na(add_umap)) {
        cds1 <- reduce_dimension(cds1)
    }
    cds1 <- cluster_cells(cds1, resolution=res)
    cds1 <- learn_graph(cds1, use_partition=use_partition, learn_graph_control=learn_graph_control)
    cds1
}

# read data

cds <- readRDS('../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.raw_autoqc.downsampled.counts.cds.rds')
str(cds)
dim(cds)

# Stroma + endo

ad <- anndata$read_h5ad('../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.processed.downsampled.h5ad')
X_pca_hm <- ad$obsm['X_pca_hm']
X_umap_bk <- ad$obsm['X_umap_bk']

cds1 <- run_monocle(cds, add_pca=X_pca_hm, res=0.01, learn_graph_control=list(minimal_branch_len=15))
print(plot_cells(cds1, color_cells_by='fsk_annot', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
cds1 <- order_cells(cds1, root_pr_nodes=c(
    get_earliest_principal_node(cds1, 'fsk_annot', 'fetal_Early fibroblast HOXC5+'),
    get_earliest_principal_node(cds1, 'fsk_annot', 'fetal_Early fibroblast FRZB+')
))
pdf(file=paste0('organoid_fibro_endo_subset.monocle3_harmony_trajectory.pdf'), height=5, width=5)
print(plot_cells(cds1, color_cells_by='fsk_annot', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
print(plot_cells(cds1, color_cells_by='fsk_annot', label_cell_groups=T,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
print(plot_cells(cds1, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
print(plot_cells(cds1, color_cells_by='day', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
print(plot_cells(cds1, color_cells_by='batch', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
dev.off()





k_fibro1 <- colData(fetal_cds)$independent_annotation_refined %in% c(
    'Dermal condensate', 'Dermal papillia', 'Early fibroblast FRZB+', 'Early fibroblast HOXC5+', 'Fibroblast WNT2+', 'Pre-dermal condensate',
    'Pericytes', 'Smooth muscle LMCD1+', 'Smooth muscle PLN+', 'Fibroblasts unknown', 'Adipocytes', 'Myofibroblasts'
)
k_fibro2 <- colData(fetal_cds)$independent_annotation_refined %in% c(
    'Dermal condensate', 'Dermal papillia', 'Early fibroblast FRZB+', 'Early fibroblast HOXC5+', 'Fibroblast WNT2+', 'Pre-dermal condensate',
    'Pericytes', 'Smooth muscle LMCD1+', 'Smooth muscle PLN+', 'Fibroblasts unknown'
)
k_fibro3 <- colData(fetal_cds)$independent_annotation_refined %in% c(
    'Dermal condensate', 'Dermal papillia', 'Early fibroblast FRZB+', 'Early fibroblast HOXC5+', 'Fibroblast WNT2+', 'Pre-dermal condensate',
    'Pericytes', 'Smooth muscle LMCD1+', 'Smooth muscle PLN+'
)
k_fibro4 <- colData(fetal_cds)$independent_annotation_refined %in% c(
    'Dermal condensate', 'Dermal papillia', 'Early fibroblast FRZB+', 'Early fibroblast HOXC5+', 'Fibroblast WNT2+', 'Pre-dermal condensate'
)
k_fibro5 <- colData(fetal_cds)$independent_annotation_refined %in% c(
    'Dermal condensate', 'Dermal papillia', 'Early fibroblast HOXC5+', 'Fibroblast WNT2+', 'Pre-dermal condensate'
)
k_fibro_list <- list(k_fibro1, k_fibro2, k_fibro3, k_fibro4, k_fibro5)

#sapply(k_fibro_list, sum)

for (i in seq_along(k_fibro_list)) {
    gc()
    k_fibro <- k_fibro_list[[i]]
    n_fibro <- sum(k_fibro)
    set.seed(0)
    rand_fibro_idx <- which(k_fibro)[sort(sample(n_fibro, size=round(n_fibro/10)))]
    fetal_fibro_cds0 <- fetal_cds[, rand_fibro_idx]

    fetal_fibro_cds1 <- fetal_fibro_cds0[, colnames(fetal_fibro_cds0) %in% fibro_cell_ids]
    fetal_fibro_cds1@colData <- droplevels(fetal_fibro_cds1@colData)
    k_cds1 <- fibro_cell_ids %in% colnames(fetal_fibro_cds1)
    fibro_X_pca_hm <- fibro_ad$obsm['X_pca_hm'][k_cds1, ]
    fibro_X_umap_bk <- fibro_ad$obsm['X_umap_bk'][k_cds1, ]

    colData(fetal_fibro_cds1)$fig3_annotation <- revalue(
        colData(fetal_fibro_cds1)$independent_annotation_refined,
        c('Fibroblasts unknown'='Fibroblast PEAR1+', 'Smooth muscle LMCD1+'='Pericyte LMCD1+', 'Smooth muscle PLN+'='Pericyte PLN+')
    )

    fetal_fibro_cds1 <- run_monocle(fetal_fibro_cds1, add_pca=fibro_X_pca_hm, res=0.01, learn_graph_control=list(minimal_branch_len=15))
    fetal_fibro_cds1 <- order_cells(fetal_fibro_cds1, root_pr_nodes=c(
        get_earliest_principal_node(fetal_fibro_cds1, 'fig3_annotation', 'Early fibroblast HOXC5+'),
        get_earliest_principal_node(fetal_fibro_cds1, 'fig3_annotation', 'Early fibroblast FRZB+')
    ))

    fetal_fibro_cds1b <- learn_graph(fetal_fibro_cds1, use_partition=T, learn_graph_control=list(minimal_branch_len=15))
    fetal_fibro_cds1b <- order_cells(fetal_fibro_cds1b, root_pr_nodes=c(
        get_earliest_principal_node(fetal_fibro_cds1b, 'fig3_annotation', 'Early fibroblast HOXC5+'),
        get_earliest_principal_node(fetal_fibro_cds1b, 'fig3_annotation', 'Early fibroblast FRZB+')
    ))

    pdf(file=paste0('fetal_fibro_subset', i, '.monocle3_harmony_trajectory.pdf'), height=5, width=5)
    print(plot_cells(fetal_fibro_cds1, color_cells_by='fig3_annotation', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='fig3_annotation', label_cell_groups=T,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, genes=c('HOXC5', 'FRZB', 'WNT', 'LPL', 'ACTA1'), label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=0.5, show_trajectory_graph=F))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='pcw', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='batch', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='chemistry', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1, color_cells_by='sorting', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    dev.off()

    pdf(file=paste0('fetal_fibro_subset', i, '.monocle3_harmony_trajectory2.pdf'), height=5, width=5)
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='fig3_annotation', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='fig3_annotation', label_cell_groups=T,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, genes=c('HOXC5', 'FRZB', 'WNT', 'LPL', 'ACTA1'), label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=0.5, show_trajectory_graph=F))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='pcw', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='batch', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='chemistry', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    print(plot_cells(fetal_fibro_cds1b, color_cells_by='sorting', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    dev.off()
}

    
rm(list=ls())
