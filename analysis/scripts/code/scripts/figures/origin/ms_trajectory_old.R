#!/usr/bin/env Rscript

# convert anndata to cds
#conda activate sceasy

loadPackage(sceasy)

cds <- sceasy:::anndata2cds('pooled_adult_fetal_skin.counts.h5ad', outFile='pooled_adult_fetal_skin.counts.cds', main_layer='X')
#cds <- sceasy:::anndata2cds('pooled_fsk_org.keratinocytes.count_with_PCA_for_monocle.20220214.h5ad', outFile='pooled_fsk_org.keratinocytes.count_with_PCA_for_monocle.20220214.cds.rds', main_layer='X')

# now run monocle
#conda activate monocle3

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

run_monocle_cds <- function(cds, batch=NA, add_pca=NA, add_umap=NA, res=1e-2, use_partition=F, learn_graph_control=NULL) {
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


run_monocle <- function(
                        input_obj,
                        batch=NULL,
                        add_pca=NULL,
                        add_umap=NULL,
                        res=1e-2,
                        use_partition=T,
                        learn_graph_control=NULL) {
  if (class(obj) == "cell_data_set") {
    cds <- input_obj
  } else if (class(obj) == "anndata._core.anndata.AnnData") {
    if (!requireNamespace("reticulate", quitely=TRUE)) {
      stop("Package `reticulate` required if `input_obj` is an `AnnData`")
    }
    builtins <- reticulate::import_builtins(convert=FALSE)
    anndata <- reticulate::import("anndata", convert=FALSE)
    scipy_sp <- reticulate::import("scipy.sparse", convert=FALSE)
    ad <- input_obj

    cds <- new_cell_data_set(
      expression_data=Matrix::t(reticulate::py_to_r(scipy_sp$csc_matrix(ad$X)))
    )
  } else {
    stop("Unsupported `input_obj` type. Please supply a CDS or AnnData (by reticulate)")
  }
}



# read data

cds <- readRDS('pooled_adult_fetal_skin.counts.cds.rds')
fetal_cds <- cds[, cds@colData@listData$dataset == 'fetal']
fetal_cds@colData <- droplevels(fetal_cds@colData[, c(1:27, 41:43)])
rm(cds)

# Stroma

fibro_ad <- anndata$read_h5ad('fetal_skin.stroma.maternal_removed.norm.20211120.h5ad')
fibro_cell_ids <- as.vector(sub('$', '-fetal', fibro_ad$obs_names$values))

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

for (i in seq_along(k_fibro_list)) {
    gc()
    i <- 5
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

    colData(fetal_fibro_cds1)$fig3_annotation <- plyr::revalue(
        colData(fetal_fibro_cds1)$independent_annotation_refined,
        c('Fibroblasts unknown'='Fibroblast PEAR1+', 'Smooth muscle LMCD1+'='Pericyte LMCD1+', 'Smooth muscle PLN+'='Pericyte PLN+')
    )

    fetal_fibro_cds1 <- run_monocle_cds(fetal_fibro_cds1, add_pca=fibro_X_pca_hm, res=0.01, learn_graph_control=list(minimal_branch_len=15))
    colnames(colData(fetal_fibro_cds1))
    fetal_fibro_cds1 <- order_cells(fetal_fibro_cds1, root_pr_nodes=c(
        get_earliest_principal_node(fetal_fibro_cds1, 'fig3_annotation', 'Early fibroblast HOXC5+')
    ))

    k_cds2 <- !(colnames(fetal_fibro_cds1) %in% outlier_cell_ids)
    fetal_fibro_cds1a <- fetal_fibro_cds1[, k_cds2]
    fetal_fibro_cds1a@colData <- droplevels(fetal_fibro_cds1a@colData)
    fibro_X_pca_hm <- fibro_X_pca_hm[k_cds2, ]
    fetal_fibro_cds1a <- run_monocle_cds(fetal_fibro_cds1a, add_pca=fibro_X_pca_hm, res=0.005, learn_graph_control=list(minimal_branch_len=15))
    fetal_fibro_cds1b <- learn_graph(fetal_fibro_cds1a, use_partition=F, learn_graph_control=list(minimal_branch_len=25))
    fetal_fibro_cds1b <- order_cells(fetal_fibro_cds1b, root_pr_nodes=c(
        get_earliest_principal_node(fetal_fibro_cds1b, 'fig3_annotation', 'Early fibroblast HOXC5+')
    ))
    #fetal_fibro_cds1b <- learn_graph(fetal_fibro_cds1, use_partition=t, learn_graph_control=list(minimal_branch_len=15))
    #fetal_fibro_cds1b <- order_cells(fetal_fibro_cds1b, root_pr_nodes=c(
    #    get_earliest_principal_node(fetal_fibro_cds1b, 'fig3_annotation', 'early fibroblast hoxc5+')
    #))
    break
}
saveRDS(fetal_fibro_cds1, "fetal_fibro_subset5.monocle3.20220314.cds.rds")
saveRDS(fetal_fibro_cds1a, "fetal_fibro_subset5.monocle3.20220322.cds.rds")

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")

fetal_fibro_cds1a <- readRDS("fetal_fibro_subset5.monocle3.20220322.cds.rds")
plot_trajectory(fetal_fibro_cds1a, color_cells_by="fig3_annotation", cell_size=1, label_graph_nodes=F, legend_loc="right")
plot_trajectory(fetal_fibro_cds1b, color_cells_by="fig3_annotation", cell_size=1, label_graph_nodes=T, legend_loc="right")
dev.copy2pdf(file="fetal_fibro_subset5.monocle3.20220322.pdf")

fibro_umap <- reducedDims(fetal_fibro_cds1)$UMAP
outlier_cell_ids <- rownames(fibro_umap)[fibro_umap[, 1] > 15]
class(fibro_X_pca_hm)
class(fibro_X_pca_hm[k_cds2, ])

fetal_p1 <- choose_graph_segments(
  fetal_fibro_cds1b,
  starting_pr_node = "Y_78",
  ending_pr_nodes = c("Y_1095"),
  return_list = T
)

fetal_p1_fibro_cds <- fetal_fibro_cds1b[, colnames(fetal_fibro_cds1b) %in% fetal_p1$cells]

table(colData(fetal_p1_fibro_cds)$fig3_annotation)/sum(table(colData(fetal_p1_fibro_cds)$fig3_annotation))

fetal_p1_celltypes <- names(which(table(colData(fetal_p1_fibro_cds)$fig3_annotation)/sum(table(colData(fetal_p1_fibro_cds)$fig3_annotation)) >= 0.05))
fetal_p1_fibro_cds1 <- fetal_p1_fibro_cds[, colData(fetal_p1_fibro_cds)$fig3_annotation %in% fetal_p1_celltypes]
colData(fetal_p1_fibro_cds1) <- droplevels(colData(fetal_p1_fibro_cds1))

plot_trajectory(fetal_p1_fibro_cds, color_cells_by="fig3_annotation", cell_size=1, label_graph_nodes=T)

fetal_p1_deg <- graph_test(fetal_p1_fibro_cds1, neighbor_graph = "principal_graph", cores = 8)

fetal_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

fetal_p1_heatmap <- plot_heatmap(
  fetal_p1_fibro_cds1,
  fetal_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(60),
  n_cells = 2000,
  cell_annot = c("fig3_annotation", "pcw"),
  show_rownames = T,
  order_genes_by = "rank"
)
fetal_p1_heatmap$plot

pdf(file="fetal_fibro_subset5.monocle3.20220322.pdf", width=7, height=6)
plot_trajectory(fetal_fibro_cds1b, color_cells_by="fig3_annotation", cell_size=1, legend_loc="right")
plot_trajectory(fetal_fibro_cds1b, color_cells_by="pseudotime", continuous_color=T, cell_size=1)
plot_trajectory(fetal_fibro_cds1b, color_cells_by="pcw", continuous_color=T, cell_size=1)
plot_trajectory(fetal_fibro_cds1b, color_cells_by="fig3_annotation", cell_size=1, label_graph_nodes=T, legend_loc=F)
plot_trajectory(fetal_p1_fibro_cds, color_cells_by="fig3_annotation", cell_size=1, label_graph_nodes=T)
dev.off()

pdf(file="fetal_fibro_subset5.fetal_hair_folicle.hoxc5_fibro_to_dermal_papillia.monocle3.trajectory_DE.20220322.pdf", width=10, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p1_heatmap$plot$gtable)
dev.off()

    #pdf(file=paste0('fetal_fibro_subset', i, '.monocle3_harmony_trajectory.pdf'), height=5, width=5)
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='fig3_annotation', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='fig3_annotation', label_cell_groups=T,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, genes=c('HOXC5', 'FRZB', 'WNT', 'LPL', 'ACTA1'), label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=0.5, show_trajectory_graph=F))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='pcw', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='batch', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='chemistry', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1, color_cells_by='sorting', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #dev.off()

    #pdf(file=paste0('fetal_fibro_subset', i, '.monocle3_harmony_trajectory2.pdf'), height=5, width=5)
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='fig3_annotation', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='fig3_annotation', label_cell_groups=T,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, genes=c('HOXC5', 'FRZB', 'WNT', 'LPL', 'ACTA1'), label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=0.5, show_trajectory_graph=F))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='pseudotime', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='pcw', label_cell_groups=F, label_leaves=F, label_branch_points=F, graph_label_size=1.5, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='batch', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='chemistry', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #print(plot_cells(fetal_fibro_cds1b, color_cells_by='sorting', label_cell_groups=F,  label_groups_by_cluster=F, label_leaves=F, label_branch_points=F, cell_size=1))
    #dev.off()
#}

    
rm(list=ls())

getwd()
cds <- readRDS('../20211214_include_orgnoid_day6/organoid_day6to133_fibro_endo.raw_autoqc.counts.cds.rds')
