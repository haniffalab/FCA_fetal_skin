#!/usr/bin/env Rscript

# convert anndata to cds

loadPackage(sceasy)

ms_cds <- sceasy:::anndata2cds(
  'pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.h5ad',
  outFile='pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.cds.rds',
  main_layer='X'
)

# now run monocle

loadPackage(monocle3)
loadPackage(tidyverse)
loadPackage(viridis)
source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")

# 20220308 analysis
ms_cds <- readRDS('pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.cds.rds')

colnames(colData(ms_cds))
colData(ms_cds)$days <- as.integer(sub("day-", "", colData(ms_cds)$day))
table(colData(ms_cds)$joint_annot)

fsk_ms_cds <- ms_cds[, colData(ms_cds)$dataset == "fsk"]
org_ms_cds <- ms_cds[, colData(ms_cds)$dataset == "org"]

ms_cds_list0 <- list(fsk=fsk_ms_cds, org=org_ms_cds, pooled=ms_cds)

run_monocle_cds <- function(cds, res=0.01, umap_n_neighbors=500L, umap_min_dist=0.3, use_partition=T,
                            learn_graph_control=list(minimal_branch_len=30, orthogonal_proj_tip=T)) {
    cds1 <- reduce_dimension(cds, umap.n_neighbors=umap_n_neighbors, umap.min_dist=umap_min_dist)
    cds1 <- cluster_cells(cds1, resolution=res)
    cds1 <- learn_graph(cds1, use_partition=use_partition, close_loop=T, learn_graph_control=learn_graph_control)
    cds1
}

ms_cds_list1 <- lapply(ms_cds_list0, run_monocle_cds)

# plot to help choose root and path
pdf(file="ms_trajectory_helper.pdf", width=16, height=9)
for (name in names(ms_cds_list1)) {
  print(
    plot_trajectory(
      ms_cds_list1[[name]], color_cells_by="joint_annot", legend_loc="right", label_graph_nodes=T, graph_label_size=2.5)
    + ggtitle(name)
  )
}
dev.off()

root_celltypes <- list(
  fsk=c('HOXC5+ early fibroblast', 'FRZB+ early fibroblast', 'Pericytes', 'Myoblasts'),
  org=c('HOXC5+ early fibroblast', 'Mural cell'),
  pooled=c('FRZB+ early fibroblast', 'Myoblasts')
)

for (name in names(root_celltypes)) {
  for (root_celltype in root_celltypes[[name]]) {
    for (node_type in c("any", "leaf", "branch")) {
      node <- get_principal_node(ms_cds_list1[[name]], 'joint_annot', root_celltype, type=node_type)
      if (is.null(node)) node <- "none"
      cat(sprintf("%s, %s, %s: %s\n", name, root_celltype, node_type, node))
    }
  }
}

for (name in names(root_celltypes)) {
  ms_cds_list1[[name]] <- order_cells(
    ms_cds_list1[[name]],
    root_pr_nodes=sapply(
      root_celltypes[[name]],
      function(root_celltype) get_principal_node(ms_cds_list1[[name]], 'joint_annot', root_celltype, type="leaf")
    )
  )
}

saveRDS(ms_cds_list1[["pooled"]], "pooled_fsk_org.mesenchymal.monocle.20220308.downsampled.rds")
saveRDS(ms_cds_list1[["fsk"]], "pooled_fsk.mesenchymal.monocle.20220308.downsampled.rds")
saveRDS(ms_cds_list1[["org"]], "pooled_org.mesenchymal.monocle.20220308.downsampled.rds")

ms_cds_list1 = list(
  pooled=readRDS("pooled_fsk_org.mesenchymal.monocle.20220308.downsampled.rds"),
  fsk=readRDS("pooled_fsk.mesenchymal.monocle.20220308.downsampled.rds"),
  org=readRDS("pooled_org.mesenchymal.monocle.20220308.downsampled.rds")
)

colorby_variables <- list(
  fsk=list(discrete=c("joint_annot"), continuous=c("pcw", "pseudotime")),
  org=list(discrete=c("joint_annot"), continuous=c("days", "pseudotime")),
  pooled=list(discrete=c('joint_annot'), continuous=c("pcw", "days", 'pseudotime'))
)

pdf(file="pooled_fsk_org.mesenchymal.monocle.20220308.pdf", width=10, height=6)
for (name in names(colorby_variables)) {
  for (v in colorby_variables[[name]][["discrete"]]) {
    print(plot_trajectory(ms_cds_list1[[name]], color_cells_by=v, legend_loc="right") + ggtitle(name))
  }
  for (v in colorby_variables[[name]][["continuous"]]) {
    print(plot_trajectory(ms_cds_list1[[name]], color_cells_by=v, legend_loc="right", continuous_color=T) + ggtitle(name))
  }
}
dev.off()

plot_trajectory(ms_cds, color_cells_by="joint_annot", legend_loc="right", label_graph_nodes=T)

# choose one path
fsk_p1 <- choose_graph_segments(
  ms_cds_list1[["fsk"]],
  starting_pr_node = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "HOXC5+ early fibroblast", type="any"),
  ending_pr_nodes = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Dermal papillia", type="leaf"),
  return_list=T
)

fsk_p1_ms_cds <- (ms_cds_list1[["fsk"]])[, colnames(ms_cds_list1[["fsk"]]) %in% fsk_p1$cells]

fsk_p1_celltypes <- names(which(table(colData(fsk_p1_ms_cds)$joint_annot)/sum(table(colData(fsk_p1_ms_cds)$joint_annot)) >= 0.01))

fsk_p1_ms_cds <- fsk_p1_ms_cds[, colData(fsk_p1_ms_cds)$joint_annot %in% fsk_p1_celltypes]
colData(fsk_p1_ms_cds) <- droplevels(colData(fsk_p1_ms_cds))

plot_trajectory(fsk_p1_ms_cds, color_cells_by="joint_annot", legend_loc="right")

# DE along path
fsk_p1_deg <- graph_test(fsk_p1_ms_cds, neighbor_graph = "principal_graph", cores = 8)

fsk_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% head(20)

fsk_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p1_heatmap <- plot_heatmap(
  fsk_p1_ms_cds,
  fsk_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(40),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "TSP"
)
fsk_p1_heatmap <- plot_heatmap(
  fsk_p1_ms_cds,
  fsk_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(40),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "corr"
)

pdf(file="pooled_fsk.mesenchymal.monocle.hoxc5_early_fibroblast_to_dermal_papillia.trajectory_DE.20220308.pdf", width=10, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap$plot$gtable)
dev.off()

##### END #####

org_p1 <- choose_graph_segments(
  ms_cds_list1[["org"]],
  starting_pr_node = get_principal_node(ms_cds_list1[["org"]], "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = get_principal_node(ms_cds_list1[["org"]], "joint_annot", "Dermal papillia", node_type="leaf"),
  return_list=T
)

org_p1_ms_cds <- (ms_cds_list1[["org"]])[, colnames(ms_cds_list1[["org"]]) %in% org_p1$cells]

org_p1_celltypes <- names(which(table(colData(org_p1_ms_cds)$joint_annot)/sum(table(colData(org_p1_ms_cds)$joint_annot)) >= 0.01))
org_p1_celltypes

org_p1_ms_cds <- org_p1_ms_cds[, colData(org_p1_ms_cds)$joint_annot %in% org_p1_celltypes]
colData(org_p1_ms_cds) <- droplevels(colData(org_p1_ms_cds))

plot_trajectory(org_p1_ms_cds, color_cells_by="joint_annot", legend_loc="right")

# DE along path
org_p1_deg <- graph_test(org_p1_ms_cds, neighbor_graph = "principal_graph", cores = 8)

org_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% head(20)

org_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
org_p1_heatmap <- plot_heatmap(
  org_p1_ms_cds,
  org_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(40),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "TSP"
)
org_p1_heatmap <- plot_heatmap(
  org_p1_ms_cds,
  org_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(40),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "corr"
)

pdf(file="pooled_org.mesenchymal.monocle.hoxc5_early_fibroblast_to_dermal_papillia.trajectory_DE.20220308.pdf", width=10, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(org_p1_heatmap$plot$gtable)
dev.off()
