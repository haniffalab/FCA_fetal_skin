#!/usr/bin/env Rscript

# convert anndata to cds

# library(sceasy)
library(monocle3)
library(tidyverse)
library(viridis)
library(grid)
source("~/sctkr/R/MonocleUtils.R")

data_root <- '/lustre/scratch126/cellgen/team205/nh3/skin'
proj_root <- '~/FCA_Fetal_Skin_priv'

# ms_cds <- sceasy:::anndata2cds(
#   'pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.h5ad',
#   outFile='pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.cds.rds',
#   main_layer='X'
# )

# now run monocle
run_monocle_cds <- function(cds, res=0.01, umap_n_neighbors=500L, umap_min_dist=0.3, use_partition=T,
                            learn_graph_control=list(minimal_branch_len=30, orthogonal_proj_tip=T)) {
    cds1 <- reduce_dimension(cds, umap.n_neighbors=umap_n_neighbors, umap.min_dist=umap_min_dist)
    cds1 <- cluster_cells(cds1, resolution=res)
    cds1 <- learn_graph(cds1, use_partition=use_partition, close_loop=T, learn_graph_control=learn_graph_control)
    cds1
}

# 20220308 analysis
ms_cds <- readRDS(file.path(data_root, '20211022_final_figures', 'monocle', 'pooled_fsk_org.mesenchymal.count_with_PCA_UMAP_for_monocle.20220308.downsampled.cds.rds'))

colnames(colData(ms_cds))
colData(ms_cds)$days <- as.integer(sub("day-", "", colData(ms_cds)$day))
table(colData(ms_cds)$joint_annot)

fsk_ms_cds <- ms_cds[, colData(ms_cds)$dataset == "fsk"]
org_ms_cds <- ms_cds[, colData(ms_cds)$dataset == "org"]

ms_cds_list0 <- list(fsk=fsk_ms_cds, org=org_ms_cds, pooled=ms_cds)

# ms_cds_list1 <- lapply(ms_cds_list0, run_monocle_cds)

# # plot to help choose root and path
# pdf(file=file.path(proj_root, "figures", "obsolete", "ms_trajectory_helper.pdf"), width=16, height=9)
# for (name in names(ms_cds_list1)) {
#   print(
#     plot_trajectory(
#       ms_cds_list1[[name]], color_cells_by="joint_annot", legend_loc="right", label_graph_nodes=T, graph_label_size=2.5)
#     + ggtitle(name)
#   )
# }
# dev.off()

# root_celltypes <- list(
#   fsk=c('HOXC5+ early fibroblast', 'FRZB+ early fibroblast', 'Pericytes', 'Myoblasts'),
#   org=c('HOXC5+ early fibroblast', 'Mural cell'),
#   pooled=c('FRZB+ early fibroblast', 'Myoblasts')
# )

# for (name in names(root_celltypes)) {
#   for (root_celltype in root_celltypes[[name]]) {
#     for (node_type in c("any", "leaf", "branch")) {
#       node <- get_principal_node(ms_cds_list1[[name]], 'joint_annot', root_celltype, type=node_type)
#       if (is.null(node)) node <- "none"
#       cat(sprintf("%s, %s, %s: %s\n", name, root_celltype, node_type, node))
#     }
#   }
# }

# for (name in names(root_celltypes)) {
#   ms_cds_list1[[name]] <- order_cells(
#     ms_cds_list1[[name]],
#     root_pr_nodes=sapply(
#       root_celltypes[[name]],
#       function(root_celltype) get_principal_node(ms_cds_list1[[name]], 'joint_annot', root_celltype, type="leaf")
#     )
#   )
# }

# saveRDS(ms_cds_list1[["pooled"]], "pooled_fsk_org.mesenchymal.monocle.20220308.downsampled.rds")
# saveRDS(ms_cds_list1[["fsk"]], "pooled_fsk.mesenchymal.monocle.20220308.downsampled.rds")
# saveRDS(ms_cds_list1[["org"]], "pooled_org.mesenchymal.monocle.20220308.downsampled.rds")

ms_cds_list1 = list(
  pooled=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_fsk_org.mesenchymal.monocle.20220308.downsampled.rds")),
  fsk=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_fsk.mesenchymal.monocle.20220308.downsampled.rds")),
  org=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_org.mesenchymal.monocle.20220308.downsampled.rds"))
)

colorby_variables <- list(
  fsk=list(discrete=c("joint_annot"), continuous=c("pcw", "pseudotime")),
  org=list(discrete=c("joint_annot"), continuous=c("days", "pseudotime")),
  pooled=list(discrete=c('joint_annot'), continuous=c("pcw", "days", 'pseudotime'))
)

pdf(file=file.path(proj_root, "figures", "obsolete", "pooled_fsk_org.mesenchymal.monocle.20220308.pdf"), width=10, height=6)
for (name in names(colorby_variables)) {
  for (v in colorby_variables[[name]][["discrete"]]) {
    print(plot_trajectory(ms_cds_list1[[name]], color_cells_by=v, legend_loc="right") + ggtitle(name))
  }
  for (v in colorby_variables[[name]][["continuous"]]) {
    print(plot_trajectory(ms_cds_list1[[name]], color_cells_by=v, legend_loc="right") + ggtitle(name))
  }
}
dev.off()

plot_trajectory(ms_cds_list1[["fsk"]], color_cells_by="joint_annot", legend_loc="right", label_graph_nodes=T)

# choose one path
fsk_p1 <- choose_graph_segments(
  ms_cds_list1[["fsk"]],
  starting_pr_node = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Dermal papillia", node_type="leaf"),
  return_list=T
)

fsk_p1_ms_cds <- (ms_cds_list1[["fsk"]])[, colnames(ms_cds_list1[["fsk"]]) %in% fsk_p1$cells]

(fsk_p1_celltypes <- names(which(table(colData(fsk_p1_ms_cds)$joint_annot)/sum(table(colData(fsk_p1_ms_cds)$joint_annot)) >= 0.01)))

fsk_p1_ms_cds <- fsk_p1_ms_cds[, colData(fsk_p1_ms_cds)$joint_annot %in% fsk_p1_celltypes]
colData(fsk_p1_ms_cds) <- droplevels(colData(fsk_p1_ms_cds))

plot_trajectory(fsk_p1_ms_cds, color_cells_by="joint_annot", legend_loc="right")

# DE along path
fsk_p1_deg <- graph_test(fsk_p1_ms_cds, neighbor_graph = "principal_graph", cores = 26)

fsk_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% head(20)

fsk_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

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

pdf(file=file.path(proj_root, "figures", "obsolete", "pooled_fsk.mesenchymal.monocle.hoxc5_early_fibroblast_to_dermal_papillia.trajectory_DE.20220308.pdf"), width=10, height=6)
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
org_p1_deg <- graph_test(org_p1_ms_cds, neighbor_graph = "principal_graph", cores = 26)

org_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% head(20)

org_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
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

pdf(file=file.path(proj_root, "figures", "obsolete", "pooled_org.mesenchymal.monocle.hoxc5_early_fibroblast_to_dermal_papillia.trajectory_DE.20220308.pdf"), width=10, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(org_p1_heatmap$plot$gtable)
dev.off()

###
# choose one path
fsk_p2 <- choose_graph_segments(
  ms_cds_list1[["fsk"]],
  starting_pr_node = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = c(
    get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Dermal papillia", node_type="leaf"),
    get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Adipocytes", node_type="leaf"),
    get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "PEAR1+ fibroblast", node_type="leaf")
  ),
  return_list=T
)

fsk_p2_ms_cds <- (ms_cds_list1[["fsk"]])[, colnames(ms_cds_list1[["fsk"]]) %in% fsk_p2$cells]

(fsk_p2_celltypes <- names(which(table(colData(fsk_p2_ms_cds)$joint_annot)/sum(table(colData(fsk_p2_ms_cds)$joint_annot)) >= 0.01)))

fsk_p2_ms_cds <- fsk_p2_ms_cds[, colData(fsk_p2_ms_cds)$joint_annot %in% fsk_p2_celltypes]
colData(fsk_p2_ms_cds) <- droplevels(colData(fsk_p2_ms_cds))

plot_trajectory(fsk_p2_ms_cds, color_cells_by="joint_annot", legend_loc="right")

###################
# 20220322 analysis
### FSK subset ####
ms_cds_list1 = list(
  pooled=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_fsk_org.mesenchymal.monocle.20220308.downsampled.rds")),
  fsk=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_fsk.mesenchymal.monocle.20220308.downsampled.rds")),
  org=readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_org.mesenchymal.monocle.20220308.downsampled.rds"))
)

# fsk_ms_cds2 <- choose_graph_segments(
#   ms_cds_list1[["fsk"]],
#   starting_pr_node = get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
#   ending_pr_nodes = c(
#     get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Dermal papillia", node_type="leaf"),
#     get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "Adipocytes", node_type="leaf"),
#     get_principal_node(ms_cds_list1[["fsk"]], "joint_annot", "PEAR1+ fibroblast", node_type="leaf")
#   ),
#   return_list=F
# )
# embeds <- SimpleList()
# embeds$PCA <- reducedDims(ms_cds_list1[["fsk"]])$PCA[colnames(ms_cds_list1[["fsk"]]) %in% colnames(fsk_ms_cds2), ]
# embeds$UMAP <- reducedDims(ms_cds_list1[["fsk"]])$UMAP[colnames(ms_cds_list1[["fsk"]]) %in% colnames(fsk_ms_cds2), ]
# reducedDims(fsk_ms_cds2) <- embeds
# rm(embeds)


# fsk_ms_cds2 <- run_monocle_cds(
#   fsk_ms_cds2, learn_graph_control = list(minimal_branch_len = 40, orthogonal_proj_tip = T)
# )

# fsk_ms_cds2 <- fsk_ms_cds2[, colData(fsk_ms_cds2)$joint_annot %in% fsk_ms2_celltypes]
# colData(fsk_ms_cds2) <- droplevels(colData(fsk_ms_cds2))

# fsk_ms_cds2 <- order_cells(
#   fsk_ms_cds2,
#   root_pr_nodes = get_principal_node(
#     fsk_ms_cds2, 'joint_annot', "HOXC5+ early fibroblast", node_type = "leaf"
#   )
# )
# saveRDS(fsk_ms_cds2, "pooled_fsk.mesenchymal_subset.monocle.20220329.downsampled.rds")

fsk_ms_cds2 <- readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_fsk.mesenchymal_subset.monocle.20220329.downsampled.rds"))

(fsk_ms2_celltypes <- names(
  which(
    table(colData(fsk_ms_cds2)$joint_annot) /
      sum(table(colData(fsk_ms_cds2)$joint_annot)) >= 0.01
  )
))

fsk_ms2_celltypes

ms_color_codes <- c(
  `HOXC5+ early fibroblast` = "#BB7784",
  `FRZB+ early fibroblast` = "#D6BDC1",
  `Pre-dermal condensate` = "#EAD3C7",
  `Dermal condensate` = "#7C87B9",
  `Dermal papillia` = "#BFC1D4",
  #`WNT2+ fibroblast` = "#4B67AF",
  `WNT2+ fibroblast` = "#29B8C9",
  `PEAR1+ fibroblast` = "#8794CA",
  `Adipocytes` = "#154496",
  `Mural cell pericytes` = "#C7DFC7",
  `Myoblasts` = "#E6AFBA",
  `Early myocytes` = "#8F143E",
  `Myocytes` = "#E07B91"
)

ms_color_palette <- sapply(fsk_ms2_celltypes, function(ct) ms_color_codes[[ct]])

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
# p0_pdf <- file.path(proj_root, "figures", "final", "pooled_fsk.mesenchymal_subset.monocle.20220329.pdf")
p0_pdf <- file.path(proj_root, "figures", "final", "pooled_fsk.mesenchymal_subset.monocle.20220823.pdf")
pdf(file = p0_pdf, width = 8, height = 6)
plot_trajectory(fsk_ms_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root")
plot_trajectory(fsk_ms_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root", color_palette = ms_color_palette)
plot_trajectory(fsk_ms_cds2, color_cells_by = "pcw", legend_loc = "right", label_graph_nodes = "root")
plot_trajectory(fsk_ms_cds2, color_cells_by = "pseudotime", legend_loc = "right", label_graph_nodes = "root")
#plot_trajectory(fsk_ms_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "all")
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p0_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# choose one path HOXC5 -> dermal papillia
fsk_p1 <- choose_graph_segments(
  fsk_ms_cds2,
  starting_pr_node = get_principal_node(fsk_ms_cds2, "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = get_principal_node(fsk_ms_cds2, "joint_annot", "Dermal papillia", node_type="leaf"),
  return_list=T
)
fsk_p1_cds2 <- fsk_ms_cds2[, colnames(fsk_ms_cds2) %in% fsk_p1$cells]
fsk_p1_cds2 <- fsk_p1_cds2[, colData(fsk_p1_cds2)$joint_annot %in% c("HOXC5+ early fibroblast", "Pre-dermal condensate", "Dermal condensate", "Dermal papillia")]
colData(fsk_p1_cds2) <- droplevels(colData(fsk_p1_cds2))

colData(fsk_p1_cds2)$joint_annot <- factor(
  colData(fsk_p1_cds2)$joint_annot,
  levels=c("HOXC5+ early fibroblast", "Pre-dermal condensate", "Dermal condensate", "Dermal papillia")
)

fsk_p1_deg <- graph_test(fsk_p1_cds2, neighbor_graph = "principal_graph", cores = 26)

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p1_heatmap1 <- plot_heatmap(
  fsk_p1_cds2,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 10),
  rank_threshold = 0.9,
  min_gene_fraction = 0.1,
  min_gene_sd = 0.2,
  smooth_heatmap = 3,
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p1_cds2)$joint_annot)
  ])
)

fsk_p1_heatmap2 <- plot_heatmap(
  fsk_p1_cds2,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.9,
  min_gene_fraction = 0.1,
  min_gene_sd = 0.2,
  smooth_heatmap = 3,
  show_rownames = c(
    "CCL2", "NDUFB8", "QPRT", "LINC01116", "F13A1",
    "CLDN11", "GPC6", "TGFBI", "CXCL12", "CREB5",
    "C2orf40", "GRP", "MTRNR2L8", "LUM", "CAV1",
    "SFRP2", "ABI3BP", "TWIST2", "GPC3", "MEF2C",
    "CCND1", "SLC26A7", "SOX4", "ALPL", "S100A6",
    "PCSK2", "RSPO3", "SOX2", "CYTL1", "CKB", "CRYM"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p1_cds2)$joint_annot)
  ])
)

fsk_p1_heatmap3 <- plot_heatmap(
  fsk_p1_cds2,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.9,
  min_gene_fraction = 0.1,
  min_gene_sd = 0.2,
  smooth_heatmap = 3,
  show_rownames = c(
    "CCL2", "NDUFB8", "QPRT", "LINC01116", "F13A1",
    "CLDN11", "GPC6", "TGFBI", "CXCL12", "CREB5",
    "C2orf40", "GRP", "MTRNR2L8", "LUM", "CAV1",
    "SFRP2", "ABI3BP", "TWIST2", "GPC3", "MEF2C",
    "CCND1", "SLC26A7", "SOX4", "ALPL", "S100A6",
    "PCSK2", "RSPO3", "SOX2", "CYTL1", "CKB", "CRYM"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p1_cds2)$joint_annot)
  ])
)

p1_pdf <- file.path(proj_root, "figures", "final", "pooled_fsk.mesenchymal.monocle.hoxc5_early_fibroblast_to_dermal_papillia.trajectory_DE.20220711.pdf")
pdf(file=p1_pdf, width=8, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap1$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap2$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap3$plot)
plot_trajectory(
  fsk_p1_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root",
  color_palette = ms_color_codes[names(ms_color_codes) %in% unique(colData(fsk_p1_cds2)$joint_annot)]
)
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p1_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# choose one path HOXC5 -> PEAR1+
fsk_p2 <- choose_graph_segments(
  fsk_ms_cds2,
  starting_pr_node = get_principal_node(fsk_ms_cds2, "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = get_principal_node(fsk_ms_cds2, "joint_annot", "PEAR1+ fibroblast", node_type="leaf"),
  return_list=T
)
fsk_p2_cds2 <- fsk_ms_cds2[, colnames(fsk_ms_cds2) %in% fsk_p2$cells]
fsk_p2_cds2 <- fsk_p2_cds2[, colData(fsk_p2_cds2)$joint_annot %in% c("HOXC5+ early fibroblast", "WNT2+ fibroblast", "PEAR1+ fibroblast")]
colData(fsk_p2_cds2) <- droplevels(colData(fsk_p2_cds2))

fsk_p2_deg <- graph_test(fsk_p2_cds2, neighbor_graph = "principal_graph", cores = 26)

colData(fsk_p2_cds2)$joint_annot <- factor(
  colData(fsk_p2_cds2)$joint_annot,
  levels=c("HOXC5+ early fibroblast", "WNT2+ fibroblast", "PEAR1+ fibroblast")
)

fsk_p2_heatmap1 <- plot_heatmap(
  fsk_p2_cds2,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(1000),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 10),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  cell_annot = c("joint_annot", "pcw"),
  smooth_heatmap = 3,
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p2_cds2)$joint_annot)
  ])
)

fsk_p2_heatmap2 <- plot_heatmap(
  fsk_p2_cds2,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(1000),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  cell_annot = c("joint_annot", "pcw"),
  smooth_heatmap = 3,
  show_rownames = c(
    "TUBB2B", "CRABP1", "WIF1", "EFNA5", "INHBA",
    "APCDD1", "QPRT", "CLDN11", "FABP7", "C2orf40",
    "ASPN", "DLK1", "POSTN", "F10", "MGP",
    "PLAC9", "GPC3", "DCN", "OGN", "MFAP5",
    "IGF2", "S100A6", "PI16", "MEOX2", "CLDN1",
    "EGFL6", "SERPINE2", "IGFBP7", "CCDC102B", "APOE"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p2_cds2)$joint_annot)
  ])
)

fsk_p2_heatmap3 <- plot_heatmap(
  fsk_p2_cds2,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(1000),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  cell_annot = c("joint_annot", "pcw"),
  smooth_heatmap = 3,
  show_rownames = c(
    "TUBB2B", "CRABP1", "WIF1", "EFNA5", "INHBA",
    "APCDD1", "QPRT", "CLDN11", "FABP7", "C2orf40",
    "ASPN", "DLK1", "POSTN", "F10", "MGP",
    "PLAC9", "GPC3", "DCN", "OGN", "MFAP5",
    "IGF2", "S100A6", "PI16", "MEOX2", "CLDN1",
    "EGFL6", "SERPINE2", "IGFBP7", "CCDC102B", "APOE"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p2_cds2)$joint_annot)
  ])
)
# maybe the following should be in "final"
p2_pdf <- file.path(proj_root, "figures", "final", "pooled_fsk.mesenchymal.monocle.hoxc5_early_fibroblast_to_pear1_fibroblast.trajectory_DE.20220823.pdf")
pdf(file=p2_pdf, width=8, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p2_heatmap1$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p2_heatmap2$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p2_heatmap3$plot)
plot_trajectory(
  fsk_p2_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root",
  color_palette = ms_color_codes[names(ms_color_codes) %in% unique(colData(fsk_p2_cds2)$joint_annot)]
)
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p2_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# choose one path HOXC5 -> adipocyte
fsk_p3 <- choose_graph_segments(
  fsk_ms_cds2,
  starting_pr_node = get_principal_node(fsk_ms_cds2, "joint_annot", "HOXC5+ early fibroblast", node_type="any"),
  ending_pr_nodes = get_principal_node(fsk_ms_cds2, "joint_annot", "Adipocytes", node_type="leaf"),
  return_list=T
)
fsk_p3_cds2 <- fsk_ms_cds2[, colnames(fsk_ms_cds2) %in% fsk_p3$cells]
fsk_p3_cds2 <- fsk_p3_cds2[, colData(fsk_p3_cds2)$joint_annot %in% c("HOXC5+ early fibroblast", "WNT2+ fibroblast", "Adipocytes")]
colData(fsk_p3_cds2) <- droplevels(colData(fsk_p3_cds2))

colData(fsk_p3_cds2)$joint_annot <- factor(
  colData(fsk_p3_cds2)$joint_annot,
  levels=c("HOXC5+ early fibroblast", "WNT2+ fibroblast", "Adipocytes")
)

fsk_p3_deg <- graph_test(fsk_p3_cds2, neighbor_graph = "principal_graph", cores = 26)

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p3_heatmap1 <- plot_heatmap(
  fsk_p3_cds2,
  fsk_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 10),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  smooth_heatmap = 3,
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p3_cds2)$joint_annot)
  ])
)

fsk_p3_heatmap2 <- plot_heatmap(
  fsk_p3_cds2,
  fsk_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  smooth_heatmap = 3,
  show_rownames = c(
    "WIF1", "CRABP1", "TUBB2B", "INHBA", "CCBE1",
    "QPRT", "APCDD1", "LINC01116", "CLDN11", "TGFBI",
    "CCL2", "ASPN", "DLK1", "MGP", "POSTN",
    "MFAP5", "DCN", "PLAC9", "OGN", "GPC3",
    "IGF2", "CFD", "NUPR1", "SRPX", "VCAN",
    "FABP4", "APOE", "MGST1", "LPL", "CD36"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p3_cds2)$joint_annot)
  ])
)
grid.draw(fsk_p3_heatmap2$plot)

fsk_p3_heatmap3 <- plot_heatmap(
  fsk_p3_cds2,
  fsk_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  cell_annot = c("joint_annot", "pcw"),
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.2,
  smooth_heatmap = 3,
  show_rownames = c(
    "WIF1", "CRABP1", "TUBB2B", "INHBA", "CCBE1",
    "QPRT", "APCDD1", "LINC01116", "CLDN11", "TGFBI",
    "CCL2", "ASPN", "DLK1", "MGP", "POSTN",
    "MFAP5", "DCN", "PLAC9", "OGN", "GPC3",
    "IGF2", "CFD", "NUPR1", "SRPX", "VCAN",
    "FABP4", "APOE", "MGST1", "LPL", "CD36"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ms_color_codes[
    names(ms_color_codes) %in% unique(colData(fsk_p3_cds2)$joint_annot)
  ])
)
# maybe the following should be in "final"
p3_pdf <- file.path(proj_root, "figures", "final", "pooled_fsk.mesenchymal.monocle.hoxc5_early_fibroblast_to_adipocytes.trajectory_DE.20220711.pdf")
pdf(file = p3_pdf, width = 8, height = 6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p3_heatmap1$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p3_heatmap2$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p3_heatmap3$plot)
plot_trajectory(
  fsk_p3_cds2, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root",
  color_palette = ms_color_codes[names(ms_color_codes) %in% unique(colData(fsk_p3_cds2)$joint_annot)]
)
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p3_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))
