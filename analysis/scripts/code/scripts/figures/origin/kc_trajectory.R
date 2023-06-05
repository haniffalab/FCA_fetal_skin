#!/usr/bin/env Rscript

# convert anndata to cds
# conda activate sceasy

loadPackage(sceasy)

cds <- sceasy:::anndata2cds(
  "pooled_fsk_org.keratinocytes.count_with_PCA_for_monocle.20220223.h5ad",
  outFile = "pooled_fsk_org.keratinocytes.count_with_PCA_for_monocle.20220223.cds.rds",
  main_layer = "X"
)

# now run monocle
# conda activate monocle3

loadPackage(monocle3)
loadPackage(tidyverse)
loadPackage(viridis)
source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")

# 20220308 analysis
kc_cds <- readRDS(
  "pooled_fsk_org.keratinocytes.count_with_PCA_for_monocle.20220223.cds.rds"
)

colnames(colData(kc_cds))
table(colData(kc_cds)$timepoint)

colData(kc_cds)$pcw <- as.integer(sub("day-", "", colData(kc_cds)$timepoint))
colData(kc_cds)$pcw[colData(kc_cds)$dataset == "org"] <- round(
  colData(kc_cds)$pcw[colData(kc_cds)$dataset == "org"] / 7
)
table(colData(kc_cds)$pcw)

colData(kc_cds)$day <- as.integer(sub("day-", "", colData(kc_cds)$timepoint))
colData(kc_cds)$day[colData(kc_cds)$dataset == "fsk"] <- (
  colData(kc_cds)$pcw[colData(kc_cds)$dataset == "fsk"] * 7
)
table(colData(kc_cds)$day)

fsk_kc_cds <- kc_cds[, colData(kc_cds)$dataset == "fsk"]
org_kc_cds <- kc_cds[, colData(kc_cds)$dataset == "org"]

kc_cds_list0 <- list(fsk = fsk_kc_cds, org = org_kc_cds, pooled = kc_cds)

run_monocle_cds <- function(cds,
                            res = 0.01,
                            umap_n_neighbors = 500L,
                            umap_min_dist = 0.3,
                            use_partition = T,
                            learn_graph_control = list(
                              minimal_branch_len = 30,
                              orthogonal_proj_tip = T
                            )) {
  cds1 <- monocle3::reduce_dimension(
    cds, umap.n_neighbors = umap_n_neighbors, umap.min_dist = umap_min_dist
  )
  cds1 <- monocle3::cluster_cells(cds1, resolution = res)
  cds1 <- monocle3::learn_graph(
    cds1,
    use_partition = use_partition,
    close_loop = T,
    learn_graph_control = learn_graph_control
  )
  cds1
}

kc_cds_list1 <- lapply(kc_cds_list0, run_monocle_cds)
kc_cds_list1[["fsk"]] <- run_monocle_cds(
  kc_cds_list0[["fsk"]], umap_n_neighbors = 100L,
  learn_graph_control = list(
    minimal_branch_len = 10, orthogonal_proj_tip = TRUE
  )
)

for (name in names(kc_cds_list1)) {
  colData(kc_cds_list1[[name]])$timepoint <- as.integer(
    sub("day-", "", colData(kc_cds_list1[[name]])$timepoint)
  )
}

# plot to help choose root and path
pdf(file = "kc_trajectory_helper.pdf", width = 16, height = 9)
for (name in names(kc_cds_list1)) {
  print(
    plot_trajectory(
      kc_cds_list1[[name]],
      color_cells_by = "joint_annot",
      legend_loc = "right",
      label_graph_nodes = T,
      graph_label_size = 2.5,
      cell_size = 50 / sqrt(ncol(kc_cds_list1[[name]]))
    )
    + ggtitle(name)
  )
  print(
    plot_trajectory(
      kc_cds_list1[[name]],
      color_cells_by = "timepoint",
      legend_loc = "right",
      graph_label_size = 2.5,
      continuous_color = T,
      cell_size = 50 / sqrt(ncol(kc_cds_list1[[name]]))
    )
    + ggtitle(name)
  )
}
dev.off()

print(
  plot_trajectory(
    kc_cds_list1[["fsk"]],
    color_cells_by = "joint_annot",
    legend_loc = "right",
    label_graph_nodes = T,
    graph_label_size = 2.5,
    cell_size = 50 / sqrt(ncol(kc_cds_list1[["fsk"]]))
  )
  + ggtitle("fsk")
)
print(
  plot_trajectory(
    kc_cds_list1[["org"]],
    color_cells_by = "joint_annot",
    legend_loc = "right",
    label_graph_nodes = T,
    graph_label_size = 2.5,
    cell_size = 50 / sqrt(ncol(kc_cds_list1[["org"]])))
  + ggtitle("org")
)
print(
  plot_trajectory(
    kc_cds_list1[["pooled"]],
    color_cells_by = "joint_annot",
    legend_loc = "right",
    label_graph_nodes = T,
    graph_label_size = 2.5,
    cell_size = 50 / sqrt(ncol(kc_cds_list1[["pooled"]]))
  )
  + ggtitle("pooled")
)

root_celltypes <- list(
  fsk = c("Immature basal", "Periderm"),
  org = c("Periderm", "Immature basal", "Basal POSTN+"),
  pooled = c("Periderm", "Immature basal", "Basal POSTN+")
)

for (name in names(root_celltypes)) {
  for (root_celltype in root_celltypes[[name]]) {
    for (node_type in c("any", "leaf", "branch")) {
      node <- get_principal_node(
        kc_cds_list1[[name]], "joint_annot", root_celltype, type = node_type
      )
      if (is.null(node)) node <- "none"
      cat(sprintf("%s, %s, %s: %s\n", name, root_celltype, node_type, node))
    }
  }
}
#fsk, Immature basal, any: Y_6
#fsk, Immature basal, leaf: Y_3*
#fsk, Immature basal, branch: none
#fsk, Periderm, any: Y_195
#fsk, Periderm, leaf: Y_195*
#fsk, Periderm, branch: none
#org, Periderm, any: Y_62
#org, Periderm, leaf: Y_62*
#org, Periderm, branch: Y_3
#org, Immature basal, any: Y_1865
#org, Immature basal, leaf: Y_1865 Y_1920*
#org, Immature basal, branch: Y_1862
#org, Basal POSTN+, any: Y_1127
#org, Basal POSTN+, leaf: Y_1127*
#org, Basal POSTN+, branch: Y_1819
#pooled, Periderm, any: Y_2087
#pooled, Periderm, leaf: Y_2087 Y_1899*
#pooled, Periderm, branch: Y_1903
#pooled, Immature basal, any: Y_227
#pooled, Immature basal, leaf: Y_227 Y_134*
#pooled, Immature basal, branch: Y_70
#pooled, Basal POSTN+, any: Y_1766
#pooled, Basal POSTN+, leaf: Y_1766 Y_1868*
#pooled, Basal POSTN+, branch: Y_1753

root_nodes <- list(
  fsk = c("Y_3", "Y_195"),
  org = c("Y_62", "Y_1127", "Y_1920"),
  pooled = c("Y_134", "Y_1868", "Y_1899")
)

for (name in names(root_celltypes)) {
  kc_cds_list1[[name]] <- order_cells(
    kc_cds_list1[[name]],
    root_pr_nodes = root_nodes[[name]]
  )
}

saveRDS(
  kc_cds_list1[["pooled"]], "pooled_fsk_org.keratinocytes.monocle.20220308.rds"
)
saveRDS(kc_cds_list1[["fsk"]], "pooled_fsk.keratinocytes.monocle.20220308.rds")
saveRDS(kc_cds_list1[["org"]], "pooled_org.keratinocytes.monocle.20220308.rds")

kc_cds_list1 <- list(
  pooled = readRDS("pooled_fsk_org.keratinocytes.monocle.20220308.rds"),
  fsk = readRDS("pooled_fsk.keratinocytes.monocle.20220308.rds"),
  org = readRDS("pooled_org.keratinocytes.monocle.20220308.rds")
)

colorby_variables <- list(
  fsk = list(
    discrete = c("joint_annot"), continuous = c("timepoint", "pseudotime")
  ),
  org = list(
    discrete = c("joint_annot"), continuous = c("timepoint", "pseudotime")
  ),
  pooled = list(
    discrete = c("joint_annot"), continuous = c("timepoint", "pseudotime")
  )
)

pdf(
  file = "pooled_fsk_org.keratinocytes.monocle.20220308.pdf",
  width = 10, height = 6
)
for (name in names(colorby_variables)) {
  for (v in colorby_variables[[name]][["discrete"]]) {
    print(plot_trajectory(
      kc_cds_list1[[name]], color_cells_by = v, legend_loc = "right"
    ) + ggtitle(name))
  }
  for (v in colorby_variables[[name]][["continuous"]]) {
    print(plot_trajectory(
      kc_cds_list1[[name]],
      color_cells_by = v,
      legend_loc = "right",
      continuous_color = T
    ) + ggtitle(name))
  }
}
dev.off()

# choose one path
fsk_p1 <- choose_graph_segments(
  kc_cds_list1[["fsk"]],
  starting_pr_node = "Y_3",
  ending_pr_nodes = "Y_99",
  return_list = T
)

fsk_p1_kc_cds <- (
  (kc_cds_list1[["fsk"]])[, colnames(kc_cds_list1[["fsk"]]) %in% fsk_p1$cells]
)

# take only cell type that makes at leat 5% of the cells of the trajectory
fsk_p1_celltypes <- names(
  which(
    table(colData(fsk_p1_kc_cds)$joint_annot) /
      sum(table(colData(fsk_p1_kc_cds)$joint_annot)) >= 0.05
  )
)

fsk_p1_kc_cds <- (
  fsk_p1_kc_cds[, colData(fsk_p1_kc_cds)$joint_annot %in% fsk_p1_celltypes]
)
colData(fsk_p1_kc_cds) <- droplevels(colData(fsk_p1_kc_cds))

plot_trajectory(
  fsk_p1_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 2
)
plot_trajectory(
  kc_cds_list1[["fsk"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
trace("calculateLW", edit = T, where = asNamespace("monocle3"))
fsk_p1_deg <- graph_test(
  fsk_p1_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


fsk_p1_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p1_heatmap <- plot_heatmap(
  fsk_p1_kc_cds,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(50),
  n_cells = 2000,
  cell_annot = c("joint_annot", "timepoint"),
  show_rownames = T,
  order_genes_by = "Spectral"
)

pdf(
  file = "pooled_fsk.mesenchymal.keratinocytes.immature_basal_to_matrix_placode.trajectory_DE.20220308.pdf",
  width = 10,
  height = 6
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap$plot$gtable)
dev.off()

##### END #####
