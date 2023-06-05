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

# 20220315
kc_cds <- kc_cds[, !colData(kc_cds)$joint_annot %in% c(
  "Immature basal", "Immature suprabasal", "Periderm", "Suprabasal IFE"
)]
colData(kc_cds) <- droplevels(colData(kc_cds))

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
  kc_cds_list0[["fsk"]],
  umap_n_neighbors = 100L,
  umap_min_dist = 0.2,
  learn_graph_control = list(
    minimal_branch_len = 5, orthogonal_proj_tip = TRUE
  )
)

#for (name in names(kc_cds_list1)) {
#  colData(kc_cds_list1[[name]])$timepoint <- as.integer(
#    sub("day-", "", colData(kc_cds_list1[[name]])$timepoint)
#  )
#}

colData(kc_cds) %>% colnames()

# plot to help choose root and path
pdf(file = "kc_trajectory_helper.20220315.pdf", width = 16, height = 9)
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
      color_cells_by = "pcw",
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
  fsk = c("Basal POSTN+"),
  org = c("Basal POSTN+"),
  pooled = c("Basal POSTN+")
)

for (name in names(root_celltypes)) {
  for (root_celltype in root_celltypes[[name]]) {
    for (node_type in c("any", "leaf", "branch")) {
      node <- get_principal_node(
        kc_cds_list1[[name]],
        "joint_annot",
        root_celltype,
        node_type = node_type
        #choose_by = "pcw",
        #ascending = TRUE
      )
      if (is.null(node) || is.na(node)) node <- "none"
      cat(sprintf("%s, %s, %s: %s\n", name, root_celltype, node_type, node))
    }
  }
}
#fsk, Basal POSTN+, any: Y_161
#fsk, Basal POSTN+, leaf: Y_85*
#fsk, Basal POSTN+, branch: none
#org, Basal POSTN+, any: Y_490
#org, Basal POSTN+, leaf: Y_490*
#org, Basal POSTN+, branch: Y_999
#pooled, Basal POSTN+, any: Y_705
#pooled, Basal POSTN+, leaf: Y_705*
#pooled, Basal POSTN+, branch: Y_853

root_nodes <- list(
  fsk = c("Y_85"),
  org = c("Y_490"),
  pooled = c("Y_705")
)

for (name in names(root_celltypes)) {
  kc_cds_list1[[name]] <- order_cells(
    kc_cds_list1[[name]],
    root_pr_nodes = root_nodes[[name]]
  )
}

saveRDS(
  kc_cds_list1[["pooled"]], "pooled_fsk_org.keratinocytes.monocle.20220315.rds"
)
saveRDS(kc_cds_list1[["fsk"]], "pooled_fsk.keratinocytes.monocle.20220315.rds")
saveRDS(kc_cds_list1[["org"]], "pooled_org.keratinocytes.monocle.20220315.rds")

kc_cds_list1 <- list(
  pooled = readRDS("pooled_fsk_org.keratinocytes.monocle.20220308.rds"),
  fsk = readRDS("pooled_fsk.keratinocytes.monocle.20220308.rds"),
  org = readRDS("pooled_org.keratinocytes.monocle.20220308.rds")
)
dim(kc_cds_list1[["fsk"]])

colorby_variables <- list(
  fsk = list(
    discrete = c("joint_annot"), continuous = c("pcw", "pseudotime")
  ),
  org = list(
    discrete = c("joint_annot"), continuous = c("day", "pseudotime")
  ),
  pooled = list(
    discrete = c("joint_annot"), continuous = c("pcw", "pseudotime")
  )
)

pdf(
  file = "pooled_fsk_org.keratinocytes.monocle.20220315.pdf",
  width = 10, height = 6
)
for (name in names(colorby_variables)) {
  for (v in colorby_variables[[name]][["discrete"]]) {
    print(plot_trajectory(
      kc_cds_list1[[name]], color_cells_by = v, legend_loc = "right",
      cell_size = 50 / sqrt(ncol(kc_cds_list1[[name]]))
    ) + ggtitle(name))
  }
  for (v in colorby_variables[[name]][["continuous"]]) {
    print(plot_trajectory(
      kc_cds_list1[[name]],
      color_cells_by = v,
      legend_loc = "right",
      continuous_color = T,
      cell_size = 50 / sqrt(ncol(kc_cds_list1[[name]]))
    ) + ggtitle(name))
  }
}
dev.off()

# choose one path
source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
get_principal_node(
  kc_cds_list1[["fsk"]], "joint_annot", "Matrix/placode", node_type = "leaf",
  choose_by = "pcw", ascending = FALSE
)
get_principal_node(
  kc_cds_list1[["fsk"]], "joint_annot", "Matrix/placode", node_type = "leaf",
  choose_one = FALSE
)

get_principal_node(
  kc_cds_list1[["fsk"]], "joint_annot", "Inner root sheath", node_type = "leaf"
)
# POSTN+ basal -> Inner root sheath
fsk_p1 <- choose_graph_segments(
  kc_cds_list1[["fsk"]],
  starting_pr_node = "Y_85",
  ending_pr_nodes = "Y_38",
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
fsk_p1_celltypes

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
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_fsk.keratinocytes.POSTN_basal_to_inner_root_sheath.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p1_heatmap$plot$gtable)
dev.off()

# POSTN+ basal -> Companion layer
fsk_p2 <- choose_graph_segments(
  kc_cds_list1[["fsk"]],
  starting_pr_node = "Y_85",
  ending_pr_nodes = "Y_213",
  return_list = T
)

fsk_p2_kc_cds <- (
  (kc_cds_list1[["fsk"]])[, colnames(kc_cds_list1[["fsk"]]) %in% fsk_p2$cells]
)
dim(fsk_p2_kc_cds)

# take only cell type that makes at leat 5% of the cells of the trajectory
fsk_p2_celltypes <- names(
  which(
    table(colData(fsk_p2_kc_cds)$joint_annot) /
      sum(table(colData(fsk_p2_kc_cds)$joint_annot)) >= 0.05
  )
)
fsk_p2_celltypes

fsk_p2_kc_cds <- (
  fsk_p2_kc_cds[, colData(fsk_p2_kc_cds)$joint_annot %in% fsk_p2_celltypes]
)
colData(fsk_p2_kc_cds) <- droplevels(colData(fsk_p2_kc_cds))

plot_trajectory(
  fsk_p2_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)
plot_trajectory(
  kc_cds_list1[["fsk"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
#trace("calculateLW", edit = T, where = asNamespace("monocle3"))
fsk_p2_deg <- graph_test(
  fsk_p2_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


fsk_p2_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p2_heatmap <- plot_heatmap(
  fsk_p2_kc_cds,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_fsk.keratinocytes.POSTN_basal_to_companion_layer.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p2_heatmap$plot$gtable)
dev.off()

# POSTN+ basal -> Matrix/placode
fsk_p3 <- choose_graph_segments(
  kc_cds_list1[["fsk"]],
  starting_pr_node = "Y_85",
  #ending_pr_nodes = get_principal_node(kc_cds_list1[["fsk"]], "joint_annot", "Matrix/placode", node_type = "leaf", choose_one = FALSE),
  ending_pr_nodes = "Y_187",
  return_list = T
)
get_principal_node(kc_cds_list1[["fsk"]], "joint_annot", "Matrix/placode", node_type = "leaf", choose_one = FALSE)

fsk_p3_kc_cds <- (
  (kc_cds_list1[["fsk"]])[, colnames(kc_cds_list1[["fsk"]]) %in% fsk_p3$cells]
)

# take only cell type that makes at leat 5% of the cells of the trajectory
fsk_p3_celltypes <- names(
  which(
    table(colData(fsk_p3_kc_cds)$joint_annot) /
      sum(table(colData(fsk_p3_kc_cds)$joint_annot)) >= 0.05
  )
)
fsk_p3_celltypes

fsk_p3_kc_cds <- (
  fsk_p3_kc_cds[, colData(fsk_p3_kc_cds)$joint_annot %in% fsk_p3_celltypes]
)
colData(fsk_p3_kc_cds) <- droplevels(colData(fsk_p3_kc_cds))

plot_trajectory(
  fsk_p3_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)
plot_trajectory(
  kc_cds_list1[["fsk"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
#trace("calculateLW", edit = T, where = asNamespace("monocle3"))
fsk_p3_deg <- graph_test(
  fsk_p3_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


fsk_p3_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fsk_p3_heatmap <- plot_heatmap(
  fsk_p3_kc_cds,
  fsk_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_fsk.keratinocytes.POSTN_basal_to_matrix_placode.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(fsk_p3_heatmap$plot$gtable)
dev.off()

### organoid ###
# choose one path
get_principal_node(
  kc_cds_list1[["org"]], "joint_annot", "Cuticle/cortex", node_type = "leaf",
  choose_by = "pcw", ascending = FALSE
)
get_principal_node(
  kc_cds_list1[["org"]], "joint_annot", "Cuticle/cortex", node_type = "leaf",
  choose_one = FALSE
)

get_principal_node(
  kc_cds_list1[["org"]], "joint_annot", "Inner root sheath", node_type = "leaf"
)
# POSTN+ basal -> Inner root sheath
org_p1 <- choose_graph_segments(
  kc_cds_list1[["org"]],
  starting_pr_node = "Y_475",
  ending_pr_nodes = c("Y_148", "Y_430"),
  return_list = T
)

org_p1_kc_cds <- (
  (kc_cds_list1[["org"]])[, colnames(kc_cds_list1[["org"]]) %in% org_p1$cells]
)

# take only cell type that makes at leat 5% of the cells of the trajectory
org_p1_celltypes <- names(
  which(
    table(colData(org_p1_kc_cds)$joint_annot) /
      sum(table(colData(org_p1_kc_cds)$joint_annot)) >= 0.05
  )
)
org_p1_celltypes

org_p1_kc_cds <- (
  org_p1_kc_cds[, colData(org_p1_kc_cds)$joint_annot %in% org_p1_celltypes]
)
colData(org_p1_kc_cds) <- droplevels(colData(org_p1_kc_cds))

plot_trajectory(
  org_p1_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)
plot_trajectory(
  kc_cds_list1[["org"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
#trace("calculateLW", edit = T, where = asNamespace("monocle3"))
org_p1_deg <- graph_test(
  org_p1_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


org_p1_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
org_p1_heatmap <- plot_heatmap(
  org_p1_kc_cds,
  org_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_org.keratinocytes.POSTN_basal_to_cuticle_cortex.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(org_p1_heatmap$plot$gtable)
dev.off()

# POSTN+ basal -> Companion layer
org_p2 <- choose_graph_segments(
  kc_cds_list1[["org"]],
  starting_pr_node = "Y_475",
  ending_pr_nodes = "Y_870",
  return_list = T
)

org_p2_kc_cds <- (
  (kc_cds_list1[["org"]])[, colnames(kc_cds_list1[["org"]]) %in% org_p2$cells]
)

# take only cell type that makes at leat 5% of the cells of the trajectory
table(colData(org_p2_kc_cds)$joint_annot)
org_p2_celltypes <- names(
  which(
    table(colData(org_p2_kc_cds)$joint_annot) /
      sum(table(colData(org_p2_kc_cds)$joint_annot)) >= 0.05
  )
)
org_p2_celltypes
org_p2_kc_cds <- (
  org_p2_kc_cds[, colData(org_p2_kc_cds)$joint_annot %in% org_p2_celltypes]
)

colData(org_p2_kc_cds) <- droplevels(colData(org_p2_kc_cds))

plot_trajectory(
  org_p2_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)
plot_trajectory(
  kc_cds_list1[["org"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
#trace("calculateLW", edit = T, where = asNamespace("monocle3"))
org_p2_deg <- graph_test(
  org_p2_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


org_p2_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
org_p2_heatmap <- plot_heatmap(
  org_p2_kc_cds,
  org_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_org.keratinocytes.POSTN_basal_to_companion_layer.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(org_p2_heatmap$plot$gtable)
dev.off()

# POSTN+ basal -> Inner root sheath
org_p3 <- choose_graph_segments(
  kc_cds_list1[["org"]],
  starting_pr_node = "Y_475",
  #ending_pr_nodes = get_principal_node(kc_cds_list1[["org"]], "joint_annot", "Matrix/placode", node_type = "leaf", choose_one = FALSE),
  ending_pr_nodes = c("Y_249", "Y_328"),
  return_list = T
)
get_principal_node(kc_cds_list1[["org"]], "joint_annot", "Inner root sheath", node_type = "leaf", choose_one = FALSE)

org_p3_kc_cds <- (
  (kc_cds_list1[["org"]])[, colnames(kc_cds_list1[["org"]]) %in% org_p3$cells]
)

# take only cell type that makes at leat 5% of the cells of the trajectory
table(colData(org_p3_kc_cds)$joint_annot)
org_p3_celltypes <- names(
  which(
    table(colData(org_p3_kc_cds)$joint_annot) /
      sum(table(colData(org_p3_kc_cds)$joint_annot)) >= 0.01
  )
)
org_p3_celltypes

org_p3_kc_cds <- (
  org_p3_kc_cds[, colData(org_p3_kc_cds)$joint_annot %in% org_p3_celltypes]
)
colData(org_p3_kc_cds) <- droplevels(colData(org_p3_kc_cds))

plot_trajectory(
  org_p3_kc_cds,
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)
plot_trajectory(
  kc_cds_list1[["org"]],
  color_cells_by = "joint_annot",
  legend_loc = "right",
  cell_size = 1
)

# DE along path
#trace("calculateLW", edit = T, where = asNamespace("monocle3"))
org_p3_deg <- graph_test(
  org_p3_kc_cds, neighbor_graph = "principal_graph", cores = 8
)


org_p3_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
  ggplot(aes(x = morans_I, y = -log10(q_value))) +
  geom_point() +
  theme_classic()

source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
org_p3_heatmap <- plot_heatmap(
  org_p3_kc_cds,
  org_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(100),
  #n_cells = 2000,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  rank_threshold = 1
)

pdf(
  file = "pooled_org.keratinocytes.POSTN_basal_to_inner_root_sheath.trajectory_DE.20220315.pdf",
  width = 10,
  height = 10
)
grid.draw(rectGrob(gp = gpar(fill = "white", lwd = 0)))
grid.draw(org_p3_heatmap$plot$gtable)
dev.off()

for (f in c(
  "pooled_fsk.keratinocytes.POSTN_basal_to_inner_root_sheath.trajectory_DE.20220315.pdf",
  "pooled_fsk.keratinocytes.POSTN_basal_to_companion_layer.trajectory_DE.20220315.pdf",
  "pooled_fsk.keratinocytes.POSTN_basal_to_matrix_placode.trajectory_DE.20220315.pdf",
  "pooled_org.keratinocytes.POSTN_basal_to_cuticle_cortex.trajectory_DE.20220315.pdf",
  "pooled_org.keratinocytes.POSTN_basal_to_companion_layer.trajectory_DE.20220315.pdf",
  "pooled_org.keratinocytes.POSTN_basal_to_inner_root_sheath.trajectory_DE.20220315.pdf"
)) {
  cmd <- sprintf('rclone copy --drive-shared-with-me %s "google:/Fetal Skin/Figures/Figs_from_Ni/"', f)
  print(cmd)
  system(cmd)
}

##### END #####
