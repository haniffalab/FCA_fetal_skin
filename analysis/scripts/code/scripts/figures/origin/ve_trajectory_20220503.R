#!/usr/bin/env Rscript

# convert anndata to cds

loadPackage(sceasy)

ve_cds <- sceasy:::anndata2cds(
  "fsk.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220503.h5ad",
  outFile = "fsk.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220503.cds.rds",
  main_layer = "X"
)

ve_cds <- sceasy:::anndata2cds(
  "pooled.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220531.h5ad",
  outFile = "pooled.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220531.cds.rds",
  main_layer = "X"
)

# now run monocle

loadPackage(monocle3)
loadPackage(tidyverse)
loadPackage(viridis)
loadPackage(grid)
loadPackage(ggnewscale)
source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")

run_monocle_cds <- function(
    cds, post_preprocess=T, batch=NULL, residual_model_formula_str=NULL,
    res=0.01, umap_n_neighbors=500L, umap_min_dist=0.3, use_partition=T,
    learn_graph_control = list(minimal_branch_len = 30,
    orthogonal_proj_tip = T), ...
) {
  if (!post_preprocess) {
    cds <- monocle3::preprocess_cds(cds, num_dim = 30)
    if (!is.null(batch)) {
      cds <- monocle3::align_cds(
        cds, alignment_group = batch,
        residual_model_formula_str = residual_model_formula_str
      )
    }
  }
  cds1 <- monocle3::reduce_dimension(
    cds, umap.n_neighbors = umap_n_neighbors, umap.min_dist = umap_min_dist
  )
  cds1 <- monocle3::cluster_cells(cds1, resolution = res)
  cds1 <- monocle3::learn_graph(
    cds1, use_partition = use_partition,
    learn_graph_control = learn_graph_control, ...
  )
  cds1
}

# 20220503 analysis
ve_cds <- readRDS(
  "fsk.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220503.cds.rds"
)
ve_cds <- readRDS(
  "pooled.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220531.cds.rds"
)
str(ve_cds)

colnames(colData(ve_cds))
colData(ve_cds)$joint_annot <- colData(ve_cds)$joint_annotation_20220202
table(colData(ve_cds)$joint_annot)
table(colData(ve_cds)$pcw)
levels(colData(ve_cds)$week)
levels(colData(ve_cds)$day)

colData(ve_cds)$pcw <- plyr::mapvalues(
  colData(ve_cds)$week,
  c(
    "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
    "4-7_fetal_wks",
    "7-10_fetal_wks",
    "14-16_fetal_wks",
    "17-20_fetal_wks"
  ),
  round(c(c(7:16)*7, 29, 48, 85, 133) / 7)
)
colData(ve_cds)$pcw <- as.integer(as.character(colData(ve_cds)$pcw))

ve_cds_list1 <- list(
  # Start from pre-computed PC
  #run1 = run_monocle_cds(
  #  ve_cds,
  #  res = 0.01,
  #  use_partition = T,
  #  umap_n_neighbors = 50L,
  #  umap_min_dist = 1,
  #  learn_graph_control = list(minimal_branch_len = 10, orthogonal_proj_tip = T)
  #),
  run2 = run_monocle_cds(
    ve_cds,
    res = 0.005,
    use_partition = F,
    umap_n_neighbors = 50L,
    umap_min_dist = 1,
    learn_graph_control = list(minimal_branch_len = 12, orthogonal_proj_tip = T)
  ),
  # Start from scratch 1
  #run3 = run_monocle_cds(
  #  ve_cds,
  #  post_preprocess = F,
  #  batch = "chemistry",
  #  res = 0.01,
  #  residual_model_formula_str = "~percent_mito + nCounts_RNA + cc_score",
  #  use_partition = T,
  #  umap_n_neighbors = 30L,
  #  umap_min_dist = 0.2,
  #  learn_graph_control = list(minimal_branch_len = 15, orthogonal_proj_tip = T)
  #),
  # Start from scratch 2
  run4 = run_monocle_cds(
    ve_cds,
    post_preprocess = F,
    batch = "sanger_id",
    res = 0.01,
    #residual_model_formula_str = "~percent_mito + nCounts_RNA + cc_score",
    use_partition = T,
    umap_n_neighbors = 10L,
    umap_min_dist = 0.1,
    learn_graph_control = list(minimal_branch_len = 10, orthogonal_proj_tip = T)
  )
)
ve_cds_list1$run2 <- run_monocle_cds(
  ve_cds,
  res = 0.01,
  use_partition = F,
  umap_n_neighbors = 50L,
  umap_min_dist = 0.8,
  learn_graph_control = list(minimal_branch_len = 12, orthogonal_proj_tip = T),
  close_loop = T
)

# plot to help choose root and path
pdf(file = "fsk.ve_trajectory_helper.pdf", width = 9, height = 9)
#pdf(file = "pooled.ve_trajectory_helper.pdf", width = 9, height = 9)
for (name in names(ve_cds_list1)) {
  print(
    plot_trajectory(
      ve_cds_list1[[name]],
      color_cells_by = "joint_annot",
      legend_loc = "right",
      label_graph_nodes = "all",
      graph_label_size = 2.5
    )
    #+ ggtitle(paste("pooled", name, sep = ", "))
    + ggtitle(paste("fsk", name, sep = ", "))
  )
  #print(
  #  plot_trajectory(
  #    ve_cds_list1[[name]],
  #    color_cells_by = "dataset",
  #    legend_loc = "right",
  #    label_graph_nodes = "all",
  #    graph_label_size = 2.5
  #  )
  #  + ggtitle(paste("pooled", name, sep = ", "))
  #)
  print(
    plot_trajectory(
      ve_cds_list1[[name]],
      color_cells_by = "pcw",
      legend_loc = "right",
      label_graph_nodes = "none",
      graph_label_size = 2.5
    )
    #+ ggtitle(paste("pooled", name, sep = ", "))
    + ggtitle(paste("fsk", name, sep = ", "))
  )
}
dev.off()

root_celltypes <- list(
  #run1 = c('Early endothelial cell', 'Tip cell (arterial)'),
  run2 = c('Early endothelial cell')
  #run3 = c('Early endothelial cell'),
  #run4 = c('Early endothelial cell')
)

# Confirming that selected nodes look good
for (name in names(root_celltypes)) {
  for (root_celltype in root_celltypes[[name]]) {
    for (node_type in c("any", "leaf", "branch")) {
      node <- get_principal_node(ve_cds_list1[[name]], 'joint_annot', root_celltype, node_type=node_type)
      if (is.null(node)) node <- "none"
      cat(sprintf("%s, %s, %s: %s\n", name, root_celltype, node_type, node))
    }
  }
}

# Calculate pseudotime
for (name in names(root_celltypes)) {
  ve_cds_list1[[name]] <- order_cells(
    ve_cds_list1[[name]],
    root_pr_nodes=sapply(
      root_celltypes[[name]],
      function(root_celltype) get_principal_node(ve_cds_list1[[name]], 'joint_annot', root_celltype, node_type="leaf")
    )
  )
}

#saveRDS(ve_cds_list1[["run1"]], "fsk.vascular_endothelium.monocle.run1.rds")
#saveRDS(ve_cds_list1[["run2"]], "fsk.vascular_endothelium.monocle.run2.rds")
#saveRDS(ve_cds_list1[["run3"]], "fsk.vascular_endothelium.monocle.run3.rds")
#saveRDS(ve_cds_list1[["run4"]], "fsk.vascular_endothelium.monocle.run4.rds")

saveRDS(ve_cds_list1[["run2"]], "pooled.vascular_endothelium.monocle.run2.rds")

ve_cds_list1 = list(
  #run1=readRDS("fsk.vascular_endothelium.monocle.run1.rds"),
  #run2=readRDS("fsk.vascular_endothelium.monocle.run2.rds"),
  #run3=readRDS("fsk.vascular_endothelium.monocle.run3.rds"),
  #run4=readRDS("fsk.vascular_endothelium.monocle.run4.rds")
  run2=readRDS("pooled.vascular_endothelium.monocle.run2.rds")
)

colorby_variables <- list(
  #run1=list(discrete=c("joint_annot", "dataset"), continuous=c("pcw", "pseudotime")),
  run2=list(discrete=c("dataset"), continuous=c("pcw", "pseudotime"))
  #run3=list(discrete=c("joint_annot", "dataset"), continuous=c("pcw", "pseudotime")),
  #run4=list(discrete=c("joint_annot", "dataset"), continuous=c("pcw", "pseudotime"))
)

ve_color_codes <- c(
  `Early endothelial cell` = "#2078B5",
  `Early LE` = "#F07E21",
  `LE` = "#2CA137",
  `Tip cell (arterial)` = "#8F67AA",
  `Arterial` = "#D72829",
  `Capillary (venular tip)` = "#8D574C",
  `Capillary/postcapillary venule` = "#D37AB0",
  `Postcapillary venule` = "#7F8080"
)

# choose one run
name <- "run2"

cds <- ve_cds_list1[[name]]
colData(cds)[, c("pcw", "dataset")] %>% cbind(reducedDims(cds)[["UMAP"]]) %>% head()
cds_df <- (
  colData(cds)[, c("pcw", "dataset")] %>%
    cbind(reducedDims(cds)[["UMAP"]]) %>%
    as_tibble() %>%
    mutate(UMAP1=V1, UMAP2=V2, .keep="unused")
)
node_df <- (
  cds@principal_graph_aux[["UMAP"]]$dp_mst %>%
    t() %>%
    as_tibble(rownames = "node_id", .repair_name = "minimal") %>%
    mutate(UMAP1 = V1, UMAP2 = V2, .keep="unused")
)

edge_df <- (
  cds@principal_graph[["UMAP"]] %>%
    igraph::as_data_frame() %>%
    as_tibble() %>%
    select(from, to) %>%
    left_join(node_df, by=c(from="node_id")) %>%
    right_join(node_df, by=c(to="node_id")) %>%
    select(from_x=`UMAP1.x`, from_y=`UMAP2.x`, to_x=`UMAP1.y`, to_y=`UMAP2.y`)
)

#p0_pdf <- "fsk.vascular_endothelium.monocle.20220527.pdf"
p0_pdf <- "pooled.vascular_endothelium.monocle.20220607.pdf"
pdf(file=p0_pdf, width=8, height=6)
for (v in colorby_variables[[name]][["discrete"]]) {
  print(plot_trajectory(ve_cds_list1[[name]], color_cells_by=v, legend_loc="right") + ggtitle(name))
}
plot_trajectory(ve_cds_list1[[name]], color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root", color_palette = color_codes)
for (v in colorby_variables[[name]][["continuous"]]) {
  print(plot_trajectory(ve_cds_list1[[name]], color_cells_by=v, legend_loc="right") + ggtitle(name))
}
(
  ggplot() +
    geom_point(
      data = cds_df %>% mutate(dataset = recode_factor(dataset, fsk = "Fetal", org = "Orgnoid")),
      aes(x = UMAP1, y = UMAP2, fill = pcw, color = dataset), shape = 21, size = 2
    ) +
    scale_fill_viridis(limits = c(4, 19)) +
    scale_color_manual(values=c(NA, "red")) +
    geom_segment(data=edge_df, aes(x = from_x, xend = to_x, y = from_y, yend = to_y), size = 1, alpha=0.75) +
    theme_classic()
)

(
  cds_df %>%
    filter(dataset == "fsk") %>%
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = pcw)) +
    scale_color_distiller(palette = "Blues", limits = c(4, 19), direction = 1, name = "Fetal pcw") +
    new_scale_color() +
    geom_point(data=cds_df %>% filter(dataset == "org"), aes(x = UMAP1, y = UMAP2, color = pcw)) +
    scale_color_distiller(palette = "Reds", limits = c(4, 19), direction = 1, name = "Organoid pcw") +
    geom_segment(data=edge_df, aes(x = from_x, xend = to_x, y = from_y, yend = to_y), size = 1) +
    theme_classic()
)

dev.off()

system(paste("rclone copy --drive-shared-with-me", p0_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

plot_trajectory(ve_cds_list1[[name]], color_cells_by="joint_annot", legend_loc="right", label_graph_nodes=T)

# choose one path: Early endothelial cell -> Arterial
fsk_p1 <- choose_graph_segments(
  ve_cds_list1[[name]],
  starting_pr_node = get_principal_node(ve_cds_list1[[name]], "joint_annot", root_celltypes[[name]], node_type="leaf"),
  ending_pr_nodes = get_principal_node(ve_cds_list1[[name]], "joint_annot", "Arterial", node_type="leaf"),
  return_list=T
)

fsk_p1_ve_cds <- (ve_cds_list1[[name]])[, colnames(ve_cds_list1[[name]]) %in% fsk_p1$cells]
table(colData(fsk_p1_ve_cds)$joint_annot)/sum(table(colData(fsk_p1_ve_cds)$joint_annot))
(fsk_p1_celltypes <- names(which(table(colData(fsk_p1_ve_cds)$joint_annot)/sum(table(colData(fsk_p1_ve_cds)$joint_annot)) >= 0.05)))

fsk_p1_ve_cds <- fsk_p1_ve_cds[, colData(fsk_p1_ve_cds)$joint_annot %in% fsk_p1_celltypes]
#colData(fsk_p1_ve_cds) <- droplevels(colData(fsk_p1_ve_cds))

colData(fsk_p1_ve_cds)$joint_annot <- factor(
  colData(fsk_p1_ve_cds)$joint_annot,
  levels = c(
    "Early endothelial cell", "Tip cell (arterial)", "Arterial",
    "Capillary (venular tip)", "Capillary/postcapillary venule", "Postcapillary venule"
  )
)

# DE along path
fsk_p1_deg <- graph_test(fsk_p1_ve_cds, neighbor_graph = "principal_graph", cores = 8)

fsk_p1_heatmap1 <- plot_heatmap(
  fsk_p1_ve_cds,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 10),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(
    joint_annot = ve_color_codes[
      names(ve_color_codes) %in% unique(colData(fsk_p1_ve_cds)$joint_annot)
    ],
    pcw = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Blues"))(100)[
      1:round(100*16/18)
    ]
  )
)

fsk_p1_heatmap2 <- plot_heatmap(
  fsk_p1_ve_cds,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = c(
    "HPGD", "ITM2A", "SLC7A5", "TMEM176B",
    "SLC38A5", "APLNR", "HAPLN3", "MFSD2A", "LGALS1",
    "MMP1", "PRSS23", "PGF", "PRSS2", "GYPC",
    "HLA-B", "HLA-C", "COL4A1", "SPARCL1", "PRND",
    "CAV1", "IGF2", "AQP1", "IGFBP3", "GJA5",
    "FBLN5", "CXCL12", "TMEM100", "SRGN"
  ),
  order_genes_by = "rank",
  annotation_colors = list(
    joint_annot = ve_color_codes[
      names(ve_color_codes) %in% unique(colData(fsk_p1_ve_cds)$joint_annot)
    ],
    pcw = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Blues"))(100)[
      1:round(100*16/18)
    ]
  ),
  drop_levels = F
)

fsk_p1_heatmap3 <- plot_heatmap(
  fsk_p1_ve_cds,
  fsk_p1_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = c(
    "HPGD", "ITM2A", "SLC7A5", "TMEM176B",
    "SLC38A5", "APLNR", "HAPLN3", "MFSD2A", "LGALS1",
    "MMP1", "PRSS23", "PGF", "PRSS2", "GYPC",
    "HLA-B", "HLA-C", "COL4A1", "SPARCL1", "PRND",
    "CAV1", "IGF2", "AQP1", "IGFBP3", "GJA5",
    "FBLN5", "CXCL12", "TMEM100", "SRGN"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ve_color_codes[
    names(ve_color_codes) %in% unique(colData(fsk_p1_ve_cds)$joint_annot)
  ])
)

p1_pdf <- "fsk.vascular_endothelium.monocle.early_endothelial_to_arterial.trajectory_DE.20220711.pdf"
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
  fsk_p1_ve_cds, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root",
  color_palette = ve_color_codes[names(ve_color_codes) %in% unique(colData(fsk_p1_ve_cds)$joint_annot)],
  cell_size=1
)
dev.off()
system(paste("rclone copy --drive-shared-with-me", p1_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# choose one path: Early endothelial cell -> Postcapillary venule
fsk_p2 <- choose_graph_segments(
  ve_cds_list1[[name]],
  starting_pr_node = get_principal_node(ve_cds_list1[[name]], "joint_annot", root_celltypes[[name]], node_type="leaf"),
  ending_pr_nodes = get_principal_node(ve_cds_list1[[name]], "joint_annot", "Postcapillary venule", node_type="leaf"),
  return_list=T
)

fsk_p2_ve_cds <- (ve_cds_list1[[name]])[, colnames(ve_cds_list1[[name]]) %in% fsk_p2$cells]
table(colData(fsk_p2_ve_cds)$joint_annot)/sum(table(colData(fsk_p2_ve_cds)$joint_annot))
(fsk_p2_celltypes <- names(which(table(colData(fsk_p2_ve_cds)$joint_annot)/sum(table(colData(fsk_p2_ve_cds)$joint_annot)) >= 0.05)))

fsk_p2_ve_cds <- fsk_p2_ve_cds[, colData(fsk_p2_ve_cds)$joint_annot %in% fsk_p2_celltypes]
colData(fsk_p2_ve_cds) <- droplevels(colData(fsk_p2_ve_cds))

colData(fsk_p2_ve_cds)$joint_annot <- factor(
  colData(fsk_p2_ve_cds)$joint_annot,
  levels=c("Early endothelial cell", "Tip cell (arterial)", "Capillary (venular tip)", "Capillary/postcapillary venule", "Postcapillary venule")
)

# DE along path
fsk_p2_deg <- graph_test(fsk_p2_ve_cds, neighbor_graph = "principal_graph", cores = 8)

fsk_p2_heatmap1 <- plot_heatmap(
  fsk_p2_ve_cds,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 7),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ve_color_codes[
    names(ve_color_codes) %in% unique(colData(fsk_p2_ve_cds)$joint_annot)
  ])
)

fsk_p2_heatmap2 <- plot_heatmap(
  fsk_p2_ve_cds,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = c(
    "ITM2A", "SLC7A5", "HPGD", "SYT1", "CDH2", "SLC38A5", "ECM1",
    "LGALS1", "ESM1", "CXCR4", "PDLIM2", "PGF", "PRSS23", "NOV",
    "COL4A1", "COL14A1", "IGF2", "FABP4", "IL32", "IGFBP7",
    "CAV1", "SPARC", "BCAM", "IL33", "LGALS3", "VWF", "SELE",
    "ACKR1", "PLVAP", "CPE", "FBLN2", "CLU", "MMRN1", "VCAN", "ELN"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ve_color_codes[
    names(ve_color_codes) %in% unique(colData(fsk_p2_ve_cds)$joint_annot)
  ])
)

fsk_p2_heatmap3 <- plot_heatmap(
  fsk_p2_ve_cds,
  fsk_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(500),
  n_cells = 1000,
  n_genes_per_level = list(joint_annot = 100),
  rank_threshold = 0.95,
  min_gene_sd = 0.2,
  min_gene_fraction = 0.1,
  cell_annot = c("joint_annot", "pcw"),
  show_rownames = c(
    "ITM2A", "SLC7A5", "HPGD", "SYT1", "CDH2", "SLC38A5", "ECM1",
    "LGALS1", "ESM1", "CXCR4", "PDLIM2", "PGF", "PRSS23", "NOV",
    "COL4A1", "COL14A1", "IGF2", "FABP4", "IL32", "IGFBP7",
    "CAV1", "SPARC", "BCAM", "IL33", "LGALS3", "VWF", "SELE",
    "ACKR1", "PLVAP", "CPE", "FBLN2", "CLU", "MMRN1", "VCAN", "ELN"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annot = ve_color_codes[
    names(ve_color_codes) %in% unique(colData(fsk_p2_ve_cds)$joint_annot)
  ])
)

p2_pdf <- "fsk.vascular_endothelium.monocle.early_endothelial_to_postcapillary_venule.trajectory_DE.20220711.pdf"
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
  fsk_p2_ve_cds, color_cells_by = "joint_annot", legend_loc = "right", label_graph_nodes = "root",
  color_palette = ve_color_codes[names(ve_color_codes) %in% unique(colData(fsk_p2_ve_cds)$joint_annot)],
  cell_size=1
)
dev.off()
system(paste("rclone copy --drive-shared-with-me", p2_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# === END ===
rm(list=ls())
gc()
