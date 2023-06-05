#!/usr/bin/env Rscript

# convert anndata to cds
#conda activate sceasy

library(sceasy)

cds <- sceasy:::anndata2cds(
  "pooled_adult_fetal_skin.counts.h5ad",
  outFile = "pooled_adult_fetal_skin.counts.cds",
  main_layer = "X"
)

# now run monocle
#conda activate monocle3
library(monocle3)
library(reticulate)
library(grid)
library(tidyverse)
source("~/sctkr/R/MonocleUtils.R")

data_root <- '/lustre/scratch126/cellgen/team205/nh3/skin'
proj_root <- '~/FCA_Fetal_Skin_priv'

trace("calculateLW", edit = T, where = asNamespace("monocle3"))

get_earliest_principal_node <- function(cds, group_by, root_group) {
  cell_ids <- which(colData(cds)[, group_by] == root_group)

  closest_vertex <- (
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  )
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]

  root_pr_nodes
}

run_monocle_cds <- function(cds, batch=NA, add_pca=NA, add_umap=NA, res=1e-2,
                            use_partition=F, learn_graph_control=NULL) {
  cds@colData <- droplevels(cds@colData)
  cds1 <- monocle3::preprocess_cds(cds, num_dim = 50)
  if (!is.na(batch)) {
    cds1 <- monocle3::align_cds(cds1, num_dim = 50, alignment_group = batch)
  }
  if (!is.na(add_pca) || !is.na(add_umap)) {
    embeds <- S4Vectors::SimpleList()
    if (!is.na(add_pca)) embeds$PCA <- add_pca
    if (!is.na(add_umap)) embeds$UMAP <- add_umap
    SingleCellExperiment::reducedDims(cds1) <- embeds
  }
  if (is.na(add_umap)) {
    cds1 <- monocle3::reduce_dimension(cds1)
  }
  cds1 <- monocle3::cluster_cells(cds1, resolution = res)
  cds1 <- monocle3::learn_graph(
    cds1, use_partition = use_partition,
    learn_graph_control = learn_graph_control
  )
  cds1
}

kc_color_codes <- c(
  `POSTN+ basal` = "#7C4F9A",
  `Basal` = "#D72029",
  `Matrix/placode` = "#CEB08D",
  `Cuticle/cortex` = "#8C584C",
  `Inner root sheath` = "#7E7E7F",
  `Outer root sheath` = "#D179AF",
  `Companion layer` = "#BABC26"
)
color_palette <- sapply(fsk_ms2_celltypes, function(ct) color_codes[[ct]])


# read data
cds <- readRDS(file.path(data_root, '20211022_final_figures', "pooled_adult_fetal_skin.counts.cds.rds"))
fetal_cds <- cds[, cds@colData@listData$dataset == "fetal"]
fetal_cds@colData <- droplevels(fetal_cds@colData[, c(1:27, 41:43)])
rm(cds)

# Keratinocytes
kc_ad <- anndata$read_h5ad(
  file.path(data_root, "20200626_make_figure_for_Muzz/pooled_keratinocytes.processed.h5ad"))
kc_cell_ids <- as.vector(sub("_skin", "", kc_ad$obs_names$values))

k_kc1 <- (
  (colData(fetal_cds)$independent_annotation_broad3 == "KC")
  & !colData(fetal_cds)$joint_annotation %in% c(
    "Periderm", "Immature basal", "Immature suprabasal", "Suprabasal IFE"
  )
  & (colData(fetal_cds)$pcw > 12)
)
fetal_kc_cds0 <- fetal_cds[, k_kc1]

fetal_kc_cds1 <- fetal_kc_cds0[, colnames(fetal_kc_cds0) %in% kc_cell_ids]
fetal_kc_cds1@colData <- droplevels(fetal_kc_cds1@colData)
k_cds1 <- kc_cell_ids %in% colnames(fetal_kc_cds1)
kc_X_pca_hm <- kc_ad$obsm["X_pca_hm"][k_cds1, ]

fetal_kc_cds1 <- run_monocle_cds(
  fetal_kc_cds1, add_pca = kc_X_pca_hm, res = 0.01,
  learn_graph_control = list(minimal_branch_len = 5)
)

levels(colData(fetal_kc_cds1)$joint_annotation)[levels(colData(fetal_kc_cds1)$joint_annotation) == "Basal SHH+"] <- "Matrix/placode"
levels(colData(fetal_kc_cds1)$joint_annotation)[levels(colData(fetal_kc_cds1)$joint_annotation) == "Basal POSTN+"] <- "POSTN+ basal"

fetal_kc_cds1 <- order_cells(fetal_kc_cds1, root_pr_nodes=get_earliest_principal_node(fetal_kc_cds1, 'joint_annotation', 'POSTN+ basal'))

# saveRDS(fetal_kc_cds1, "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.monocle.cds.rds")
# fetal_kc_cds1 <- readRDS(file.path(data_root, '20211022_final_figures', 'monocle', "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.monocle.cds.rds"))

pdf(file=file.path(proj_root, "figures", "final", "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.monocle3.20220329.pdf"), width=7, height=6)
plot_trajectory(fetal_kc_cds1, color_cells_by="joint_annotation", cell_size=1, legend_loc="right")
plot_trajectory(fetal_kc_cds1, color_cells_by="joint_annotation", cell_size=1, legend_loc="right", color_palette = color_codes)
plot_trajectory(fetal_kc_cds1, color_cells_by="pseudotime", cell_size=1)
plot_trajectory(fetal_kc_cds1, color_cells_by="pcw", cell_size=1)
plot_trajectory(fetal_kc_cds1, color_cells_by="joint_annotation", cell_size=1, label_graph_nodes="all", legend_loc=F)
dev.off()
# system("rclone copy --drive-shared-with-me pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.monocle3.20220329.pdf 'google:/Fetal Skin/Figures/Figs_from_Ni/'")

plot_cells(
  fetal_kc_cds1,
  genes=c(
    "KRT85", "SHH", "SPON2", "ICA1", "SCG5", "SOSTDC1", "AGR2", "DAPL1", "SOX9",
    "NFATC1", "ASS1", "KRT79", "DKK4", "LHX2", "KRT75", "INHBA", "HSD17B2",
    "SLC26A2", "LMCD1", "PAQR5", "CALML5", "GRHL1", "LYPD3", "RHOV", "FRMD4A"),
  norm_method="log",
  min_expr=0.5,
  scale_to_range=T,
  show_trajectory_graph=F, label_cell_groups=F, cell_size=0.5
) + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())

# Basal POSTN+ -> Matrix/placode
fetal_p1 <- choose_graph_segments(
  fetal_kc_cds1,
  starting_pr_node = "Y_279",
  ending_pr_nodes = c("Y_249", "Y_233", "Y_166", "Y_105"),
  return_list = T
)

fetal_p1_kc_cds <- fetal_kc_cds1[, colnames(fetal_kc_cds1) %in% fetal_p1$cells]

table(colData(fetal_p1_kc_cds)$joint_annotation)/sum(table(colData(fetal_p1_kc_cds)$joint_annotation))

fetal_p1_celltypes <- names(which(table(colData(fetal_p1_kc_cds)$joint_annotation)/sum(table(colData(fetal_p1_kc_cds)$joint_annotation)) >= 0.05))
fetal_p1_kc_cds1 <- fetal_p1_kc_cds[, colData(fetal_p1_kc_cds)$joint_annotation %in% fetal_p1_celltypes]
colData(fetal_p1_kc_cds1) <- droplevels(colData(fetal_p1_kc_cds1))

plot_trajectory(fetal_p1_kc_cds, color_cells_by="joint_annotation", cell_size=1, label_graph_nodes=T)

fetal_p1_deg <- graph_test(fetal_p1_kc_cds1, neighbor_graph = "principal_graph", cores = 26)

fetal_p1_deg %>% filter(status=="OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

colData(fetal_p1_kc_cds)$joint_annotation <- factor(
  colData(fetal_p1_kc_cds)$joint_annotation, levels=c(
    "POSTN+ basal", "Basal", "Matrix/placode", "Cuticle/cortex"
  )
)

custom_kc_genes = c(
  "IGFBP6", "SERPINF1", "FABP5", "SNCG", "IGFBP5", "IGF2", "SERPINB7", "ATP5O",
  "MEIS2", "CRABP2", "DAPL1", "RARRES2", "CREB5", "CARD16", "SPON2", "LRRCC1",
  "TNC", "SEMA3B", "BDNF", "CASP1", "ODC1", "WIPI1", "IL11RA", "SCUBE3",
  "LAMA3", "GALT", "PRKCG", "LTB", "DPEP1", "NOTUM", "HSPA6", "LYPD3", "SGF29",
  "MX1", "RFX5", "KIF9", "ZNF549", "OPHN1", "USPL1", "ZFYVE9"
)

fetal_p1_heatmap_a <- plot_heatmap(
  fetal_p1_kc_cds,
  custom_kc_genes,
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  scale_by_gene = F,
  show_rownames = T,
  order_genes_by = F
)

fetal_p1_heatmap_b <- plot_heatmap(
  fetal_p1_kc_cds,
  fetal_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(60),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  show_rownames = T,
  order_genes_by = "rank"
)

fetal_p1_heatmap_c <- plot_heatmap(
  fetal_p1_kc_cds1,
  fetal_p1_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% arrange(q_value, -morans_I) %>% rownames() %>% head(500),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  n_genes_per_level = list(joint_annotation = 15),
  show_rownames = T,
  order_genes_by = "rank"
)

pdf(file=file.path(proj_root, "figures", "obsolete", "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.basal_postn_to_matrix_placode.monocle3.trajectory_DE.20220314.pdf"), width=10, height=6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p1_heatmap_a$plot$gtable)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p1_heatmap_b$plot$gtable)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p1_heatmap_c$plot$gtable)
dev.off()

# Basal POSTN+ -> Outer root sheath/Companion layer
fetal_p2 <- choose_graph_segments(
  fetal_kc_cds1,
  starting_pr_node = "Y_279",
  ending_pr_nodes = c("Y_44"),
  return_list = T
)

fetal_p2_kc_cds <- fetal_kc_cds1[, colnames(fetal_kc_cds1) %in% fetal_p2$cells]

table(colData(fetal_p2_kc_cds)$joint_annotation)/sum(table(colData(fetal_p2_kc_cds)$joint_annotation))
fetal_p2_kc_cds <- fetal_p2_kc_cds[, colData(fetal_p2_kc_cds)$joint_annotation != "Inner root sheath"]
colData(fetal_p2_kc_cds) <- droplevels(colData(fetal_p2_kc_cds))

fetal_p2_deg <- graph_test(fetal_p2_kc_cds, neighbor_graph = "principal_graph", cores = 26)
fetal_p2_deg$gene_ids.fetal <- NULL
fetal_p2_deg$gene_ids.SKN8090524.adult <- NULL
fetal_p2_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% head()

fetal_p2_heatmap1 <- plot_heatmap(
  fetal_p2_kc_cds,
  fetal_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  show_rownames = T,
  min_gene_fraction = 0.15,
  rank_threshold = 0.95,
  n_genes_per_level = list(joint_annotation = 8),
  order_genes_by = "rank",
  smooth_heatmap = 3,
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p2_kc_cds)$joint_annotation)
  ])
)

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fetal_p2_heatmap2 <- plot_heatmap(
  fetal_p2_kc_cds,
  fetal_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  min_gene_fraction = 0.15,
  rank_threshold = 0.95,
  n_genes_per_level = list(joint_annotation = 100),
  smooth_heatmap = 3,
  show_rownames = c(
    "IGFBP6", "ABRACL", "KRT15", "NPR3", "POSTN",
    "IGF2", "WNT3", "DCN", "VCAN", "COL14A1",
    "GJB6", "CALB2", "DAPL1", "TM4SF1", "FRZB",
    "SPON2", "CREB5", "TINCR", "CRES4", "DHCR7",
    "HLA-C", "DPP4", "DSG3", "PALMD", "DMKN",
    "LYPD3", "AQP3", "SULT2B1", "CLDN1", "ALCAM",
    "DSP", "NUPR1", "DBI", "CLDN4"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p2_kc_cds)$joint_annotation)
  ])
)

fetal_p2_heatmap3 <- plot_heatmap(
  fetal_p2_kc_cds,
  fetal_p2_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  min_gene_fraction = 0.15,
  rank_threshold = 0.95,
  n_genes_per_level = list(joint_annotation = 100),
  smooth_heatmap = 3,
  show_rownames = c(
    "IGFBP6", "ABRACL", "KRT15", "NPR3", "POSTN",
    "IGF2", "WNT3", "DCN", "VCAN", "COL14A1",
    "GJB6", "CALB2", "DAPL1", "TM4SF1", "FRZB",
    "SPON2", "CREB5", "TINCR", "CRES4", "DHCR7",
    "HLA-C", "DPP4", "DSG3", "PALMD", "DMKN",
    "LYPD3", "AQP3", "SULT2B1", "CLDN1", "ALCAM",
    "DSP", "NUPR1", "DBI", "CLDN4"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p2_kc_cds)$joint_annotation)
  ])
)

p2_pdf <- file.path(proj_root, "figures", "final", "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.basal_postn_to_companion_layer.monocle3.trajectory_DE.20220711.pdf")
pdf(file = p2_pdf, width = 8, height = 6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p2_heatmap1$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p2_heatmap2$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p2_heatmap3$plot)
plot_trajectory(fetal_p2_kc_cds, color_cells_by = "joint_annotation", cell_size = 1)
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p2_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

# Basal POSTN+ -> Inner root sheath
fetal_p3 <- choose_graph_segments(
  fetal_kc_cds1,
  starting_pr_node = "Y_279",
  ending_pr_nodes = c("Y_29"),
  return_list = T
)

fetal_p3_kc_cds <- fetal_kc_cds1[, colnames(fetal_kc_cds1) %in% fetal_p3$cells]

table(colData(fetal_p3_kc_cds)$joint_annotation) / sum(table(colData(fetal_p3_kc_cds)$joint_annotation))
fetal_p3_kc_cds <- fetal_p3_kc_cds[, colData(fetal_p3_kc_cds)$joint_annotation != "Cuticle/cortex"]
colData(fetal_p3_kc_cds) <- droplevels(colData(fetal_p3_kc_cds))

fetal_p3_deg <- graph_test(fetal_p3_kc_cds, neighbor_graph = "principal_graph", cores = 26)
fetal_p3_deg1 <- fetal_p3_deg[complete.cases(fetal_p3_deg), ]

fetal_p3_deg %>% filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01) %>% ggplot(aes(x=morans_I, y=-log10(q_value))) + geom_point() + theme_classic()

fetal_p3_deg %>%
  filter(status == "OK", !mito, !ribo, !hb, q_value < 0.01, !grepl("[.]", rownames(.))) %>%
  arrange(q_value, -morans_I) %>%
  mutate(rank=seq_along(q_value)) %>%
  head()

fetal_p3_heatmap1 <- plot_heatmap(
  fetal_p3_kc_cds,
  fetal_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  n_genes_per_level = list(joint_annotation = 10),
  min_gene_fraction = 0.05,
  min_gene_sd = 0.15,
  rank_threshold = 0.95,
  smooth_heatmap = 3,
  show_rownames = T,
  order_genes_by = "rank",
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p3_kc_cds)$joint_annotation)
  ])
)

# source("~/src/github/Teichlab/sctkr/R/MonocleUtils.R")
fetal_p3_heatmap2 <- plot_heatmap(
  fetal_p3_kc_cds,
  fetal_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  n_genes_per_level = list(joint_annotation = 100),
  min_gene_fraction = 0.05,
  min_gene_sd = 0.2,
  rank_threshold = 0.95,
  smooth_heatmap = 3,
  show_rownames = c(
    "IGFBP6", "KRT14", "DLK1", "ABRACL", "KRT15",
    "CRABP2", "POSTN", "SOX6", "IGF2", "IL32",
    "VCAN", "TNC", "DAPL1", "TM4SF1", "CLDN1",
    "RARRES2", "WIF1", "FRZB", "AQP3",
    "NRP2", "CXADR", "SAT1", "GJA1", "DSP",
    "PVALB", "PVRL4", "APOE", "SSFA2", "DSC3"
  ),
  order_genes_by = "rank",
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p3_kc_cds)$joint_annotation)
  ])
)

fetal_p3_heatmap3 <- plot_heatmap(
  fetal_p3_kc_cds,
  fetal_p3_deg %>%
    filter(status == "OK", !mito, !ribo, !hb, q_value < 0.05) %>%
    arrange(q_value, -morans_I) %>%
    rownames() %>%
    head(3000),
  n_cells = 2000,
  cell_annot = c("joint_annotation", "pcw"),
  n_genes_per_level = list(joint_annotation = 100),
  min_gene_fraction = 0.05,
  min_gene_sd = 0.2,
  rank_threshold = 0.95,
  smooth_heatmap = 3,
  show_rownames = c(
    "IGFBP6", "KRT14", "DLK1", "ABRACL", "KRT15",
    "CRABP2", "POSTN", "SOX6", "IGF2", "IL32",
    "VCAN", "TNC", "DAPL1", "TM4SF1", "CLDN1",
    "RARRES2", "WIF1", "FRZB", "AQP3",
    "NRP2", "CXADR", "SAT1", "GJA1", "DSP",
    "PVALB", "PVRL4", "APOE", "SSFA2", "DSC3"
  ),
  show_arrows = TRUE,
  order_genes_by = "rank",
  annotation_colors = list(joint_annotation = kc_color_codes[
    names(kc_color_codes) %in% unique(colData(fetal_p3_kc_cds)$joint_annotation)
  ])
)
# maybe the following should be in "final"
p3_pdf <- file.path(proj_root, "figures", "final", "pooled_keratinocytes.processed.fetal_hair_folicle_pcw12.basal_postn_to_inner_root_sheath.monocle3.trajectory_DE.20220711.pdf")
pdf(file = p3_pdf, width = 8, height = 6)
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p3_heatmap1$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p3_heatmap2$plot)
grid.newpage()
grid.draw(rectGrob(gp=gpar(fill = "white", lwd = 0)))
grid.draw(fetal_p3_heatmap3$plot)
plot_trajectory(fetal_p3_kc_cds, color_cells_by = "joint_annotation", cell_size = 1)
dev.off()
# system(paste("rclone copy --drive-shared-with-me", p3_pdf, "'google:/Fetal Skin/Figures/Figs_from_Ni/'"))

##### END #####
# Sys.setenv(DISPLAY="localhost:11.0")

# fetal_p3_heatmap2$plot$layout$name
# ?pheatmap::pheatmap
