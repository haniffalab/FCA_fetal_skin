#!/usr/bin/env Rscript

#####################
# convert anndata to sce
# conda activate zellkonverter
# 
# loadPackage(zellkonverter)
# loadPackage(reticulate)
# reticulate::use_condaenv('zellkonverter')
# 
# anndata <- reticulate::import('anndata')
# ad <- anndata$read_h5ad('fetal_skin.scvi_donor.norm.h5ad')
# ad
# ad$raw
# sce <- AnnData2SCE(ad, X_name='logcounts', layers=F, uns=T, var=T, obs=T, varm=F, obsm=T, obsp=T, hdf5_backed=F)
# 
# saveRDS(sce, 'fetal_skin.scvi_donor.norm.sce.rds')


#####################
# now run milo
# conda activate milo

loadPackage(miloR)
loadPackage(SingleCellExperiment)
loadPackage(scater)
loadPackage(dplyr)
loadPackage(Matrix)


# new analysis 20220202
sce <- readRDS('fetal_skin.scvi_donor.norm.sce.rds')

obs <- read.csv("../../20210611_final_object/fetal_skin.norm.maternal_removed.20220202.obs.csv.gz", row.names=1)

sce <- sce[, colnames(sce) %in% rownames(obs)]

milo <- Milo(sce)
milo <- buildGraph(milo, k=15, d=30, reduced.dim = "X_scvi_donor")
milo <- makeNhoods(milo, prop=0.05, k=15, d=30, refined=T, reduced_dims='X_scvi_donor')
milo <- calcNhoodDistance(milo, d=30, reduced.dim='X_scvi_donor')
milo <- countCells(milo, meta.data=as.data.frame(colData(milo)), sample='sanger_id')
saveRDS(milo, 'fetal_skin.scvi_donor.milo.20220202.rds')
milo <- readRDS("fetal_skin.scvi_donor.milo.20220202.rds")

milo_design <- data.frame(colData(milo)[, c('donor', 'pcw', 'sorting', 'sanger_id')])
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$sanger_id

results <- testNhoods(milo, design=~sorting + pcw, design.df=milo_design, norm.method='TMM')
saveRDS(results, 'fetal_skin.scvi_donor.milo_results.20220202.rds')
results <- readRDS("fetal_skin.scvi_donor.milo_results.20220202.rds")

## fig1b_annotation

fig1b_annotation <- obs$fig1b_annotation_20220202 %>% as.character()
fig1b_annotation[fig1b_annotation=='Progenitor'] = 'Haem progenitors'
fig1b_annotation[fig1b_annotation=='Neuronal cells'] = 'Neuronal cell'
fig1b_order <- c(
    'Haem progenitors', 'ILC', 'T cell', 'B cell', 'pDC', 'cDC', 'Langerhans cell',
    'Macrophage', 'Monocyte', 'Neutrophil', 'Mast cell', 'Megakaryocyte', 'Erythroid', 'Vascular endothelium',
    'Lymphatic endothelium', 'Mural cell', 'Skeletal muscle', 'Myofibroblast', 'Fibroblast',
    "Dermal papillia", 'Adipocyte', 'Keratinocyte', 'Melanocyte', 'Schwann cell', 'Neuronal cell'
)
colData(milo)$fig1b_annotation <- factor(fig1b_annotation, levels=fig1b_order)

colData(milo)$independent_annotation_refined <- obs$independent_annotation_refined_20220202
colData(milo)$joint_annotation <- obs$joint_annotation_20220202

ct_results <- annotateNhoods(milo, results, coldata_col='fig1b_annotation')
ct_results$nhood_annotation <- ifelse(ct_results$fig1b_annotation_fraction < 0.7, 'Mixed', ct_results$fig1b_annotation)
ordered_ct <- (ct_results %>% group_by(fig1b_annotation) %>% summarise(mean_logFC=-mean(logFC)) %>% arrange(mean_logFC))$fig1b_annotation

p1 <- (
    plotDAbeeswarm(ct_results %>% filter(nhood_annotation != 'Mixed'), group.by='nhood_annotation', alpha=0.05)
    + theme(axis.text.y=element_text(size=8))
)
p1a <- p1 + scale_x_discrete(limits=rev(fig1b_order))
q1a <- ggplot_build(p1a)
q1a$data[[1]]$size <- 0.3
gp1a <- ggplot_gtable(q1a)

p1b <- p1 + scale_x_discrete(limits=ordered_ct)
q1b <- ggplot_build(p1b)
q1b$data[[1]]$size <- 0.3
gp1b <- ggplot_gtable(q1b)

pdf(file='fetal_skin.scvi_donor.milo.fig1b_annotation_vs_pcw.20221010.pdf')
plot(gp1a)
plot(gp1b)
dev.off()

system("rclone copy --drive-shared-with-me fetal_skin.scvi_donor.milo.fig1b_annotation_vs_pcw.20221010.pdf 'google:/Fetal Skin/Figures/Figs_from_Ni/'")

## joint_annotation

jct_results <- annotateNhoods(milo, results, coldata_col='joint_annotation')
jct_results$nhood_annotation <- ifelse(jct_results$joint_annotation_fraction < 0.7, 'Mixed', jct_results$joint_annotation)

p2 <- (
    plotDAbeeswarm(jct_results %>% filter(nhood_annotation != 'Mixed'), group.by='nhood_annotation', alpha=0.05)
    + theme(axis.text.y=element_text(size=8))
)
q2 <- ggplot_build(p2)
q2$data[[1]]$size <- 0.5
gp2 <- ggplot_gtable(q2)
plot(gp2)

pdf(file='fetal_skin.scvi_donor.milo.joint_annotation_vs_pcw.20220202.pdf', height=10)
plot(gp2)
dev.off()

## independent_annotation_refined

rfd_results <- annotateNhoods(milo, results, coldata_col='independent_annotation_refined')
rfd_results$nhood_annotation <- ifelse(rfd_results$independent_annotation_refined_fraction < 0.7, 'Mixed', rfd_results$independent_annotation_refined)

p3 <- (
    plotDAbeeswarm(rfd_results %>% filter(nhood_annotation != 'Mixed'), group.by='nhood_annotation', alpha=0.05)
    + theme(axis.text.y=element_text(size=8))
)

p3_order <- c(
    "Periderm", "Early KC (stem cell)", "Basal KC", "Hair follicle", "Suprabasal",
    "Melanoblast", "Melanocyte",
    "FRZB+ early fibroblast", "HOXC5+ early fibroblast", "Pre-dermal condensate", "Dermal condensate", "Dermal papillia",
    "WNT2+ fibroblast", "Myofibroblasts", "PEAR1+ fibroblast", "Adipocytes",
    "Pericytes", "PLN+ mural cell", "LMCD1+ mural cell",
    "Myoblasts", "Early myocytes", "Myocytes",
    "Neuron progenitors", "SPP1+ proliferating neuron proneitors", "Schwann/Schwann precursors",
    "PID1+ schwann cellls", "Myelinating Schwann cells", "Neuroendocrine",
    "Early LE", "LE", "Early endothelial cell", "Tip cell (arterial)", "Arterial",
    "Capillary (venular tip)", "Capillary/postcapillary venule", "Postcapillary venule",
    "HSC", "MEMP - Early erythroid", "MEMP - Megak", "Megakaryocyte",
    "Eo/baso/mast cell progenitor", "Mast cell (earliest)", "Mast cell (medium)", "Mast cell (most mature)",
    "Early erythroid (embryonic)", "Erythroid (embryonic)", "Early erythroid", "Erythroid (fetal)",
    "Lymphoid progenitor", "ILC2", "ILC3", "LTi cell", "Innate T type3",
    "Treg", "CD4 T cell", "CD8 T cell", "Innate T type1", "NK cell",
    "Pre pro B cell", "Pro B cell", "Pre B cell", "B cell",
    "Monocyte precursor", "Monocyte", "Monocyte (activated/differentiating)",
    "TREM2+ macrophage", "LYVE1++ macrophage", "MHCII+ macrophage", "Iron-recycling macrophage",
    "LC", "DC1", "DC2", "Inflammatory DC", "ASDC", "pDC",
    "Neutrophil1", "Neutrophil2"
)

p3o <- p3 + scale_x_discrete(limits=rev(p3_order))

q3o <- ggplot_build(p3o)
q3o$data[[1]]$size <- 0.5
gp3o <- ggplot_gtable(q3o)
plot(gp3o)

pdf(file='fetal_skin.scvi_donor.milo.independent_annotation_refined_vs_pcw.20220627.pdf', height=10)
plot(gp3o)
dev.off()

#####################
# get early/late cells

early_nhoods <- tibble(ct_results) %>% filter(FDR < 0.05, logFC < 0) %>% select(Nhood) %>% unlist(use.names = F)
late_nhoods <- tibble(ct_results) %>% filter(FDR < 0.05, logFC > 0) %>% select(Nhood) %>% unlist(use.names = F)

early_cells <- colnames(milo)[Matrix::rowSums(milo@nhoods[, early_nhoods]) >= 1]
late_cells <- colnames(milo)[Matrix::rowSums(milo@nhoods[, late_nhoods]) >= 1]

write.table(early_cells, "fetal_skin.scvi_donor.milo.fig1b_annotation_vs_pcw.early_cells.list", row.names = F, col.names = F, quote = F)

write.table(late_cells, "fetal_skin.scvi_donor.milo.fig1b_annotation_vs_pcw.late_cells.list", row.names = F, col.names = F, quote = F)
