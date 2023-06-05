loadPackage(tidyverse)
cpdb_out <- read_tsv('../20201107_organoid_cellphonedb/fetal_skin_all_vs_all/cellphonedb_summary.tsv')
dim(cpdb_out)
ct_pairs <- read_tsv('cpdb_ligand_receptor_celltypes.tsv')
gene_pairs <- read_tsv('cpdb_ligand_receptor_genes.tsv')
colnames(ct_pairs) <- c('ligand', 'receptor')
colnames(gene_pairs) <- c('pair')
ct_pairs$pair <- paste(ct_pairs$ligand, ct_pairs$receptor, sep='|')

df1 <- cpdb_out %>% filter(gene_pair %in% gene_pairs$pair, celltype_pair %in% ct_pairs$pair)
print(
      df1 %>%
          mutate(gene_pair=str_replace(gene_pair, '\\|', ' | '), celltype_pair=str_replace(celltype_pair, '\\|', ' | ')) %>%
          ggplot(aes(x=celltype_pair, y=gene_pair, color=mean, size=padj)) +
          geom_point() + scale_size_continuous(trans='reverse') +
          scale_color_distiller(palette='OrRd') +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), panel.grid=element_blank())
)
dev.copy2pdf(file='cpdb_stroma_for_hudaa_20211208a.pdf', width=12, height=6)
savehistory()
