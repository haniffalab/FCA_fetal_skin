loadPackage(tidyverse)
loadPackage(pheatmap)
loadPackage(viridis)

make_celltype_matrix <- function(df, var) {
    df1 <- (df
        %>% group_by(celltype_pair)
        %>% summarise(n=length(.data[[var]]))
        %>% separate(col='celltype_pair', into=c('celltype_a', 'celltype_b'), sep='[|]')
        %>% spread(key='celltype_b', value='n', fill=0))
    mat <- as.matrix(df1[, -1])
    rownames(mat) <- df1$celltype_a
    mat
}

make_gene_matrix <- function(df, var) {
    df1 <- (df
        %>% group_by(gene_pair)
        %>% summarise(n=length(.data[[var]]))
        %>% separate(col='gene_pair', into=c('gene_a', 'gene_b'), sep='[|]')
        %>% spread(key='gene_b', value='n', fill=0))
    mat <- as.matrix(df1[, -1])
    rownames(mat) <- df1$gene_a
    mat
}

read_cellphonedb_output <- function(outdir) {
    pvalues_tsv <- paste(outdir, 'pvalues.tsv', sep='/')
    if (file.exists(pvalues_tsv)) {
        means_tsv <- paste(outdir, 'means.tsv', sep='/')
    } else {
        means_tsv <- paste(outdir, 'significant_means.tsv', sep='/')
    }

    expr_df <- read_tsv(means_tsv)
    expr_df1 <- (expr_df
               %>% filter(secreted, (receptor_a | receptor_b), !is_integrin, gene_a!="", gene_b!="")
               %>% gather(key='celltype_pair', value='mean', 13:ncol(expr_df))
               %>% filter(!is.na(mean))
               %>% select(interacting_pair, gene_a, gene_b, secreted, receptor_a, receptor_b, celltype_pair, mean)
               %>% separate(col='celltype_pair', into=c('celltype_a', 'celltype_b'), sep='[|]')
               %>% mutate(
                   gene_pair=ifelse(receptor_b, paste(gene_a, gene_b, sep='|'), paste(gene_b, gene_a, sep='|')),
                   celltype_pair=ifelse(receptor_b, paste(celltype_a, celltype_b, sep='|'), paste(celltype_b, celltype_a, sep='|')))
               %>% select(gene_pair, celltype_pair, mean))

    if (file.exists(pvalues_tsv)) {
        p_df <- read_tsv(pvalues_tsv)
        padj_df <- p_df
        for (cname in colnames(padj_df)[13:ncol(padj_df)]) {
            padj_df[[cname]] <- p.adjust(padj_df[[cname]], method='fdr')
        }
        p_df1 <- (p_df
                   %>% filter(secreted, (receptor_a | receptor_b), !is_integrin, gene_a!="", gene_b!="")
                   %>% gather(key='celltype_pair', value='p', 13:ncol(p_df))
                   %>% select(interacting_pair, gene_a, gene_b, secreted, receptor_a, receptor_b, celltype_pair, p)
                   %>% separate(col='celltype_pair', into=c('celltype_a', 'celltype_b'), sep='[|]')
                   %>% mutate(
                       gene_pair=ifelse(receptor_b, paste(gene_a, gene_b, sep='|'), paste(gene_b, gene_a, sep='|')),
                       celltype_pair=ifelse(receptor_b, paste(celltype_a, celltype_b, sep='|'), paste(celltype_b, celltype_a, sep='|')))
                   %>% select(gene_pair, celltype_pair, p))
        padj_df1 <- (padj_df
                   %>% filter(secreted, (receptor_a | receptor_b), !is_integrin, gene_a!="", gene_b!="")
                   %>% gather(key='celltype_pair', value='p', 13:ncol(padj_df))
                   %>% select(interacting_pair, gene_a, gene_b, secreted, receptor_a, receptor_b, celltype_pair, p)
                   %>% separate(col='celltype_pair', into=c('celltype_a', 'celltype_b'), sep='[|]')
                   %>% mutate(
                       gene_pair=ifelse(receptor_b, paste(gene_a, gene_b, sep='|'), paste(gene_b, gene_a, sep='|')),
                       celltype_pair=ifelse(receptor_b, paste(celltype_a, celltype_b, sep='|'), paste(celltype_b, celltype_a, sep='|')))
                   %>% select(gene_pair, celltype_pair, p))

        expr_df1$p <- p_df1$p
        expr_df1$padj <- padj_df1$p
    }

    if (file.exists(pvalues_tsv)) {
        expr_df1 <- expr_df1 %>% filter(p < 0.05)
    }
    expr_df1
}

get_cutoff <- function(x) {
    dst <- density(x)
    mode.x <- dst$x[which.max(dst$y)]
    mad.x <- median(abs(x-mode.x))
    mode.x + 3*mad.x
}

args <- commandArgs(trailingOnly=T)
if (length(args) > 0) {
    outdir <- args[1]
} else {
    outdir <- 'all'
}

if (!dir.exists(outdir)) {
    stop(paste(outdir, 'not found'))
}

expr_tsv <- paste0(outdir, '/cellphonedb_summary.tsv')
if (file.exists(expr_tsv)) {
    expr_df1 <- read_tsv(expr_tsv)
} else {
    expr_df1 <- read_cellphonedb_output(outdir)
    expr_df1 %>% write_tsv(expr_tsv)
}

min_mean <- get_cutoff(expr_df1$mean)

expr_df2 <- expr_df1 %>% filter(padj > 0.05)
expr_df3 <- expr_df1 %>% filter(mean > min_mean)
expr_df4 <- expr_df2 %>% filter(mean > min_mean)

expr_ct_mat1 <- expr_df1 %>% make_celltype_matrix('mean')
expr_g_mat1 <- expr_df1 %>% make_gene_matrix('mean')
expr_ct_mat2 <- expr_df2 %>% make_celltype_matrix('mean')
expr_g_mat2 <- expr_df2 %>% make_gene_matrix('mean')
expr_ct_mat3 <- expr_df3 %>% make_celltype_matrix('mean')
expr_g_mat3 <- expr_df3 %>% make_gene_matrix('mean')
expr_ct_mat4 <- expr_df4 %>% make_celltype_matrix('mean')
expr_g_mat4 <- expr_df4 %>% make_gene_matrix('mean')

pdf(file=paste0(outdir, '/cellphonedb_summary.pdf'))
hist(expr_df1$mean, 100, freq=F, main='significant co-expression distribution')
lines(density(expr_df1$mean))
abline(v=min_mean)
mtext(text=paste('cutoff =', round(min_mean, 3)), side=3, adj=1)
pheatmap(expr_ct_mat1, fontsize=5, main='p < 0.05')
pheatmap(log1p(expr_g_mat1), fontsize=5, main='p < 0.05')
pheatmap(expr_ct_mat2, fontsize=5, main='padj < 0.05')
pheatmap(log1p(expr_g_mat2), fontsize=5, main='padj < 0.05')
pheatmap(expr_ct_mat3, fontsize=5, main='p < 0.05, mean > cutoff')
pheatmap(log1p(expr_g_mat3), fontsize=5, main='p < 0.05, mean > cutoff')
pheatmap(expr_ct_mat4, fontsize=5, main='padj < 0.05, mean > cutoff')
pheatmap(log1p(expr_g_mat4), fontsize=5, main='padj < 0.05, mean > cutoff')
dev.off()
