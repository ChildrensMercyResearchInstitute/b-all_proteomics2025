library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(proxy)
library(proDA)
library(EnhancedVolcano)
library(UpSetR)
library(gprofiler2)
library(org.Hs.eg.db)
library(Seurat) 
library(qdapRegex)

#read in all data
phospho <- read.csv('inputs/phospho_input.csv')

#filter out contaminants
#since there is a gap between signal and contaminants, look for it
phospho[phospho$Master.Protein.Accessions == '', 1:4]
#rows 1801-3 are blank
phospho <- phospho[1:1800, ]


### working with sample abundances -----

#pull out abundances (not normalized)
ph_abundance <- phospho[c(3, 5, 6, 9, 55:74)]

#ER = ETV6::RUNX1, Ph = Ph-like
#D = diagnosis, N = normal (remission)
colnames(ph_abundance) <- c('Accession', 'GeneSymbol', 'Position', 'Modifications',
                            'ER_D1', 'ER_D2', 'ER_D3', 'ER_D4', 'ER_D5', 
                            'ER_N1', 'ER_N2', 'ER_N3', 'ER_N4', 'ER_N5', 
                            'Ph_D1', 'Ph_D2', 'Ph_D3', 'Ph_D4', 'Ph_D5', 
                            'Ph_N1', 'Ph_N2', 'Ph_N3', 'Ph_N4', 'Ph_N5')

ph_matrix <- ph_abundance[ph_abundance$Accession != '', ]
ph_matrix$uAccession <- make.unique(ph_matrix$Accession)
rownames(ph_matrix) <- ph_matrix$uAccession
ph_matrix <- ph_matrix[5:24]

ph_log <- log2(ph_matrix)

## QC & normalization using proDA package -----

#number of missing values 
barplot(colSums(is.na(ph_log)), ylab = "# missing values", xlab = "sample #")

#intensity distribution
boxplot(ph_log, ylab = "Intensity distribution", xlab = "sample #")

#applying 'conservative' median normalization
norm_ph_log <- median_normalization(as.matrix(ph_log))
boxplot(norm_ph_log, ylab = "Intensity distribution", xlab = "sample #")

#save w/(unique) Accession as row names
write.csv(norm_ph_log, "phospho_log_norm_Accession.csv", row.names = TRUE)

## Correlation matrix -----

# use gene names (made unique) as row names for downstream processing

norm_ph_log2 <- norm_ph_log
rownames(norm_ph_log2) <- make.unique(ph_abundance$GeneSymbol)

#save w/GeneSymbol as row names
write.csv(norm_ph_log2, "phospho_log_norm_GeneSymbol.csv", row.names = TRUE)

#correlation using base R cor function
#replace NAs with 0s
norm_ph_log3 <- norm_ph_log2
norm_ph_log3[is.na(norm_ph_log3)] <- 0

ph_cor <- cor(norm_ph_log3)
pheatmap(ph_cor, scale = 'none')

ann_data <- as.data.frame(colnames(norm_ph_log))
colnames(ann_data) <- 'sample'
rownames(ann_data) <- ann_data$sample
ann_data$patient <- paste0(substr(ann_data$sample, 1, 2), substr(ann_data$sample, 5, 5))
ann_data$subtype <- c(rep('ETV6::RUNX1', 10), rep('Ph-like', 10))
ann_data$timepoint <- c(rep('Diagnosis', 5), rep('Remission', 5), 
                        rep('Diagnosis', 5), rep('Remission', 5))
ann_data$subtype_timepoint <- paste0(ann_data$subtype, '-', ann_data$timepoint)

ann_colors = list(
  subtype = c(`ETV6::RUNX1` = '#F8766D', `Ph-like` = '#00BFC4'), 
  timepoint = c(Diagnosis = 'white', Remission = 'black')
)

pdf('phospho_corr_matrix.pdf', width = 5, height = 4)
pheatmap(ph_cor, scale = 'none', 
         annotation_col = ann_data[c('subtype', 'timepoint')], 
         annotation_colors = ann_colors,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

pdf('phospho_corr_matrix_label.pdf', width = 5, height = 4)
pheatmap(ph_cor, scale = 'none', 
         annotation_col = ann_data[c('subtype', 'timepoint')], 
         annotation_colors = ann_colors)
dev.off()


## PCA -----
#transpose for PCA
#purposely use matrix with 0s since prcomp doesn't handle NAs

t_ph <- as.data.frame(t(norm_ph_log3))
colnames(t_ph) <- rownames(norm_ph_log3)

ph_pca <- prcomp(t_ph, center = TRUE)

#scree plot
pdf('phospho_PCA_Scree.pdf', width = 4, height = 4)
fviz_eig(ph_pca)
dev.off()

#PCA plot, dimensions 1 and 3
pdf('phospho_PCA_PC13.pdf', width = 8, height = 6)
fviz_pca_ind(ph_pca, habillage = ann_data$subtype_timepoint, geom = "point",
             mean.point = FALSE, pointsize = 4, axes = c(1, 3)) + 
  scale_shape_manual(values = rep(16, 4)) +
  geom_line(aes(group = ann_data$patient), color = 'grey')
dev.off()

#pull out top drivers of PCs
res_var <- get_pca_var(ph_pca)

res_var2 <- as.data.frame(res_var$contrib)

vars_PC1 <- res_var2[order(res_var2$Dim.1, decreasing = TRUE), ]
vars_PC3 <- res_var2[order(res_var2$Dim.3, decreasing = TRUE), ]

vars_plot <- unique(c(rownames(vars_PC1)[1:5], rownames(vars_PC3)[1:5]))

#top 5 contributing proteins to PC1 and PC3
pdf('phospho_PCA_5vars_PC13.pdf', width = 5, height = 5)
fviz_pca_var(ph_pca,
             col.var = "contrib", 
             repel = TRUE,
             select.var = list(name = vars_plot), 
             axes = c(1, 3)) + 
  scale_color_gradient2(low = 'grey', high = 'black')
dev.off()


### Set up for kinase activity Ph vs. ER at diagnosis using KSEA App -----

KSEA_stats <- phospho[c('Master.Protein.Accessions', 'Gene.Symbol', 
                        'Modifications.in.Master.Proteins..all.Sites.',
                        'Annotated.Sequence',
                        'Abundance.Ratio...PLL..Diseased.....BLL..Diseased.',
                        'Abundance.Ratio..log2....PLL..Diseased.....BLL..Diseased.',
                        'Abundance.Ratio.P.Value...PLL..Diseased.....BLL..Diseased.')]
colnames(KSEA_stats) <- c('Accession', 'GeneSymbol', 'Modifications', 'Sequence', 
                          'AR', 'log2AR', 'pval')

#unlist any modifications that have ; between them
KSEA_stats$Mods2 <- strsplit(as.character(KSEA_stats$Modifications), '];')

stats2 <- unnest(KSEA_stats, Mods2)

#fix any [] that got broken by just adding a ] to all Mods2
stats2$Mods2 <- paste0(stats2$Mods2, ']')
stats2$Mods2 <- trimws(stats2$Mods2)

#pull out accession for any proteins that have ';' in Accession column
stats2$Accession2 <- vapply(strsplit(as.character(stats2$Mods2), ' '), '[', '', 1)

#pull out just content in [] in Mods2 column
stats2$ModPos <- rm_between(stats2$Mods2, '[', ']', extract = TRUE)

#remove whitespace
stats2$ModPos <- gsub(' ', '', stats2$ModPos)

#now fix peptide sequence formatting
#select for content between '.'s
stats2$Sequence2 <- rm_between(stats2$Sequence, '.', '.', extract = TRUE)

#note that KSEA input should be FC NOT log2FC
KSEA_DX <- stats2[c('Accession2', 'GeneSymbol', 'Sequence2', 'ModPos',
                    'pval', 'AR')]
colnames(KSEA_DX) <- c('Protein', 'Gene', 'Peptide', 'Residue.Both', 'p', 'FC')

#remove lines with no gene name, or NA FC
KSEA_DX <- KSEA_DX[KSEA_DX$Gene != '', ]
KSEA_DX <- KSEA_DX[!is.na(KSEA_DX$FC), ]

write.csv(KSEA_DX, 'KSEA_in_DX.csv', row.names = FALSE, quote = FALSE)

shiny::runGitHub('KSEA', 'casecpb')

## Visualize results from KSEA App -----

DX_phospho_all <- read.csv('KSEA_DXout_kinase_table.csv')

#filter down to columns of interest
DX_phospho <- DX_phospho_all[c('Kinase.Gene', 'z.score', 'p.value')]

sig_kinases <- DX_phospho[DX_phospho$p.value < 0.1, ]
sig_kinases <- sig_kinases[order(sig_kinases$z.score, decreasing = FALSE), ]
sig_kinases$Kinase.Gene <- factor(sig_kinases$Kinase.Gene, levels = sig_kinases$Kinase.Gene)

pdf('phospho_DX_kinases.pdf', width = 8, height = 3)
ggplot(sig_kinases, aes(x = Kinase.Gene, y = z.score, fill = -log10(p.value))) + 
  scale_fill_continuous(type = 'viridis', limits = c(0, 5.5)) + 
  geom_bar(stat = 'identity') + theme_bw() + labs(fill = '-log10(p-value)') & RotatedAxis()
dev.off()



