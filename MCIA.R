library(tidyverse)
library(ggrepel)
library(omicade4)
library(ggvenn)
library(pheatmap)
library(ggpubr)

## Read in all datasets -----
# use gene names for everything in order to correlate across the datasets

#RNA

RNA <- read.csv('inputs/RNAseq_counts_log2normalized_bc.csv', row.names = 1)
t_RNA <- as.data.frame(t(RNA))

#filter down to 9 samples w/proteo & phospho data
RNA9 <- t_RNA[c('TB-000018-D1WBLL', 'TB-000023-D1WBLL', 'TB-000025-D1WBLL', 
                'TB-000042-D1BBLL', 'TB-000203-D1WBLL', 
                'TB-000005-D1WPLL', 'TB-000071-D1WPLL', 'TB-000115-D1WPLL', 
                'TB-000174-D1WPLL')]
colnames(RNA9) <- c('ER1', 'ER2', 'ER3', 'ER4', 'ER5', 
                    'Ph1', 'Ph2', 'Ph3', 'Ph5')
#convert to numeric dataframe 
RNA9m <- as.matrix(RNA9)
mode(RNA9m) <- "numeric"
RNA9m[is.na(RNA9m)] <- 0
RNA9 <- as.data.frame(RNA9m)
#remove any genes with all zeros
RNA9b <- RNA9[rowSums(RNA9) != 0, ]

#proteo

proteo <- read.csv('inputs/proteo_log_norm_GeneSymbol.csv', row.names = 1)
proteo <- proteo[rownames(proteo) != "", ]
#filter down to 9 samples w/RNA data
proteo9 <- proteo[c('ER_D1', 'ER_D2', 'ER_D3', 'ER_D4', 'ER_D5', 
                    'Ph_D1', 'Ph_D2', 'Ph_D3', 'Ph_D5')]
colnames(proteo9) <- c('ER1', 'ER2', 'ER3', 'ER4', 'ER5', 
                       'Ph1', 'Ph2', 'Ph3', 'Ph5')
proteo9[is.na(proteo9)] <- 0
#remove proteins with all zeros
proteo9b <- proteo9[rowSums(proteo9) != 0, ]

#phospho

phospho <- read.csv('inputs/phospho_log_norm_GeneSymbol.csv', row.names = 1)
phospho <- phospho[rownames(phospho) != "", ]
#filter down to 9 samples w/RNA data
phospho9 <- phospho[c('ER_D1', 'ER_D2', 'ER_D3', 'ER_D4', 'ER_D5', 
                      'Ph_D1', 'Ph_D2', 'Ph_D3', 'Ph_D5')]
colnames(phospho9) <- c('ER1', 'ER2', 'ER3', 'ER4', 'ER5', 
                        'Ph1', 'Ph2', 'Ph3', 'Ph5')
phospho9[is.na(phospho9)] <- 0
#remove phosphopeptides with all zeros
phospho9b <- phospho9[rowSums(phospho9) != 0, ]


## Run MCIA -----

#list assays
assays <- list(RNA9b, proteo9b, phospho9b)

#run MCIA
mcoin <- mcia(assays, cia.nf = 8)

cancer <- c(rep('ETV6::RUNX1', 5), rep('Ph-like', 4))

#plot first two PCs
pdf('MCIA.pdf', height = 6, width = 8)
plot(mcoin, axes = 1:2, 
     phenovec = cancer,
     sample.lab = TRUE,
     sample.legend = FALSE,
     sample.color = c(rep('#F8766D', 5), rep('#00BFC4', 4)),
     df.color = 'black')
dev.off()

pdf('MCIA_nolabel.pdf', height = 6, width = 8)
plot(mcoin, axes = 1:2, 
     phenovec = cancer,
     sample.lab = FALSE,
     sample.legend = FALSE,
     sample.color = c(rep('#F8766D', 5), rep('#00BFC4', 4)),
     df.color = 'black')
dev.off()

pdf('MCIA_nolabel_colors.pdf', height = 6, width = 8)
plot(mcoin, axes = 1:2, 
     phenovec = cancer,
     sample.lab = FALSE,
     sample.legend = FALSE,
     sample.color = c(rep('#F8766D', 5), rep('#00BFC4', 4)),
     df.color = c('#D75828', '#28D758', '#5828D7'))
dev.off()

# pull out contribution of each assay to PC variance

contribs <- mcoin$mcoa$lambda
contribs12 <- contribs[1:2]
colnames(contribs12) <- c('Dim1', 'Dim2')
contribs12$dataset <- c('RNA', 'proteo', 'phospho')

contribs_long <- pivot_longer(contribs12, cols = Dim1:Dim2, 
                              names_to = 'Dimension', 
                              values_to = 'Contribution')
contribs_long$dataset <- factor(contribs_long$dataset, levels = c('RNA', 'proteo', 'phospho'))

pdf('MCIA_dataset_contribution.pdf', width = 5, height = 5)
ggplot(contribs_long, aes(x = Dimension, y = Contribution, fill = dataset)) + 
  geom_bar(stat = 'identity', position = position_dodge()) + theme_classic() + 
  scale_fill_manual(values = c('#D75828', '#28D758', '#5828D7'))
dev.off()


## Look at features driving subtype separation -----

#pull out features driving separation across PC2
Ph_features <- selectVar(mcoin, a2.lim = c(1, Inf))
ER_features <- selectVar(mcoin, a2.lim = c(-Inf, -1))

#look at overlap across assays
#manually removing HLA genes because they get swept up into one entry in MCIA processing

combined_sets <- list(RNA = unique(c(Ph_features[Ph_features$df1 == TRUE, 'var'], 
                                     ER_features[ER_features$df1 == TRUE & ER_features$var != 'HLA', 'var'])),
                      proteo = unique(c(Ph_features[Ph_features$df2 == TRUE, 'var'], 
                                        ER_features[ER_features$df2 == TRUE & ER_features$var != 'HLA', 'var'])),
                      phospho = unique(c(Ph_features[Ph_features$df3 == TRUE, 'var'], 
                                         ER_features[ER_features$df3 == TRUE & ER_features$var != 'HLA', 'var'])))

pdf('MCIA_assayVenn.pdf', height = 4, width = 4)
ggvenn(combined_sets, show_percentage = FALSE, 
       fill_color = c('#D75828', '#28D758', '#5828D7'))
dev.off()

# pull out overlapping features for each subtype
Ph_multi <- Ph_features[rowSums(1*Ph_features[2:4]) > 1, ]
ER_multi <- ER_features[rowSums(1*ER_features[2:4]) > 1 & ER_features$var != 'HLA', ]

Ph_multi2 <- 1*(Ph_multi[2:4])
rownames(Ph_multi2) <- Ph_multi$var
colnames(Ph_multi2) <- c('RNA', 'proteo', 'phospho')
Ph_multi2$proteo <- 2*Ph_multi2$proteo
Ph_multi2$phospho <- 3*Ph_multi2$phospho

pdf('MCIA_Ph_sharedFeatures.pdf', height = 4, width = 4)
pheatmap(Ph_multi2, cluster_rows = FALSE,
         cluster_cols = FALSE, color = c('white', '#D75828', '#28D758', '#5828D7'),
         angle_col = "0", legend = FALSE)
dev.off()

ER_multi2 <- 1*(ER_multi[2:4])
rownames(ER_multi2) <- ER_multi$var
colnames(ER_multi2) <- c('RNA', 'proteo', 'phospho')
ER_multi2$proteo <- 2*ER_multi2$proteo
ER_multi2$phospho <- 3*ER_multi2$phospho

pdf('MCIA_ER_sharedFeatures.pdf', height = 4, width = 4)
pheatmap(ER_multi2, cluster_rows = FALSE,
         cluster_cols = FALSE, color = c('white', '#D75828', '#28D758', '#5828D7'),
         angle_col = "0", legend = FALSE)
dev.off()


## Plotting feature expression across assays -----

# add subtype information to RNA assay
RNA_sampleTable <- read.csv('RNAseq_sampleTable.csv')

#features of interest (using GeneSymbol)
fig_goi <- c('IGF2BP1', 'MS4A1', 'SH3BP1', 'IQGAP1', 'ADD2', 'BCLAF1', 
             'HSPB1', 'TMPO', 'CHD3', 'ICAM3')
#phosphopeptides of interest
fig_poi <- c('IGF2BP1', 'MS4A1.1', 'SH3BP1', 'IQGAP1', 
             'ADD2.7', 'BCLAF1.7', 'BCLAF1.11', 
             'BCLAF1.13', 'HSPB1', 'CHD3')

RNA_goi <- RNA[ , colnames(RNA) %in% fig_goi]
RNA_goi$subtype <- RNA_sampleTable$subtype

t_proteo <- as.data.frame(t(proteo[c(1:5, 11:15)]))
#set NAs to global minimum
t_proteo[is.na(t_proteo)] <- min(t_proteo[!is.na(t_proteo)])
proteo_goi <- t_proteo[ , colnames(t_proteo) %in% fig_goi]
proteo_goi$subtype <- c(rep('ETV6::RUNX1', 5), rep('Ph-like', 5))

t_phospho <- as.data.frame(t(phospho[c(1:5, 11:15)]))
#set NAs to global minimum
t_phospho[is.na(t_phospho)] <- min(t_phospho[!is.na(t_phospho)])
phospho_goi <- t_phospho[ , colnames(t_phospho) %in% fig_poi]
phospho_goi$subtype <- c(rep('ETV6::RUNX1', 5), rep('Ph-like', 5))

colnames(phospho_goi) <- c('MS4A1-pS35,S36', 'BCLAF1-pS177', 'BCLAF1-pS658', 'BCLAF1-pS512', 
                           'ADD2-pS455', 'CHD3-pS1660,S1664', 'HSPB1-pS15', 'IGF2BP1-pS181', 
                           'IQGAP1-pS1443', 'SH3BP1-pS544', 'subtype')


#plots colored according to subtype
for(gene in colnames(RNA_goi[1:10])){
  filename <- paste0('features_subtype/RNA_', gene, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(RNA_goi, aes(x = subtype, y = .data[[gene]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(gene) + ylab('log(counts)'))
  dev.off()
}

for(protein in colnames(proteo_goi[1:10])){
  filename <- paste0('features_subtype/proteo_', protein, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(proteo_goi, aes(x = subtype, y = .data[[protein]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)'))
  dev.off()
}

for(phosphopeptide in colnames(phospho_goi[1:10])){
  print(phosphopeptide)
  filename <- paste0('features_subtype/phospho_', phosphopeptide, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(phospho_goi, aes(x = subtype, y = .data[[phosphopeptide]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(phosphopeptide) + ylab('log(counts)'))
  dev.off()
}

#plots colored according to assay
for(gene in colnames(RNA_goi[1:10])){
  filename <- paste0('features_assay/RNA_', gene, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(RNA_goi, aes(x = subtype, y = .data[[gene]])) +
          geom_boxplot(fill = '#D75828') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(gene) + ylab('log(counts)'))
  dev.off()
}

for(protein in colnames(proteo_goi[1:10])){
  filename <- paste0('features_assay/proteo_', protein, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(proteo_goi, aes(x = subtype, y = .data[[protein]])) +
          geom_boxplot(fill = '#28D758') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)'))
  dev.off()
}

for(phosphopeptide in colnames(phospho_goi[1:10])){
  filename <- paste0('features_assay/phospho_', phosphopeptide, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(phospho_goi, aes(x = subtype, y = .data[[phosphopeptide]])) +
          geom_boxplot(fill = '#5828D7') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(phosphopeptide) + ylab('log(counts)'))
  dev.off()
}

## using ggpubr to add p-values
#here computing wilcox comparisons
# * = p < 0.05
# ** = p < 0.01
# *** = p < 0.001
# **** = p < 0.0001

#plots colored according to subtype w/p-values
for(gene in colnames(RNA_goi[1:10])){
  filename <- paste0('features_subtype_pval/RNA_', gene, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(RNA_goi, aes(x = subtype, y = .data[[gene]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(gene) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}

for(protein in colnames(proteo_goi[1:10])){
  filename <- paste0('features_subtype_pval/proteo_', protein, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(proteo_goi, aes(x = subtype, y = .data[[protein]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}

for(phosphopeptide in colnames(phospho_goi[1:10])){
  filename <- paste0('features_subtype_pval/phospho_', phosphopeptide, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(phospho_goi, aes(x = subtype, y = .data[[phosphopeptide]], fill = subtype)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(phosphopeptide) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}

#plots colored according to assay w/p-values
for(gene in colnames(RNA_goi[1:10])){
  filename <- paste0('features_assay_pval/RNA_', gene, '.pdf')
  pdf(filename, height = 4, width = 2)
  print(ggplot(RNA_goi, aes(x = subtype, y = .data[[gene]])) +
          geom_boxplot(fill = '#D75828') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(gene) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}

for(protein in colnames(proteo_goi[1:10])){
  filename <- paste0('features_assay_pval/proteo_', protein, '.pdf')
  pdf(filename, height = 4, width = 2)
  print(ggplot(proteo_goi, aes(x = subtype, y = .data[[protein]])) +
          geom_boxplot(fill = '#28D758') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}

for(phosphopeptide in colnames(phospho_goi[1:10])){
  filename <- paste0('features_assay_pval/phospho_', phosphopeptide, '.pdf')
  pdf(filename, height = 4, width = 2)
  print(ggplot(phospho_goi, aes(x = subtype, y = .data[[phosphopeptide]])) +
          geom_boxplot(fill = '#5828D7') + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(phosphopeptide) + ylab('log(counts)') + 
          stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c('ETV6::RUNX1', 'Ph-like')), 
                             label = "p.signif"))
  dev.off()
}


## Combining all datasets for Cytoscape visualization -----

# RNA

RNA_in <- read.csv('RNAseq_DE_allGenes.csv')
RNA_stats <- RNA_in[c('gene', 'log2FoldChange', 'padj')]
colnames(RNA_stats) <- c('GeneSymbol', 'RNA_log2FC', 'RNA_padj')

#cap log2FC at +/- 6.64 (to match phospho & proteo data)
RNA_stats$RNA_log2FC[RNA_stats$RNA_log2FC < -6.64] <- -6.64
RNA_stats$RNA_log2FC[RNA_stats$RNA_log2FC > 6.64] <- 6.64

# proteo

proteo_in <- read.csv('24-04-10_Hs_PLL_BLL_LFQ_UnEnr.csv')
proteo_in <- proteo_in[1:2988, ]

proteo_stats <- proteo_in[c('Accession', 'Gene.Symbol',
                         'Abundance.Ratio..log2....PLL..Diseased.....BLL..Diseased.', 
                         'Abundance.Ratio.P.Value...PLL..Diseased.....BLL..Diseased.')]
colnames(proteo_stats) <- c('Accession', 'GeneSymbol', 'proteo_log2FC', 'proteo_pval')

# KSEA output

#kinase targets
ksea <- read.csv('KSEA_DXout_kinase_targets.csv')
ksea$Gene_Mod <- paste(ksea$Substrate.Gene, ksea$Substrate.Mod, sep = '-')

#kinase activity
ksea_kinases <- read.csv('KSEA_DXout_kinase_table.csv')

ksea <- merge(ksea, ksea_kinases)


##make new table for everything

all_table <- ksea[c('Kinase.Gene', 'Substrate.Gene', 'Substrate.Mod', 'Gene_Mod', 
                    'log2FC', 'Enrichment', 'z.score')]
colnames(all_table) <- c('Kinase', 'GeneSymbol', 'Modification', 'Gene_Mod', 
                         'kinase_log2FC', 'kinase_enrichment', 'kinase_zscore')


#check how many kinases are not listed as targets
kinase_list <- unique(ksea$Kinase.Gene)
length(kinase_list[!(kinase_list %in% ksea$Substrate.Gene)])
#130

#new data frame to add kinase genes (not in Substrate.Gene) to all_table
add_df <- data.frame(rep('', 130), kinase_list[!(kinase_list %in% ksea$Substrate.Gene)], 
                     rep('', 130), kinase_list[!(kinase_list %in% ksea$Substrate.Gene)], 
                     rep('', 130), rep('', 130), rep('', 130))
colnames(add_df) <- c('Kinase', 'GeneSymbol', 'Modification', 'Gene_Mod', 
                      'kinase_log2FC', 'kinase_enrichment', 'kinase_zscore')

all_table <- rbind(all_table, add_df)


#add proteo data
all_table2 <- merge(all_table, proteo_stats, all.x = TRUE)
#add RNA data
all_table3 <- merge(all_table2, RNA_stats, all.x = TRUE)

#check which kinase genes don't have abundance data
kinase_missing <- all_table3[all_table3$kinase_log2FC == '' & is.na(all_table3$proteo_log2FC) & 
                               is.na(all_table3$RNA_log2FC), 'GeneSymbol']

#remove any kinases with no abundance data
all_table4 <- all_table3[!(all_table3$GeneSymbol %in% kinase_missing), ]

#save this table
write.csv(all_table4, 'kinase_multiome_table.csv', row.names = FALSE, quote = FALSE)
