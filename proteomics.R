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

#read in all data
proteo <- read.csv('inputs/proteo_input.csv')

#filter out contaminants
#since there is a gap between signal and contaminants, look for it
proteo[proteo$Master == '', 1:4]
#rows 2989-91 are blank
proteo <- proteo[1:2988, ]


### working with sample abundances -----

#pull out abundances (not normalized)
proteo_abundance <- proteo[c(3, 4, 6, 7, 48:67)]

#ER = ETV6::RUNX1, Ph = Ph-like
#D = diagnosis, N = normal (remission)
colnames(proteo_abundance) <- c('Accession', 'GeneSymbol', 'NumPeptides', 'NumPSMs',
                                'ER_D1', 'ER_D2', 'ER_D3', 'ER_D4', 'ER_D5', 
                                'ER_N1', 'ER_N2', 'ER_N3', 'ER_N4', 'ER_N5', 
                                'Ph_D1', 'Ph_D2', 'Ph_D3', 'Ph_D4', 'Ph_D5', 
                                'Ph_N1', 'Ph_N2', 'Ph_N3', 'Ph_N4', 'Ph_N5')

proteo_matrix <- proteo_abundance[proteo_abundance$Accession != '', ]
rownames(proteo_matrix) <- proteo_matrix$Accession
proteo_matrix <- proteo_matrix[5:24]

proteo_log <- log2(proteo_matrix)

## QC & normalization using proDA package -----

#number of missing values 
barplot(colSums(is.na(proteo_log)), ylab = "# missing values", xlab = "sample #")

#intensity distribution
boxplot(proteo_log, ylab = "Intensity distribution", xlab = "sample #")

#applying 'conservative' median normalization
norm_proteo_log <- median_normalization(as.matrix(proteo_log))
boxplot(norm_proteo_log, ylab = "Intensity distribution", xlab = "sample #")

#save w/Accession as row names
write.csv(norm_proteo_log, "proteo_log_norm_Accession.csv", row.names = TRUE)

## Correlation matrix -----

# use gene names (made unique) as row names for downstream processing

norm_proteo_log2 <- norm_proteo_log
rownames(norm_proteo_log2) <- make.unique(proteo_abundance$GeneSymbol)

#save w/GeneSymbol as row names
write.csv(norm_proteo_log2, "proteo_log_norm_GeneSymbol.csv", row.names = TRUE)

#correlation using base R cor function
#replace NAs with 0s
norm_proteo_log3 <- norm_proteo_log2
norm_proteo_log3[is.na(norm_proteo_log3)] <- 0

proteo_cor <- cor(norm_proteo_log3)
pheatmap(proteo_cor, scale = 'none')

ann_data <- as.data.frame(colnames(norm_proteo_log))
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

pdf('proteo_corr_matrix.pdf', width = 5, height = 4)
pheatmap(proteo_cor, scale = 'none', 
         annotation_col = ann_data[c('subtype', 'timepoint')], 
         annotation_colors = ann_colors,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

pdf('proteo_corr_matrix_label.pdf', width = 5, height = 4)
pheatmap(proteo_cor, scale = 'none', 
         annotation_col = ann_data[c('subtype', 'timepoint')], 
         annotation_colors = ann_colors)
dev.off()


## PCA -----
#transpose for PCA
#purposely use matrix with 0s since prcomp doesn't handle NAs

t_proteo <- as.data.frame(t(norm_proteo_log3))
colnames(t_proteo) <- rownames(norm_proteo_log3)

proteo_pca <- prcomp(t_proteo, center = TRUE)

#scree plot
pdf('proteo_PCA_Scree.pdf', width = 4, height = 4)
fviz_eig(proteo_pca)
dev.off()

#PCA plot, dimensions 1 and 3
pdf('proteo_PCA_PC13.pdf', width = 8, height = 6)
fviz_pca_ind(proteo_pca, habillage = ann_data$subtype_timepoint, geom = "point",
             mean.point = FALSE, pointsize = 4, axes = c(1, 3)) + 
  scale_shape_manual(values = rep(16, 4)) +
  geom_line(aes(group = ann_data$patient), color = 'grey')
dev.off()

#pull out top drivers of PCs
res_var <- get_pca_var(proteo_pca)

res_var2 <- as.data.frame(res_var$contrib)

vars_PC1 <- res_var2[order(res_var2$Dim.1, decreasing = TRUE), ]
vars_PC3 <- res_var2[order(res_var2$Dim.3, decreasing = TRUE), ]

vars_plot <- unique(c(rownames(vars_PC1)[1:5], rownames(vars_PC3)[1:5]))

#top 5 contributing proteins to PC1 and PC3
pdf('proteo_PCA_5vars_PC13.pdf', width = 5, height = 5)
fviz_pca_var(proteo_pca,
             col.var = "contrib", 
             repel = TRUE,
             select.var = list(name = vars_plot), 
             axes = c(1, 3)) + 
  scale_color_gradient2(low = 'grey', high = 'black')
dev.off()



### Differential expression

## pull out log2FC and p-value for all comparisons ----

proteo_stats <- proteo[c('Accession', 'Gene.Symbol', 
                         'Abundance.Ratio..log2....BLL..Diseased.....BLL..Normal.',
                         'Abundance.Ratio..log2....PLL..Diseased.....BLL..Diseased.',
                         'Abundance.Ratio..log2....PLL..Normal.....BLL..Normal.',
                         'Abundance.Ratio..log2....PLL..Diseased.....PLL..Normal.',
                         'Abundance.Ratio.P.Value...BLL..Diseased.....BLL..Normal.',
                         'Abundance.Ratio.P.Value...PLL..Diseased.....BLL..Diseased.',
                         'Abundance.Ratio.P.Value...PLL..Normal.....BLL..Normal.',
                         'Abundance.Ratio.P.Value...PLL..Diseased.....PLL..Normal.')]
colnames(proteo_stats) <- c('Accession', 'GeneSymbol', 'log2AR_ER_DvN', 
                            'log2AR_D_PhvER', 'log2AR_N_PhvER', 'log2AR_Ph_DvN',
                            'pval_ER_DvN', 'pval_D_PhvER', 'pval_N_PhvER', 
                            'pval_Ph_DvN')

## format DE results for UpSet plots

#filters: |log2AR| > 1 & p-value < 0.05
upset_pval_df <- proteo_stats |> 
  mutate(ER_DvN_D = 1*(log2AR_ER_DvN > 1 & pval_ER_DvN < 0.05)) |>
  mutate(ER_DvN_N = 1*(log2AR_ER_DvN < -1 & pval_ER_DvN < 0.05)) |>
  mutate(D_PhvER_Ph = 1*(log2AR_D_PhvER > 1 & pval_D_PhvER < 0.05)) |>
  mutate(D_PhvER_ER = 1*(log2AR_D_PhvER < -1 & pval_D_PhvER < 0.05)) |>
  mutate(N_PhvER_Ph = 1*(log2AR_N_PhvER > 1 & pval_N_PhvER < 0.05)) |>
  mutate(N_PhvER_ER = 1*(log2AR_N_PhvER < -1 & pval_D_PhvER < 0.05)) |>
  mutate(Ph_DvN_D = 1*(log2AR_Ph_DvN > 1 & pval_Ph_DvN < 0.05)) |>
  mutate(Ph_DvN_N = 1*(log2AR_Ph_DvN < -1 & pval_Ph_DvN < 0.05))

upset_pval_df <- upset_pval_df[c(1, 2, 11:18)]

#replace all NAs with 0s
upset_pval_df[is.na(upset_pval_df)] <- 0

#save it
write.csv(upset_pval_df, 'DE_proteo_pval05_UpSet_df.csv', 
          row.names = FALSE, quote = FALSE)

#diagnosis vs. remission by subtype UpSet plot
pdf('UpSet_pval05_DvN.pdf', width = 4, height = 4)
upset(upset_pval_df, sets = c('ER_DvN_D', 'ER_DvN_N', 'Ph_DvN_D', 'Ph_DvN_N'), keep.order = T, order.by = 'freq')
dev.off()


## Pathway analysis (using gProfiler) for all comparisons -----
#using gProfiler multiquery to get p-values for all terms for all comparisons

g_query <- list()

for(i in 1:8){
  name <- colnames(upset_pval_df)[i + 2]
  print(name)
  goi <- unique(upset_pval_df[upset_pval_df[name] == 1, 'GeneSymbol'])
  print(length(goi))
  g_query[[i]] <- goi
}

names(g_query) <- colnames(upset_pval_df)[3:10]

gostres <- gost(query = g_query, organism = "hsapiens", ordered_query = FALSE,
                sources = c("GO:BP", "GO:MF", "KEGG", "REAC", "CORUM"),
                multi_query = TRUE, significant = TRUE, user_threshold = 0.05, 
                correction_method = "g_SCS", domain_scope = "annotated")

res_df <- as.data.frame(gostres$result)

res_unnest <- unnest_wider(res_df, c(p_values, significant, query_sizes, intersection_sizes), names_sep = '_')

#save this (removing 'parents' column (39) to flatten to list)
write.csv(res_unnest[1:38], 'gProfiler_proteo_multiquery.csv', row.names = FALSE)

#turn 'significant' column into UpSet df
upset_pathways <- as.data.frame(res_unnest[c(1, 36, 10:17)])
upset_pathways[3:10] <- 1*upset_pathways[3:10] 
colnames(upset_pathways)[3:10] <- colnames(upset_pval_df)[3:10]

#save it
write.csv(upset_pathways, 'proteo_pathways_UpSet_df.csv', 
          row.names = FALSE, quote = FALSE)

pdf('UpSet_pathways_DvN.pdf', width = 4, height = 4)
upset(upset_pathways, sets = c('ER_DvN_D', 'ER_DvN_N', 'Ph_DvN_D', 'Ph_DvN_N'), keep.order = T, order.by = 'freq')
dev.off()


### Individual patient diagnosis vs. remission differential expression -----

#pull out normalized abundance
norm_abundance <- proteo[c(3, 4, 6, 7, 28:47)]

colnames(norm_abundance) <- c('Accession', 'GeneSymbol', 'NumPeptides', 'NumPSMs',
                              'ER1_D', 'ER2_D', 'ER3_D', 'ER4_D', 'ER5_D', 
                              'ER1_N', 'ER2_N', 'ER3_N', 'ER4_N', 'ER5_N', 
                              'Ph1_D', 'Ph2_D', 'Ph3_D', 'Ph4_D', 'Ph5_D', 
                              'Ph1_N', 'Ph2_N', 'Ph3_N', 'Ph4_N', 'Ph5_N')

#add 0s to replace NA
norm_abundance[is.na(norm_abundance)] <- 0

#make pseudocount abundance
pseudocount_abundance <- norm_abundance
pseudocount_abundance[5:24] <- norm_abundance[5:24] + 1

#patients
pats <- c('ER1', 'ER2', 'ER3', 'ER4', 'ER5',
          'Ph1', 'Ph2', 'Ph3', 'Ph4', 'Ph5')

pat_DA <- norm_abundance[1:2]

for(pat in pats){
  pat_D <- paste0(pat, '_D')
  pat_N <- paste0(pat, '_N')
  pat_DA[pat] <- pseudocount_abundance[pat_D] / pseudocount_abundance[pat_N]
}

#cap fold change at 100 (following Proteome Discoverer convention)
pat_DA[3:12] <- pat_DA[3:12] %>% mutate(across(everything(), ~ ifelse(. > 100, 100, .)))
#correspondingly, set lower bound at 0.01
pat_DA[3:12] <- pat_DA[3:12] %>% mutate(across(everything(), ~ ifelse(. < 0.01, 0.01, .)))

#log2 transform
log_pat_DA <- pat_DA
log_pat_DA[3:12] <- log2(log_pat_DA[3:12])


#make UpSet plot of individual patient MD/DX fold changes
upset_DA <- log_pat_DA %>%
  mutate(ER1_D = 1*(ER1 >= 1)) %>%
  mutate(ER2_D = 1*(ER2 >= 1)) %>%
  mutate(ER3_D = 1*(ER3 >= 1)) %>%
  mutate(ER4_D = 1*(ER4 >= 1)) %>%
  mutate(ER5_D = 1*(ER5 >= 1)) %>%
  mutate(Ph1_D = 1*(Ph1 >= 1)) %>%
  mutate(Ph2_D = 1*(Ph2 >= 1)) %>%
  mutate(Ph3_D = 1*(Ph3 >= 1)) %>%
  mutate(Ph4_D = 1*(Ph4 >= 1)) %>%
  mutate(Ph5_D = 1*(Ph5 >= 1)) %>%
  mutate(ER1_N = 1*(ER1 <= -1)) %>%
  mutate(ER2_N = 1*(ER2 <= -1)) %>%
  mutate(ER3_N = 1*(ER3 <= -1)) %>%
  mutate(ER4_N = 1*(ER4 <= -1)) %>%
  mutate(ER5_N = 1*(ER5 <= -1)) %>%
  mutate(Ph1_N = 1*(Ph1 <= -1)) %>%
  mutate(Ph2_N = 1*(Ph2 <= -1)) %>%
  mutate(Ph3_N = 1*(Ph3 <= -1)) %>%
  mutate(Ph4_N = 1*(Ph4 <= -1)) %>%
  mutate(Ph5_N = 1*(Ph5 <= -1))

upset_DA <- upset_DA[c(1,2, 13:32)]

write.csv(upset_DA, 'patientDA_UpSet_df.csv', row.names = FALSE)

## individual patient diagnosis vs. remission pathway analysis -----

#list of terms (empty)
pathway_terms <- data.frame()

for(geneset in colnames(upset_DA)[3:22]){
  print(geneset)
  goi <- upset_DA[upset_DA[geneset] == 1, ]
  genelist <- unique(goi$GeneSymbol)
  print(length(genelist))
  gostres <- gost(query = genelist, organism = "hsapiens", ordered_query = FALSE, 
                  sources = c("GO:BP", "GO:MF", "KEGG", "REAC", "CORUM"),
                  multi_query = FALSE, significant = TRUE, user_threshold = 0.05, 
                  correction_method = "g_SCS", domain_scope = "annotated")
  res_df <- as.data.frame(gostres$result)
  print(dim(res_df))
  #assign to a variable
  var_name <- paste0(geneset, '_pathways')
  assign(var_name, res_df)
  term_ids <- res_df[c('term_id', 'term_name')]
  pathway_terms <- unique(rbind(pathway_terms, term_ids))
}

#set up the dataframe
pathways_df <- pathway_terms

for(geneset in colnames(upset_DA)[3:22]){
  df_name <- paste0(geneset, '_pathways')
  res_df <- get(df_name)
  pathways_df[geneset] <- 1*pathways_df$term_id %in% res_df$term_id
}

#save this 
write.csv(pathways_df, 'patientDA_pathway_UpSet_df.csv', row.names = FALSE)

#visualize as heatmap

pathways_heatmap <- pheatmap(pathways_df[3:22])
pathways_clust <- data.frame(cluster = cutree(pathways_heatmap$tree_row, k = 6))
pathways_clust$cluster <- as.character(pathways_clust$cluster)

clust_cols <- list(cluster = c('1' = '#F564E3', '2' = '#F8766D', '3' =  '#B79F00',
                               '4' = '#00BA38', '5' = '#619CFF', '6' = '#00BFC4'))

pdf('patientDA_pathway_heatmap.pdf', height = 8, width = 6)
pheatmap(pathways_df[3:22], scale = "none", cluster_cols = FALSE, 
         show_rownames = FALSE, color = c('white', 'black'),
         cutree_rows = 6, annotation_row = pathways_clust, 
         annotation_colors = clust_cols)
dev.off()

clustered_pathways <- cbind(pathways_df, pathways_clust)

#cluster legend
# 1: pretty much everything; ignore
# 2: ETV6-RUNX1 diagnosis!
# 3: patient-specific stuff; ignore
# 4: diagnosis regardless of subtype
# 5: remission regardless of subtype
# 6: Ph-like diagnosis!

#pull out subtype-specific diagnosis-enriched pathways
pathways_clust2 <- clustered_pathways[clustered_pathways$cluster == '2', ]
pathways_clust6 <- clustered_pathways[clustered_pathways$cluster == '6', ]


## Visualize top differentially expressed proteins

#filter UpSet dataframe to Ph vs. ER DE proteins at diagnosis
upset_df_DX <- upset_pval_df[upset_pval_df$D_PhvER_Ph == 1 | upset_pval_df$D_PhvER_ER == 1, ]
#553 proteins

#filter normalized log abundance dataframe down to DX only
DX_proteo <- as.data.frame(norm_proteo_log2)[c(1:5, 11:15)]

#pull out DE proteins only
sub_proteo <- DX_proteo[rownames(DX_proteo) %in% upset_df_DX$GeneSymbol, ]
#at least 6 non-zero (non-NA) values
sub_proteo2 <- sub_proteo[rowSums(is.na(sub_proteo)) < 7, ]
#216 proteins

#replace NAs with 0s
sub_proteo2[is.na(sub_proteo2)] <- 0

heatmap_cols = list(subtype = c(`ETV6::RUNX1` = '#F8766D', `Ph-like` = '#00BFC4'))

pdf('proteo_heatmap_DEproteins.pdf', height = 3, width = 8)
pheatmap(as.data.frame(t(sub_proteo2)), scale = 'none', 
         annotation_row = ann_data[c('subtype')], 
         annotation_colors = heatmap_cols, 
         show_colnames = FALSE)
dev.off()