library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(DESeq2)
library(EnhancedVolcano)
library(UpSetR)
library(gprofiler2)
library(org.Hs.eg.db)
library(Seurat) 

#read in sample information
sampleTable <- read.csv('inputs/RNAseq_sampleTable.csv')

#read in counts
counts <- read.csv('inputs/RNAseq_counts.csv', row.names = 1)
#replace '.' in colnames with '-'
colnames(counts) <- gsub('\\.', '-', colnames(counts))

## Set up DESeq2 object -----

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sampleTable,
                              design = ~ batch + subtype)
#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~subtype, sampleTable)
mat <- limma::removeBatchEffect(mat, vsd$batch, design = design0)
assay(vsd) <- mat

DESeq2::plotPCA(vsd, intgroup = "subtype")

## Differential expression -----

dds <- DESeq(dds)

res <- results(dds, contrast = c("subtype", "Ph-like", "ETV6::RUNX1"))

DESeq2::plotMA(res) 

#convert to tibble
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

#save full output
write.csv(res_tbl,'RNAseq_DE_allGenes.csv', row.names = FALSE, quote = FALSE)

#filter down to the significant results
sig_res <- dplyr::filter(res_tbl, padj < 0.05) %>%
  dplyr::arrange(padj)

dim(sig_res) #810 genes

#save sig genes
write.csv(sig_res,'RNAseq_DE_p05sigGenes.csv', row.names = FALSE, quote = FALSE)


## Gene expression visualization -----

#read in log2 normalized, batch corrected counts
counts_data <- read.csv('inputs/RNAseq_counts_log2normalized_bc.csv', row.names = 1)
                        
#visualize CRLF2 expression
pdf('RNAseq_CRLF2exp.pdf', height = 4, width = 4)
ggplot(counts_data, aes(x = sampleTable$subtype, y = CRLF2)) + 
  geom_boxplot(aes(fill = sampleTable$subtype)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) + 
  theme_classic() + ggtitle('CRLF2 expression') + ylab('log(counts)')
dev.off()

#visualize top 50 DE genes

top50_sigres <- sig_res %>% head(n = 50)

top50_counts <- counts_data[ , colnames(counts_data) %in% top50_sigres$gene]

rownames(sampleTable) <- sampleTable$ID

heatmap_cols = list(subtype = c(`ETV6::RUNX1` = '#F8766D', `Ph-like` = '#00BFC4'))

pdf('RNAseq_heatmap_top50DGE.pdf', height = 5, width = 10)
pheatmap(top50_counts[1:50], scale = 'none', annotation_row = sampleTable['subtype'],
         annotation_colors = heatmap_cols)
dev.off()


## Pathway analysis using gProfiler -----

#break up significantly DE genes by direction
Ph <- sig_res[sig_res$log2FoldChange > 0, ]
ER <- sig_res[sig_res$log2FoldChange < 0, ]

#filter by log2FC (for consistency w/proteomics)
sigPh <- Ph[Ph$log2FoldChange > 1, ] #344
sigER <- ER[ER$log2FoldChange < -1, ] #366

#run gprofiler *multiquery* to get p-values for all pathways for all lists
g_query_RNA <- list(sigPh$gene, sigER$gene)

names(g_query_RNA) <- c('Ph-like', 'ETV6::RUNX1')

gostres_RNA <- gost(query = g_query_RNA, organism = "hsapiens", ordered_query = FALSE,
                    sources = c("GO:BP", "GO:MF", "KEGG", "REAC", "CORUM"),
                    multi_query = TRUE, significant = TRUE, user_threshold = 0.05, 
                    correction_method = "g_SCS", domain_scope = "annotated")

res_df_RNA <- as.data.frame(gostres_RNA$result)

res_unnest_RNA <- unnest_wider(res_df_RNA, c(p_values, significant, query_sizes, intersection_sizes), names_sep = '_')

#save this, cutting out 'parents' column which doesn't flatten nicely
write.csv(res_unnest_RNA[1:14], 'gProfiler_RNA_multiquery.csv', row.names = FALSE)

