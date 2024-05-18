library(ggrepel)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(org.Sc.sgd.db)
setwd("/media/barbitoff/DATA/Working issues/FizGen/Translation and Genome/rnaseq/paper/sup35n_expression_analysis/")

# Count matrix import
countData <- read.csv('featureCounts_sup35.tsv', header = TRUE, sep = "\t")
dim(countData)
сutData <- countData[, c(7:16)]
rownames(сutData) = countData$Geneid
count_matrix <- as.matrix(сutData)

# Condition file import
coldata <- read.csv('Condition_35.csv', row.names=1, header = TRUE, sep = "\t")
coldata$condition <- factor(coldata$condition)
coldata$batch <- factor(coldata$batch)
coldata$Extraction = factor(coldata$Extraction)
all(rownames(coldata) %in% colnames(count_matrix))

# PCA for all samples
dds_all <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition )
rlog_dds_all <- rlog(dds_all, blind = T,fitType = "mean")
write.table(assay(rlog_dds_all), "rlog_all.csv",quote = F,sep = "\t",row.names = T)
plotPCA(rlog_dds_all, intgroup = "condition", ntop = 500)  


# DE analysis stage 1 vs stage 3
count_matrix_st1vs3 <- count_matrix[, c(1:3, 8:10)]
coldata_st1vs3 <- coldata[c(1:3, 8:10), ]
dds_st1vs3 <- DESeqDataSetFromMatrix(countData = count_matrix_st1vs3,
                              colData = coldata_st1vs3,
                              design = ~ condition )
dds_st1vs3$condition <- relevel(dds_st1vs3$condition, ref = "Stage_1_wt")
dds_st1vs3 <- DESeq(dds_st1vs3)
res_st1vs3  <- results(dds_st1vs3, contrast=c("condition","Stage_3_mut","Stage_1_wt"))
res_df_st1vs3 <- as.data.frame(res_st1vs3)


# DE analysis stage 3 vs stage 3
count_matrix_st3vs3 <- count_matrix[, c(4:10)]
coldata_st3vs3 <- coldata[c(4:10), ]
dds_st3vs3 <- DESeqDataSetFromMatrix(countData = count_matrix_st3vs3,
                                     colData = coldata_st3vs3,
                                     design = ~ condition )
dds_st3vs3$condition <- relevel(dds_st3vs3$condition, ref = "Stage_3_wt")
dds_st3vs3 <- DESeq(dds_st3vs3)
res_st3vs3  <- results(dds_st3vs3, contrast=c("condition","Stage_3_mut","Stage_3_wt"))
res_df_st3vs3 <- as.data.frame(res_st3vs3)


# Counting DEGs
deg_count <- data.frame(regulation = rep(c('up', 'down'), 6),
                        reference = rep(c('stage1', 'stage3'), each=6),
                        logFC_cutoff = rep(rep(c(0, 0.5, 1), each=2), 2))
count_degs <- function(regulation, ref_stage, cutoff) {
  cutoff = as.numeric(cutoff)
  if (ref_stage == 'stage1') { 
    res_df = res_df_st1vs3 
  } else {
    res_df = res_df_st3vs3
  }
  res_df <- na.omit(res_df[res_df$padj < 0.05, ])
  if (regulation == 'up') {
    return(sum(res_df$log2FoldChange > cutoff))
  } else {
    return(sum(res_df$log2FoldChange < -1 * cutoff))
  }
}
deg_count$nGenes = apply(deg_count, 1, function(x) count_degs(x[1], x[2], x[3]))
ggplot(deg_count, aes(x=reference, y=nGenes, fill=regulation)) +
  geom_bar(stat='identity', position='dodge', col='black') + 
  theme_bw() + facet_wrap(~logFC_cutoff, nrow=1) + 
  theme(legend.position = 'top')


# Drawing volcano plot - logFC > 0.5 cutoff
res_df_st3vs3$regulation <- "none"
res_df_st3vs3$regulation[res_df_st3vs3$log2FoldChange > 0.5 & 
                           res_df_st3vs3$padj < 0.05] <- "up"
res_df_st3vs3$regulation[res_df_st3vs3$log2FoldChange < -0.5 & 
                           res_df_st3vs3$padj < 0.05] <- "down"

res_genemap <- select(org.Sc.sgd.db, keys=rownames(res_df_st3vs3), 
                      columns=c('ORF', 'GENENAME'), keytype="ORF")
res_df_st3vs3$final_id <- ifelse(is.na(res_genemap$GENENAME), res_genemap$ORF,
                                      res_genemap$GENENAME)
res_df_st3vs3$volcano_label = res_df_st3vs3$final_id
res_df_st3vs3$volcano_label[res_df_st3vs3$padj > 0.001 & 
                              abs(res_df_st3vs3$log2FoldChange) < 1.5] = NA

ggplot(data=res_df_st3vs3, aes(x=log2FoldChange, y=-log10(padj), col=regulation)) + 
  geom_point() + theme_bw() + geom_text_repel(aes(label = volcano_label)) +
  scale_color_manual(values=c("none" = "#999999", "up" = "#E86850", "down" = "#587498"))
write.table(res_df_st3vs3,"0_Lfc_35_3wt3wt.csv",quote = F,sep = "\t",row.names = T)


# Exporting DEGs for external analysis - no logFC cutoff!
res_df_st3vs3$gene_id = rownames(res_df_st3vs3)
down_degs <- res_df_st3vs3[res_df_st3vs3$log2FoldChange < 0 & res_df_st3vs3$padj < 0.05, ]
up_degs <- res_df_st3vs3[res_df_st3vs3$log2FoldChange > 0 & res_df_st3vs3$padj < 0.05, ]
write.table(down_degs, "0_DOWN_35_3wt3wt.csv",quote = F,sep = "\t",row.names = F)
write.table(up_degs, "0_UP_35_3wt3wt.csv",quote = F,sep = "\t",row.names = F)

  
#### GO term enrichment

# Biological process
library(clusterProfiler)
gobp_down <- enrichGO(gene = down_degs$gene_id,
                     OrgDb = "org.Sc.sgd.db",
                     keyType = "ORF",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)
goplot(gobp_down)
dotplot(gobp_down) 
write.table(gobp_down,"DOWN_35_lfc0_BP.csv",quote = F,sep = "\t",row.names = T)

gobp_up <- enrichGO(gene = up_degs$gene_id,
                     OrgDb = "org.Sc.sgd.db",
                     keyType = "ORF",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)
goplot(gobp_up)
dotplot(gobp_up) 
write.table(gobp_up,"UP_35_lfc0_BP.csv",quote = F,sep = "\t",row.names = T)

# Cellular component
library(clusterProfiler)
gocc_down <- enrichGO(gene = down_degs$gene_id,
                      OrgDb = "org.Sc.sgd.db",
                      keyType = "ORF",
                      ont = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05)
goplot(gocc_down)
dotplot(gocc_down) 
write.table(gocc_down,"DOWN_35_lfc0_CC.csv",quote = F,sep = "\t",row.names = T)

gocc_up <- enrichGO(gene = up_degs$gene_id,
                    OrgDb = "org.Sc.sgd.db",
                    keyType = "ORF",
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05)
goplot(gocc_up)
dotplot(gocc_up) 
write.table(gocc_up,"UP_35_lfc0_CC.csv",quote = F,sep = "\t",row.names = T)

#Molecular function
library(clusterProfiler)
gomf_down <- enrichGO(gene = down_degs$gene_id,
                      OrgDb = "org.Sc.sgd.db",
                      keyType = "ORF",
                      ont = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05)
goplot(gomf_down)
dotplot(gomf_down) 
write.table(gomf_down,"DOWN_35_lfc0_MF.csv",quote = F,sep = "\t",row.names = T)

gomf_up <- enrichGO(gene = up_degs$gene_id,
                    OrgDb = "org.Sc.sgd.db",
                    keyType = "ORF",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05)
goplot(gomf_up)
dotplot(gomf_up) 
write.table(gomf_up,"UP_35_lfc0_MF.csv",quote = F,sep = "\t",row.names = T)


### Figure 4
qpcr_df <- read.csv('qPCR_data.csv', header = TRUE, sep = "\t")
ggplot(qpcr_df, aes(x=Gene, y= Copy_number, fill=Sample))+ 
  scale_fill_brewer(palette = "Dark2")+ 
  geom_boxplot() + geom_point() + theme_bw() + 
  facet_wrap(vars(Primers), scales = "free")

wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('CDC20', qpcr_df$Gene), ])
wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('CDC23', qpcr_df$Gene), ])
wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('CDC28', qpcr_df$Gene), ])
wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('FKH1', qpcr_df$Gene), ])
wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('SUP35', qpcr_df$Gene), ])
wilcox.test(Copy_number ~ Gene, qpcr_df[grepl('SUP45', qpcr_df$Gene), ])


### Figure 5
proteomics_de <- read.table('DA_Sup35_proteomics.tsv', sep='\t', row.names=1, header=T)
head(proteomics_de)

proteomics_de$regulation <- "none"
proteomics_de$regulation[proteomics_de$logFC > 0.5 & 
                           proteomics_de$adj.P.Val < 0.05] <- "up"
proteomics_de$regulation[proteomics_de$logFC < -0.5 & 
                           proteomics_de$adj.P.Val < 0.05] <- "down"
proteomics_de$volcano_label <- rownames(proteomics_de)
proteomics_de$volcano_label[proteomics_de$regulation == 'none'] = NA

ggplot(data=proteomics_de, aes(x=logFC, y=-log10(adj.P.Val), col=regulation)) + 
  geom_point() + theme_bw() + geom_text_repel(aes(label = volcano_label)) +
  scale_color_manual(values=c("none" = "#999999", "up" = "#E86850", "down" = "#587498"))

merged_de <- merge(res_df_st3vs3, proteomics_de, by.x='final_id', by.y='row.names')
head(merged_de)

ggplot(merged_de, aes(x=log2FoldChange, y=logFC)) + geom_point() + 
  theme_bw()
cor.test(merged_de$log2FoldChange, merged_de$logFC, method='kendall')

for (lfc_cut in c(0, 0.5, 1, 2)){
  prot_sign <- (abs(merged_de$logFC) >= lfc_cut & merged_de$adj.P.Val < 0.05)
  rna_sign <- (abs(merged_de$log2FoldChange) >= lfc_cut & merged_de$padj < 0.05)
  both_sign <- (prot_sign & rna_sign)
  print(paste('For cutoff', lfc_cut, ' (proteome, RNA, both):', 
              sum(prot_sign, na.rm=T), sum(rna_sign, na.rm=T), sum(both_sign, na.rm=T)))
}


