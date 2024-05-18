library(ggrepel)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(clusterProfiler)
library(ggpubr)
library(corrgram)
library(RColorBrewer)
library(matrixStats)
library(ggpmisc)
library(dplyr)
library(tidyverse)
library(cowplot)
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
# A - volcano plot
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


# Some additional exploration
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

# B - correlation plots (by @lavrentiidanilov)
head(countData)

all_genemap <- select(org.Sc.sgd.db, keys=countData$Geneid, 
                    columns=c('ORF', 'GENENAME'), keytype="ORF")
countData$gene_name <- ifelse(is.na(all_genemap$GENENAME), all_genemap$ORF,
                                          all_genemap$GENENAME)
head(countData)
countData_new <- as_data_frame(countData)

wt35_1.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(wt35_1.3st_sum = sum(wt35_1.3st, na.rm = TRUE))
wt35_2.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(wt35_2.3st_sum = sum(wt35_2.3st, na.rm = TRUE))
wt35_3.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(wt35_3.3st_sum = sum(wt35_3.3st, na.rm = TRUE))
wt35_4.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(wt35_4.3st_sum = sum(wt35_4.3st, na.rm = TRUE))

Mut218_1.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(Mut218_1.3st_sum = sum(Mut218_1.3st, na.rm = TRUE))
Mut218_2.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(Mut218_2.3st_sum = sum(Mut218_2.3st, na.rm = TRUE))
Mut218_3.3st_sum <- countData_new %>% group_by(gene_name) %>% summarise(Mut218_3.3st_sum = sum(Mut218_3.3st, na.rm = TRUE))

#put all data frames into list
df_list <- list(wt35_1.3st_sum, wt35_2.3st_sum, wt35_3.3st_sum, wt35_4.3st_sum, 
                Mut218_1.3st_sum, Mut218_2.3st_sum, Mut218_3.3st_sum)      
#merge all data frames together
total_trascriptome <- bind_cols(df_list)[, (1:7)*2]
total_trascriptome <- apply(total_trascriptome, 2, function(x) x/countData_new$Length)
total_trascriptome <- apply(total_trascriptome, 2, function(x) x/sum(x) * 1000000)
total_trascriptome <- bind_cols(wt35_1.3st_sum$gene_name, total_trascriptome)
colnames(total_trascriptome)[1] = 'Protein'
head(total_trascriptome)
total_trascriptome <- replace(total_trascriptome, total_trascriptome==0, NA)

total_trans_finall <- na.omit(total_trascriptome)
total_trans_finall$Median_WT_T = rowMedians(as.matrix(total_trans_finall[,c(2:5)]))
total_trans_finall$Median_MUT_T = rowMedians(as.matrix(total_trans_finall[,c(6:8)]))

data_prot <- read.csv('Sup35_report_unique_genes_matrix.tsv', sep = '\t')
data_prot[is.na(data_prot)] <-  0
data_prot$Median_WT_P = rowMedians(as.matrix(data_prot[,c(6:9)]))
data_prot$Median_MUT_P = rowMedians(as.matrix(data_prot[,c(2:5)]))

TranscriptomeData <-  total_trans_finall[,c(1, 9:10)]
ProteomeData = data_prot[, c(1,10:11)]
colnames(ProteomeData) <- c('Protein', 'Median_WT_P', 'Median_MUT_P')
Full_table <- merge(TranscriptomeData, ProteomeData, by="Protein")

Full_table = Full_table[apply(Full_table[, 2:5], 1, function(x) sum(x == 0)) == 0, ]
Full_table$log_wt_t = log(Full_table$Median_WT_T + 1, 2)
Full_table$log_mut_t = log(Full_table$Median_MUT_T + 1, 2)
Full_table$log_wt_p = log(Full_table$Median_WT_P + 1, 2)
Full_table$log_mut_p = log(Full_table$Median_MUT_P + 1, 2)

cor.test(Full_table$log_wt_t, Full_table$log_wt_p, method = 'kendall')
cor.test(Full_table$log_mut_t, Full_table$log_mut_p, method = 'kendall')

cor_wt_p <- ggplot(Full_table, aes(x = log_wt_t, y = log_wt_p)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se=F) +
  theme_bw() +
  xlab('log(normalized read count)') +
  ylab('log(protein intensity)') +
  ggtitle('Correlation for strains with wild-type Sup35') +
  annotate("text", x = 12, y = 4, label = 'tau = 0.256')

cor_mut_p <- ggplot(Full_table, aes(x = log_mut_t, y = log_mut_p)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se=F) +
  theme_bw() +
  xlab('log(normalized read count)') +
  ylab('log(protein intensity)') +
  ggtitle('Correlation for strains with  sup35-218 mutation') +
  annotate("text", x = 12, y = 4, label = 'tau = 0.238')

plot_grid(cor_wt_p, cor_mut_p, labels = c("B", "C"),
          ncol = 2, nrow = 1)
