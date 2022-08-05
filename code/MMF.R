#MMF reversal
library("dplyr")
library("octad")
library("ggplot2")
library("ggrepel")

load("raw/DIPG_12samples.RData")
write.csv(DIPG_log2.read.count.matrix, "DIPG_log2.read.count.matrix.csv")
write.csv(DIPG_log2.tpm.matrix, "DIPG_log2.tpm.matrix.csv")

dz_signature = read.csv("meta_dz_signature.csv")
dz_signature <- dz_signature %>% filter(abs(log2FoldChange) > 1)

DIPG_log2.tpm.matrix_averge <- apply(DIPG_log2.tpm.matrix, 1, mean)
DIPG_log2.tpm.matrix <- DIPG_log2.tpm.matrix[DIPG_log2.tpm.matrix_averge > 4 , ]

DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, DMSO = log(apply(2^DIPG_log2.tpm.matrix[, 1:3], 1, median)),2)
DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, MMF = log(apply(2^DIPG_log2.tpm.matrix[, 4:6], 1, median)),2)
DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, TPL = log(apply(2^DIPG_log2.tpm.matrix[, 7:9], 1, median)),2)
DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, TRIAM = log(apply(2^DIPG_log2.tpm.matrix[, 10:12], 1, median)),2)

DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, MMF_DMSO = DIPG_log2.tpm.matrix[, "MMF"] - DIPG_log2.tpm.matrix[, "DMSO"])
DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, TPL_DMSO = DIPG_log2.tpm.matrix[, "TPL"] - DIPG_log2.tpm.matrix[, "DMSO"])
DIPG_log2.tpm.matrix <- cbind(DIPG_log2.tpm.matrix, TRIAM_DMSO = DIPG_log2.tpm.matrix[, "TRIAM"] - DIPG_log2.tpm.matrix[, "DMSO"])

dz_drug <- merge(DIPG_log2.tpm.matrix,dz_signature, by.x = 0, by.y = "identifier.x" )
dz_drug <- dz_drug[abs(dz_drug$log2FoldChange) > 2, ]
dz_drug <- dz_drug[order(abs(dz_drug$log2FoldChange)), ]
cor.test(dz_drug$MMF_DMSO, dz_drug$log2FoldChange, method = "spearman")

##can we label important genes
plot(dz_drug$MMF_DMSO, dz_drug$log2FoldChange)

cor.test(dz_drug$TPL_DMSO, dz_drug$log2FoldChange, method = "spearman")
cor.test(dz_drug$TRIAM_DMSO, dz_drug$log2FoldChange, method = "spearman")
cor.test(dz_drug$TRIAM_DMSO, dz_drug$MMF_DMSO, method = "spearman")

drug_dz_signature <- subset(dz_drug, select=c("log2FoldChange", "MMF_DMSO", "TRIAM_DMSO")) #"TPL_DMSO",
drug_dz_signature_rank <- drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature_rank[, i] <- rank(drug_dz_signature[,i])
}

drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank$log2FoldChange), ]

pdf("MMF_reversal.pdf")
colPal <- redblue(100)
par(mar=c(6, 4, 2, 0.5))
#image(t(druggable_targets_pos), col=redblue(2))
image(t(drug_dz_signature_rank), col= colPal,   axes=F, srt=45)
axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature_rank ) ), labels= F)
text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
     labels = c( "DIPG", "MMF",  "Triametrene"), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4)
dev.off()

#########
#DE genes
#
dz_signature=diffExp(c("MMF_1","MMF_2", "MMF_3"), c("DMSO_1", "DMSO_2", "DMSO_3"),source='',expSet = DIPG_log2.read.count.matrix, normalize_samples=F)
dz_signature$diffexpressed='NO'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange>1,'diffexpressed']='UP'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange<(-1),'diffexpressed']='DOWN'
table(dz_signature$diffexpressed)
#label DE genes
dz_signature$delabel <- NA
dz_signature$delabel[dz_signature$diffexpressed != "NO" & dz_signature$log2FoldChange ] <- as.character(dz_signature$Symbol[dz_signature$diffexpressed != "NO"])
#fancy volcano plot
p=ggplot(data=dz_signature, aes(x=log2FoldChange , y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p

write.csv(dz_signature, "MMF_sig.csv")

#
dz_signature=diffExp(c("TPL_1","TPL_1", "TPL_1"), c("DMSO_1", "DMSO_2", "DMSO_3"),source='',expSet = DIPG_log2.read.count.matrix, normalize_samples=F)
dz_signature$diffexpressed='NO'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange>1,'diffexpressed']='UP'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange<(-1),'diffexpressed']='DOWN'
table(dz_signature$diffexpressed)
#label DE genes
dz_signature$delabel <- NA
dz_signature$delabel[dz_signature$diffexpressed != "NO" & dz_signature$log2FoldChange ] <- as.character(dz_signature$Symbol[dz_signature$diffexpressed != "NO"])
#fancy volcano plot
p=ggplot(data=dz_signature, aes(x=log2FoldChange , y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p

write.csv(dz_signature, "TPL_sig.csv")

#
dz_signature=diffExp(c("TRIAM_1","TRIAM_2", "TRIAM_3"), c("DMSO_1", "DMSO_2", "DMSO_3"),source='',expSet = DIPG_log2.read.count.matrix, normalize_samples=F)
dz_signature$diffexpressed='NO'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange>1,'diffexpressed']='UP'
dz_signature[dz_signature$padj<0.05 & dz_signature$log2FoldChange<(-1),'diffexpressed']='DOWN'
table(dz_signature$diffexpressed)
#label DE genes
dz_signature$delabel <- NA
dz_signature$delabel[dz_signature$diffexpressed != "NO" & dz_signature$log2FoldChange ] <- as.character(dz_signature$Symbol[dz_signature$diffexpressed != "NO"])
#fancy volcano plot
p=ggplot(data=dz_signature, aes(x=log2FoldChange , y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p

write.csv(dz_signature, "TRIAM_sig.csv")

#identify reservsed genes
dz_signature = read.csv("meta_dz_signature.csv")
mmf_signature = read.csv("MMF_sig.csv")

dz_mmf = merge(dz_signature[, c("log2FoldChange", "padj", "Symbol")], mmf_signature[, c("log2FoldChange", "padj", "Symbol")], by = "Symbol")
dz_mmf = dz_mmf[abs(dz_mmf$log2FoldChange.x) > 1 & abs(dz_mmf$log2FoldChange.y) > 1, ]
dz_mmf$genelabel = ""

dz_mmf$genelabel[(dz_mmf$log2FoldChange.x > 1 & dz_mmf$log2FoldChange.y < -1) |
                (dz_mmf$log2FoldChange.x < -1 & dz_mmf$log2FoldChange.y > 1)  ] = as.character(dz_mmf$Symbol[(dz_mmf$log2FoldChange.x > 1 & dz_mmf$log2FoldChange.y < -1) |
                                                                                                               (dz_mmf$log2FoldChange.x < -1 & dz_mmf$log2FoldChange.y > 1)] )

p=ggplot(data=dz_mmf, aes(x=log2FoldChange.x , y= log2FoldChange.y,   label=genelabel)) +
  geom_point() +
  geom_text_repel() 
  
p
  

library("enrichplot")
library("clusterProfiler")
library("GO.db")
library(org.Hs.eg.db)
mmf_signature = read.csv("MMF_sig.csv")
gene_ids = mmf_signature$GeneID[mmf_signature$log2FoldChange  < -1]
gene_ids = gene_ids[!is.na(gene_ids)]

ego <- enrichGO(gene          = gene_ids,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
pdf("MMF_down_BP.pdf")
dotplot(ego, showCategory=10) + ggtitle("down-regulated GO Biological Processes")
dev.off()
