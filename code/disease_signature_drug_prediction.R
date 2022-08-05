#disease signature creation
#The goal is to find a robust disease signature that leads to a prediction that mostly correlates with experimental data
#Q1: How to find the best reference tissues since most of datasets have no DIPG adjacent tissues?
#Q2: Does the case selection matter? We could select DIPG only, HGG all (DIPG and non brain stem), H3 mutant DIPG. 
#Q3: Do the parameters matter, e.g., number of control tissues, fold change, normalization, max gene size in drug prediction?
#Q4: Would more datasets help the prediction? at least two datasets (St. Jude; CBTTC) have enough tumor samples
#Q5: Can single cell RNASeq lead to a better prediction?

######
library(DESeq2)
library(pheatmap)
library(vegan)
library(ape)
library(rgl)
library(Rtsne)
library(ggrepel)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(dplyr)
library(ggsignif)
library(octad)
library(octad.db)
library(stringr)
library(matrixTests)

###########
#need to download octad.counts.and.tpm.h5 under ~/Download
#https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5

#drug sensitivity from NCAT matrix 
#it include four cell lines; CCLASS and CCLASS2 are used to filter low quality data
sens.df <- readxl::read_xlsx("raw/DIPG_cellActivity_summarized from Matrix by Billy.xlsx", sheet = "Combined Table", skip = 1)
sens.df = sens.df[sens.df$CCLASS2 < 4 & sens.df$CCLASS < 4, ]
sens.df$pert_iname = sens.df$Name

#SU-DIPG-IV 4-1
#JHH-DIPG-1 1-1
#SU-DIPG-XIII 13
#SU-DIPG-VI 6

#st Jude's data mutation status
#DIPG mutation; EGA mutation is downloaded from Wu paper supplementary where ACVR1 and H3K27 mutation are highlighted
EGA_sample_mapping = read.csv("raw/Sample_File.map.txt", sep="\t", header=F)
EGA_sample_mapping$sample = str_remove_all(EGA_sample_mapping$V3, ".bam.cip")
DIPG_mutation = read.csv("raw/EGA_mutation.csv")
ACVR1 = DIPG_mutation[DIPG_mutation$ACVR1 != "", ]
H3 = DIPG_mutation[DIPG_mutation$Histone3 != "", ]
ACVR1_sample = EGA_sample_mapping$sample[EGA_sample_mapping$sample %in% ACVR1$Sample ]
H3_sample = EGA_sample_mapping$sample[EGA_sample_mapping$sample %in% H3$Sample ]

#########
##########
dz = subset(phenoDF,cancer=='diffuse intrinsic pontine glioma') #or diffuse intrinsic pontine glioma ; non-brainstem high-grade glioma; brain lower grade glioma  Non-Brainstem High-Grade Glioma Diffuse Intrinsic Pontine Glioma Brain Lower Grade Glioma
case_id = dz$sample.id
#case_id = case_id[case_id %in% H3_sample] #select a subset of cases from a mutation

control_id=computeRefTissue(case_id,outputFolder='',output=T,adjacent=F,source = "octad",control_size = 100)

#from simulation, BRAIN - HIPPOCAMPUS gives the best performance
#control_id = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c("BRAIN - CAUDATE (BASAL GANGLIA)")]

res=diffExp(case_id,control_id,source='octad.whole',DE_method='edgeR', output=T,n_topGenes=5000,file='~/Downloads/octad.counts.and.tpm.h5', normalize_samples=T)

srges=runsRGES(res[abs(res$log2FoldChange)>1 & res$padj < 0.05, ],choose_fda_drugs = T, max_gene_size=100,permutations=100) #
srges = srges[srges$sRGES < -0.1, ] #removed low confidence RGES
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
##can we create a better visualization?
plot(srges_ac50$sRGES, srges_ac50$AC50)
#expected correlation: 0.47 for diffuse intrinsic pontine glioma; 0.42 non-brainstem high-grade glioma; 0.30 brain lower grade glioma
#H3 DIPG 0.422
cor.test(srges_ac50$sRGES, srges_ac50$AC50) 
t = cor.test(srges_ac50$sRGES, srges_ac50$AC50) 

p=ggplot(data=srges_ac50, label.size=15,aes(x=sRGES , y=AC50,label=pert_iname)) +
  geom_smooth(method=lm,color='red') +
  geom_point(size=3) +
  geom_text_repel() +
  theme_minimal()+
  theme(plot.title=element_text(size=20,face="bold",hjust = 0.5),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))+
  ggtitle("St Jude (bulk RNAseq)")+ #change this to modify title
  annotate('text',label=c(
    paste('r=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$estimate,2),sep='')),
    x=min(srges_ac50$sRGES)+0.004,y=max(srges_ac50$AC50),size=8, fontface = 'bold')+ #hardcoded
  annotate('text',label=c(
    paste('p=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$p.value,4),sep='')),
    x=min(srges_ac50$sRGES)+0.012,y=max(srges_ac50$AC50)-2.5,size=8, fontface = 'bold')

p
ggsave('AC50_vs_sRGES_DIPG.pdf',height=10,width=10)

res_dipg_jude = res
write.csv(res_dipg_jude, "st_jude_dz_signature.csv")

########test pubchem activity


#explore different controls
dz = subset(phenoDF,cancer=='diffuse intrinsic pontine glioma') #or diffuse intrinsic pontine glioma ; non-brainstem high-grade glioma; brain lower grade glioma  Non-Brainstem High-Grade Glioma Diffuse Intrinsic Pontine Glioma Brain Lower Grade Glioma
case_id = dz$sample.id

all_tissues = unique(phenoDF$biopsy.site)
brain_tissue = all_tissues[grep("BRAIN ", all_tissues, ignore.case = T)][1:13] #14 and 15 do not have normal tissues

cors = NULL
for (i in 1:length(brain_tissue)){
control_id = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c(brain_tissue[i])]

res=diffExp(case_id,control_id,source='octad.whole',DE_method='edgeR', output=T,n_topGenes=5000,file='~/Downloads/octad.counts.and.tpm.h5', normalize_samples=T)

srges=runsRGES(res[abs(res$log2FoldChange)>1 & res$padj < 0.05, ],choose_fda_drugs = T, max_gene_size=100,permutations=100) #
srges = srges[srges$sRGES < -0.1, ] #removed low confidence RGES
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
cor_test  = cor.test(srges_ac50$sRGES, srges_ac50$AC50) 
print(cor_test)
cors = rbind(cors, data.frame(tissue = brain_tissue[i], cor = cor_test$estimate, p = cor_test$p.value))
}

write.csv(cors, "reference_tissue_scores.csv")

#################
##can we differentiate tumor vs normal in the plot?
dipg_sample = subset(phenoDF,cancer=='diffuse intrinsic pontine glioma')$sample.id
NHGG_sample = subset(phenoDF,cancer=='non-brainstem high-grade glioma')$sample.id

tsne$type <- "Tumors"
tsne$type[tsne$sample.id %in% phenoDF$sample.id[phenoDF$sample.type %in% c("normal", "adjacent")]] <- "Normal"
tsne$type[tsne$sample.id %in% dipg_sample] <- "DIPG"
tsne$type[tsne$sample.id %in% H3_sample] <- "DIPG (H3K27M)"
tsne$type[tsne$sample.id %in% NHGG_sample] <- "Non-brainstem HGG"

#plot
(p2 <- ggplot(tsne, aes(X, Y, color = type)) + geom_point(alpha = 0.4)+
    labs(title = paste ('TNSE PLOT'), x= 'TSNE Dim1', y='TSNE Dim2', caption="OCTAD")+
    theme_bw())

## can we summerize the performance of selecting individual brain tissue samples as control (cortex?)


###########
#CBTTC data
#
load("raw/DIPG_CBTTC.RData")
case_id = colnames(log2.read.count.matrix)

#because CBTTC is not included in OCTAD, impossible to compute reference so use the reference from EGA instead
control_id=computeRefTissue(subset(phenoDF,cancer=='diffuse intrinsic pontine glioma')$sample.id,outputFolder='',output=T,adjacent=F,source = "octad",control_size = 100)

#control_id = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c("BRAIN - CAUDATE (BASAL GANGLIA)")]

#find GTEX control samples counts; using the same control as HGG
expression_log2_control=loadOctadCounts(c( control_id),type='counts',file='~/Downloads/octad.counts.and.tpm.h5')
expression_log2_case = log2.read.count.matrix

#for some reason, it cannot run
res=diffExp(case_id,control_id,source='',expSet = cbind(expression_log2_control, expression_log2_case), normalize_samples=T)

srges=runsRGES(res[abs(res$log2FoldChange)>1, ],max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$AC50)
#expected correlation: 
cor.test(srges_ac50$sRGES, srges_ac50$AC50)

res_CBTTC = res
write.csv(res_CBTTC, "CBTTC_dz_signature.csv")


p=ggplot(data=srges_ac50, label.size=15,aes(x=sRGES , y=AC50,label=pert_iname)) +
  geom_smooth(method=lm,color='red') +
  geom_point(size=3) +
  geom_text_repel() +
  theme_minimal()+
  theme(plot.title=element_text(size=20,face="bold",hjust = 0.5),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))+
  ggtitle("CBTTC (bulk RNAseq)")+ #change this to modify title
  annotate('text',label=c(
    paste('r=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$estimate,2),sep='')),
    x=min(srges_ac50$sRGES)+0.004,y=max(srges_ac50$AC50),size=8, fontface = 'bold')+ #hardcoded
  annotate('text',label=c(
    paste('p=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$p.value,4),sep='')),
    x=min(srges_ac50$sRGES)+0.012,y=max(srges_ac50$AC50)-2.5,size=8, fontface = 'bold')

p
ggsave('AC50_vs_sRGES_CBTTC.pdf',height=10,width=10)

#explore different controls
load("raw/DIPG_CBTTC.RData")
case_id = colnames(log2.read.count.matrix)

all_tissues = unique(phenoDF$biopsy.site)
brain_tissue = all_tissues[grep("BRAIN ", all_tissues, ignore.case = T)][1:13] #14 and 15 do not have normal tissues

cors = NULL
for (i in 1:length(brain_tissue)){
  control_id = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c(brain_tissue[i])]
  
  expression_log2_control=loadOctadCounts(c( control_id),type='counts',file='~/Downloads/octad.counts.and.tpm.h5')
  expression_log2_case = log2.read.count.matrix
  
  #for some reason, it cannot run
  res=diffExp(case_id,control_id,source='',expSet = cbind(expression_log2_control, expression_log2_case), normalize_samples=T)
  
  srges=runsRGES(res[abs(res$log2FoldChange)>1 & res$padj < 0.05, ],choose_fda_drugs = T, max_gene_size=100,permutations=100) #
  srges = srges[srges$sRGES < -0.1, ] #removed low confidence RGES
  pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
  srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
  cor_test  = cor.test(srges_ac50$sRGES, srges_ac50$AC50) 
  print(cor_test)
  cors = rbind(cors, data.frame(tissue = brain_tissue[i], cor = cor_test$estimate, p = cor_test$p.value))
}

write.csv(cors, "CBTTC_reference_tissue_scores.csv")



############
#combine Jude and CBTTC 
meta_sig = read.csv("meta_dz_signature.csv")
meta_sig = merge(res_dipg_jude[abs(res_dipg_jude$log2FoldChange) > 1 & res_dipg_jude$padj < 0.05,], 
                 res_CBTTC[abs(res_CBTTC$log2FoldChange) > 1 & res_CBTTC$padj < 0.05,], by = "Symbol")

meta_sig$log2FoldChange = (meta_sig$log2FoldChange.x + meta_sig$log2FoldChange.y)/2
meta_sig$padj = (meta_sig$padj.x + meta_sig$padj.y)/2
meta_sig$identifier = meta_sig$identifier.x

plot(meta_sig$log2FoldChange.y, meta_sig$log2FoldChange.x)

srges=runsRGES(data.frame(Symbol = meta_sig$Symbol, log2FoldChange = meta_sig$log2FoldChange), choose_fda_drugs = T, max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$AC50)
#expected correlation: 0.62
cor.test(srges_ac50$sRGES, srges_ac50$AC50)

p=ggplot(data=srges_ac50, label.size=15,aes(x=sRGES , y=AC50,label=pert_iname)) +
  geom_smooth(method=lm,color='red') +
  geom_point(size=3) +
  geom_text_repel() +
  theme_minimal()+
  theme(plot.title=element_text(size=20,face="bold",hjust = 0.5),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))+
  ggtitle("Combined")+ #change this to modify title
  annotate('text',label=c(
    paste('r=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$estimate,2),sep='')),
    x=min(srges_ac50$sRGES)+0.004,y=max(srges_ac50$AC50),size=8, fontface = 'bold')+ #hardcoded
  annotate('text',label=c(
    paste('p=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$p.value,5),sep='')),
    x=min(srges_ac50$sRGES)+0.012,y=max(srges_ac50$AC50)-2.5,size=8, fontface = 'bold')

p
ggsave('AC50_vs_sRGES_Meta.pdf',height=10,width=10)

##################
pubchem = read.csv("~/Downloads/PubChem_DIPG_efficacy_agg.csv")
pubchem_srges = merge(pubchem, srges, by.x = "cmap_name", by.y =  "pert_iname")
plot(pubchem_srges$sRGES, pubchem_srges$Potency)
cor.test(pubchem_srges$sRGES, pubchem_srges$Potency)
pubchem_pert_AC50 = merge(pubchem, pert_AC50, by.x = "cmap_name", by.y =  "pert_iname")
plot(pubchem_pert_AC50$PUBCHEM_ACTIVITY_SCORE, pubchem_pert_AC50$AC50)
plot(pubchem_pert_AC50$Potency, pubchem_pert_AC50$AC50)
plot(pubchem_pert_AC50$Potency, pubchem_pert_AC50$AC50)

###############

#how about only H3K27 mutation
HM_down = read.csv("~/Downloads/PATHWAY/DIPG_HM_downregulated.csv", stringsAsFactors = F)
HM_up= read.csv("~/Downloads/PATHWAY/DIPG_HM_upregulated.csv", stringsAsFactors = F)
HM_down = unique(unlist(strsplit(paste(HM_down$Genes, collapse = ";"), ";")))
HM_up = unique(unlist(strsplit(paste(HM_up$Genes, collapse = ";"), ";")))

meta_sig = read.csv("meta_dz_signature.csv")
meta_sig = meta_sig[meta_sig$Symbol %in% c(HM_up, HM_down),]

srges=runsRGES(data.frame(Symbol = meta_sig$Symbol, log2FoldChange = meta_sig$log2FoldChange), choose_fda_drugs = T, max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$AC50)
#expected to recover epigenetic inhibtiros
cor.test(srges_ac50$sRGES, srges_ac50$AC50)


##################
#can we visualize signature pathways?
GO_up = geneEnrich(meta_sig$Symbol[meta_sig$log2FoldChange.x  > 1],database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
GO_down = geneEnrich(meta_sig$Symbol[meta_sig$log2FoldChange.x  < -1],database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)

octadDrugEnrichment(sRGES = srges) #compute drug enrichment for obtained scores

write.csv(srges, "meta_sRGES.csv")
write.csv(meta_sig, "meta_dz_signature.csv")

srges_fda = read.csv("sRGES_FDAapproveddrugs.csv")
lincs_cmpd = unique(read.csv("raw/lincs_cmpd_info.csv")[, c("pert_iname", "canonical_smiles")])
lincs_fda = merge(srges_fda, lincs_cmpd, by = "pert_iname")
write.csv(lincs_fda, "sRGES_FDAapproveddrugs_smiles.csv")

#visualized enriched pathways
library("enrichplot")
library("clusterProfiler")
library("GO.db")
library(org.Hs.eg.db)

gene_ids = meta_sig$GeneID.x[meta_sig$log2FoldChange.x  > 2]
gene_ids = gene_ids[!is.na(gene_ids)]

ego <- enrichGO(gene          = gene_ids,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
pdf("up_BP.pdf")
dotplot(ego, showCategory=15) + ggtitle("up-regulated GO Biological Processes")
dev.off()

gene_ids = meta_sig$GeneID.x[meta_sig$log2FoldChange.x  < 0]
gene_ids = gene_ids[!is.na(gene_ids)]

ego <- enrichGO(gene          = gene_ids,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
pdf("down_BP.pdf")
dotplot(ego, showCategory=15) + ggtitle("down-regulated GO Biological Processes")
dev.off()

#can we visualize the reversal of the following hits
hits = c("triptolide", "mycophenolic-acid", "mycophenolate-mofetil", "triamterene", "relcovaptan", "zatebradine", "trioxsalen", "hexachlorophene")

#########
#single cell
#dowload https://singlecell.broadinstitute.org/single_cell/study/SCP147/single-cell-analysis-in-pediatric-midline-gliomas-with-histone-h3k27m-mutation
#paper https://www-science-org.proxy2.cl.msu.edu/doi/full/10.1126/science.aao4750
sc = read.csv("raw/K27Mproject.RSEM.vh20170621.txt", sep = "\t", check.names = F)
sc_meta = read.csv("raw/PortalK27M_Metadata.vh20180223.txt", sep = "\t", stringsAsFactors = F)
cancer_cells = sc_meta$NAME[sc_meta$Type == "Malignant" & sc_meta$Cellcycle > 2  ] #BCH836 MUV5 #& sc_meta$Sample == "BCH836"
normal_cells = sc_meta$NAME[sc_meta$Type == "Oligodendrocyte"  ] # & sc_meta$Sample == "BCH836"
t = row_t_equalvar(log(sc[, cancer_cells] + 0.1), log(sc[, normal_cells] + 0.1))
t = cbind(t, Symbol =  sc$Gene)
t = t[t$var.pooled > 0.1, ]
t$padj = p.adjust(t$pvalue)
t$log2FoldChange = t$mean.diff

res = t[abs(t$log2FoldChange) > 1 & t$padj < 0.05, ]

#compare with published signatures
H3_sc = read.csv("raw/aao4750_filbin_sm_tables6.csv")
H3_sc$H3K27M.high_new = sapply(H3_sc$H3K27M.high, function(x){
  str_replace_all(x, " ", "")
})

srges=runsRGES(res,max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$AC50)

p=ggplot(data=srges_ac50, label.size=15,aes(x=sRGES , y=AC50,label=pert_iname)) +
  geom_smooth(method=lm,color='red') +
  geom_point(size=3) +
  geom_text_repel() +
  theme_minimal()+
  theme(plot.title=element_text(size=20,face="bold",hjust = 0.5),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))+
  ggtitle("Broad (scRNASeq)")+ #change this to modify title
  annotate('text',label=c(
    paste('r=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$estimate,2),sep='')),
    x=min(srges_ac50$sRGES)+0.004,y=max(srges_ac50$AC50),size=8, fontface = 'bold')+ #hardcoded
  annotate('text',label=c(
    paste('p=',round(cor.test(srges_ac50$sRGES, srges_ac50$AC50)$p.value,4),sep='')),
    x=min(srges_ac50$sRGES)+0.012,y=max(srges_ac50$AC50)-2.5,size=8, fontface = 'bold')

p
ggsave('AC50_vs_sRGES_sc.pdf',height=10,width=10)



#expected correlation: 0.44
cor.test(srges_ac50$sRGES, srges_ac50$AC50)

res_sc = res
write.csv(res_sc, "sc_dz_signature.csv")

#it seems combining sc and bulk does not increase the correlation
sc_bulk_sig = merge(meta_sig[, c("Symbol","log2FoldChange",  "padj")], res_sc, by = "Symbol")

sc_bulk_sig$log2FoldChange = (sc_bulk_sig$log2FoldChange.x + sc_bulk_sig$log2FoldChange.y)/2
sc_bulk_sig$padj = (sc_bulk_sig$padj.x + sc_bulk_sig$padj.y)/2

plot(sc_bulk_sig$log2FoldChange.y, sc_bulk_sig$log2FoldChange.x)

srges=runsRGES(data.frame(Symbol = meta_sig$Symbol, log2FoldChange = meta_sig$log2FoldChange), choose_fda_drugs = T, max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
pert_AC50 = aggregate(AC50 ~ pert_iname, sens.df, median)
srges_ac50 = merge(srges, pert_AC50, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$AC50)
cor.test(srges_ac50$sRGES, srges_ac50$AC50)

up_common = intersect(meta_sig$Symbol[meta_sig$log2FoldChange > 1], res_sc$Symbol[res_sc$log2FoldChange >1])
GO_up_common = geneEnrich(up_common,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_up_common$GO_Biological_Process_2017)
#
up_bulk = meta_sig$Symbol[meta_sig$log2FoldChange > 1][!meta_sig$Symbol[meta_sig$log2FoldChange > 1] %in% up_common]
GO_up_bulk = geneEnrich(up_bulk,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_up_bulk$GO_Biological_Process_2017)
head(GO_up_bulk$GO_Cellular_Component_2017)

up_sc = res_sc$Symbol[res_sc$log2FoldChange > 1][!res_sc$Symbol[res_sc$log2FoldChange > 1] %in% up_common]
GO_up_sc = geneEnrich(up_sc,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_up_sc$GO_Biological_Process_2017)
head(GO_up_sc$GO_Cellular_Component_2017)


down_common = intersect(meta_sig$Symbol[meta_sig$log2FoldChange < -1], res_sc$Symbol[res_sc$log2FoldChange >1])
GO_down_common = geneEnrich(down_common,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_down_common$GO_Biological_Process_2017)
#
down_bulk = meta_sig$Symbol[meta_sig$log2FoldChange < -1][!meta_sig$Symbol[meta_sig$log2FoldChange < -1] %in% down_common]
GO_down_bulk = geneEnrich(down_bulk,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_down_bulk$GO_Biological_Process_2017)
head(GO_down_bulk$GO_Cellular_Component_2017)

down_sc = res_sc$Symbol[res_sc$log2FoldChange < -1][!res_sc$Symbol[res_sc$log2FoldChange < -1] %in% down_common]
GO_down_sc = geneEnrich(down_sc,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO_down_sc$GO_Biological_Process_2017)
head(GO_down_sc$GO_Cellular_Component_2017)
