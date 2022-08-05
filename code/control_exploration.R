#explore different controls
dz = subset(phenoDF,cancer=='diffuse intrinsic pontine glioma') #or diffuse intrinsic pontine glioma ; non-brainstem high-grade glioma; brain lower grade glioma  Non-Brainstem High-Grade Glioma Diffuse Intrinsic Pontine Glioma Brain Lower Grade Glioma
case_id1 = dz$sample.id

load("raw/DIPG_CBTTC.RData")
case_id2 = colnames(log2.read.count.matrix)

all_tissues = unique(phenoDF$biopsy.site)
brain_tissue = all_tissues[grep("BRAIN ", all_tissues, ignore.case = T)][1:13] #14 and 15 do not have normal tissues

cors = NULL
for (i in 1:length(brain_tissue)){
  control_id1 = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c(brain_tissue[i])]
  res_dipg_jude=diffExp(case_id1,control_id1,source='octad.whole',DE_method='edgeR', output=T,n_topGenes=5000,file='~/Downloads/octad.counts.and.tpm.h5', normalize_samples=T)
  
  control_id2 = phenoDF$sample.id[phenoDF$sample.type %in% c("normal") & phenoDF$biopsy.site %in% c(brain_tissue[i])]
  
  expression_log2_control=loadOctadCounts(c( control_id2),type='counts',file='~/Downloads/octad.counts.and.tpm.h5')
  expression_log2_case = log2.read.count.matrix
  
  #for some reason, it cannot run
  res_CBTTC=diffExp(case_id2,control_id2,source='',expSet = cbind(expression_log2_control, expression_log2_case), normalize_samples=T)
  
  
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
  
  cor_test =  cor.test(srges_ac50$sRGES, srges_ac50$AC50)
  print(cor_test)
  cors = rbind(cors, data.frame(tissue = brain_tissue[i], cor = cor_test$estimate, p = cor_test$p.value))
}

write.csv(cors, "reference_tissue_scores.csv")

