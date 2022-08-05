#evaluate published DIPG signatures
#none can compete with  our meta signatures
setwd("~/Documents/msu/DIPG/DIPG/DIPG_release/data")

library("readxl")
library("octad")
library("octad.db")
dz_meta = read.csv("meta_dz_signature_GPS.csv")[, c("GeneSymbol", "Value")]
dim(dz_meta)

dz_Anastas_KDM1A = read_excel("raw/Anastas_CancerCell_2019_mmc5.xlsx", sheet = 1, skip = 3)[, c(1, 3)]
colnames(dz_Anastas_KDM1A) = c("GeneSymbol", "Value")
dim(dz_Anastas_KDM1A)

dz_Anastas_survival= read_excel("raw/Anastas_CancerCell_2019_mmc5.xlsx", sheet = 2, skip = 3)[, c(1, 3)]
colnames(dz_Anastas_survival) = c("GeneSymbol", "Value")
dim(dz_Anastas_survival)

dz_Pathania_K27MpHGG= read_excel("raw/Pathania_CancerCell_2017ccell_2544_mmc5 2.xlsx", sheet = 8, skip = 2)[, c(8, 2)]
colnames(dz_Pathania_K27MpHGG) = c("GeneSymbol", "Value")
dim(dz_Pathania_K27MpHGG)

dz_Paugh_DIPGvsHGG= read_excel("raw/Paugh_JCO_2011_Table_S4_online_only.xls", sheet = 1, skip = 0)[, c(8, 2)]
colnames(dz_Paugh_DIPGvsHGG) = c("GeneSymbol", "Value")
dim(dz_Paugh_DIPGvsHGG)

dz_Saratsis_DIPGvsOther = read_excel("raw/Saratsis_Acta_2013_401_2013_1218_MOESM7_ESM.xlsx", sheet = 1, skip = 0)[, c(1, 8)]
colnames(dz_Saratsis_DIPGvsOther) = c("GeneSymbol", "Value")
dim(dz_Saratsis_DIPGvsOther)

dz_Paugh_DIPGvsNormal= read_excel("raw/Saratsis_Acta_2013_401_2013_1218_MOESM13_ESM.xlsx", sheet = 1, skip = 0)[, c(1, 5)]
colnames(dz_Paugh_DIPGvsNormal) = c("GeneSymbol", "Value")
dim(dz_Paugh_DIPGvsNormal)

############
sens.df <- readxl::read_xlsx("raw/DIPG_cellActivity_summarized from Matrix by Billy.xlsx", sheet = "Combined Table", skip = 1)
sens.df = sens.df[sens.df$CCLASS2 < 4 & sens.df$CCLASS < 4, ]
sens.df$pert_iname = sens.df$Name
pert_efficacy1 = aggregate(AC50 ~ pert_iname, sens.df, median)
colnames(pert_efficacy1) = c("pert_iname", "efficacy")
##############
#from pubchem
pubchem = read.csv("raw/PubChem_DIPG_efficacy_auc.csv")
lincs_cmpd_info = read.csv("raw/lincs_cmpd_info.csv")
pubchem_auc = merge(pubchem, lincs_cmpd_info, by.x = "PUBCHEM_CID", by.y = "pubchem_cid"  )
pert_efficacy2 = aggregate(AUC ~ pert_iname, pubchem_auc, median)
colnames(pert_efficacy2) = c("pert_iname", "efficacy")
##########

pert_efficacy = pert_efficacy2
dz = dz_meta 
srges=runsRGES(data.frame(Symbol = dz$GeneSymbol, log2FoldChange = dz$Value), choose_fda_drugs = T, max_gene_size=100,permutations=100) #conpute sRGES
srges = srges[srges$sRGES < -0.1, ]
srges_ac50 = merge(srges, pert_efficacy, by = "pert_iname")
plot(srges_ac50$sRGES, srges_ac50$efficacy)
cor.test(srges_ac50$sRGES, srges_ac50$efficacy, method = "spearman")
#AC50, srges<-0.1   all srges, AUC srges<-0.1, AUC
#dz_meta 0.61 0.05 0.27
#dz_Anastas_KDM1A  -0.12 0.12 0.16
#dz_Anastas_survival 0.03 0.17 0.01
#dz_Pathania_K27MpHGG 0.1 0.24 0.18
#dz_Paugh_DIPGvsHGG 0.21 0.03 0.11
#dz_Saratsis_DIPGvsOther 0.03 0.06 0.1
#dz_Paugh_DIPGvsNormal -0.16 -0.05 0.04



