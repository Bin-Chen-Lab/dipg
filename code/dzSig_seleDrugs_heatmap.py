#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:20:04 2021

@author: jingxing
"""

import pandas as pd
import seaborn as sns
import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 14

wrkdir = './'
lincs_sig = pd.read_hdf('/Users/jingxing/Documents/OneDrive - Michigan State University/Data_Sources/LINCS_profiles/L1000_P1P2_Lv5_landmark_highQ_cmpd_median.hdf5')
sele_drugs = ['triptolide', 'mycophenolic-acid', 'mycophenolate-mofetil', 'triamterene', 'relcovaptan', 'zatebradine', 'trioxsalen']
#sele_drug_names = ['IMD-0354 (1)', 'Puromycin (22)', 'Methotrexate (110)', 'Methylene Blue (148)', 'Dasatinib (391)', 'Ethambutol (1299)', 'Betulinic Acid (1355)']
sele_drug_names = sele_drugs
metasig = pd.read_csv(wrkdir + 'meta_dz_signature-101821.csv', index_col='Symbol')
metasig = metasig[metasig.index.isin(lincs_sig.index)]

lincs_sig = lincs_sig[sele_drugs]
metasig = metasig.sort_values(by='log2FoldChange')
lincs_sig = lincs_sig.reindex(metasig.index)
lincs_sig.columns = sele_drug_names

FIG, axes = plt.subplots(3, 1, figsize=(8, 5), gridspec_kw={'height_ratios': [1, 4, 3]})
sns.heatmap(metasig[['log2FoldChange']].T, center=0, vmin=-5, vmax=5, cmap='coolwarm', \
            xticklabels=False, cbar=False, ax=axes[0])
sns.heatmap(lincs_sig.iloc[:,:4].T, center=0, vmin=-1.5, vmax=1.5, cmap='coolwarm', cbar=False, \
            xticklabels=False, ax=axes[1])
sns.heatmap(lincs_sig.iloc[:,-3:].T, center=0, vmin=-1.5, vmax=1.5, cmap='coolwarm', cbar=False, \
            xticklabels=False, ax=axes[2])
axes[0].set_title('DIPG Signature', fontsize=14)
axes[0].set_xlabel(None)
axes[0].set_yticklabels(['Log2 Fold Change'], rotation=0)
axes[1].set_title('Selected Candidates', fontsize=14)
axes[1].set_xlabel(None)
axes[1].set_yticklabels(axes[1].get_yticklabels())
axes[2].set_title('Control Compounds', fontsize=14)
axes[2].set_xlabel(None)
plt.tight_layout()
FIG.savefig(wrkdir + 'DIPG_candidates_heatmap.pdf', transparent=True)