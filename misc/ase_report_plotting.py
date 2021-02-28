#!/usr/bin/env python3
# -*- coding : utf-8 -*_

import pandas as pd
import matplotlib.pyplot as plt

wkdir = '/home/umcg-zzhang/Documents/projects/wp_ase_dlp/outputs/aseReport'
fpath = wkdir + '/ase_report-pval_matrix_sorted.csv'
mtrx = pd.read_csv(fpath, low_memory=False, index_col=0)
cols = [x for x in mtrx.columns if x not in ['pos', 'chrom', 'geneid']]
pvmtrx = mtrx.loc[:, cols]

ase_samples_per_gene = pvmtrx[pvmtrx.loc[:, cols] < 5e-6].count(axis=1)
print(sum(ase_samples_per_gene.sort_values() > 200))
print("The number of individuals in which the gene is ASE:")
print("Maximum: ", ase_samples_per_gene.max())
print("Median:  ", ase_samples_per_gene.median())
print()

axe = ase_samples_per_gene.plot.hist(bins=200)
fig = axe.get_figure()
fig.savefig(wkdir + "/ase_samples_per_gene.pdf")
plt.close("all")

ase_genes_per_sample = pvmtrx[pvmtrx.loc[:, cols] < 5e-6].count(axis=0)
print("The number of ASE genes in one individual:")
print("Maximum: ", ase_genes_per_sample.max())
print("Median:  ", ase_genes_per_sample.median())

axe = ase_genes_per_sample.plot.hist(bins=200)
fig = axe.get_figure()
fig.savefig(wkdir + "/ase_gene_per_samples.pdf")
plt.close("all")
