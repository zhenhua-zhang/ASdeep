'''Generate ASE report.'''

import os
import sys
import math
import logging

from typing import Union
from itertools import accumulate

import pandas as pd
import gffutils as gut
import matplotlib.pyplot as plt

from .zutils import CHRLEN, fdr_bh

# Start position of each chromosome when binned by genomic coordination
xmax_list = accumulate(CHRLEN.values())
xmin_list = accumulate([0] + list(CHRLEN.values())[:-1])
CHRSPAN = dict(zip(CHRLEN.keys(), zip(xmin_list, xmax_list)))


class ASEReport:
    '''A class to process ASE report.'''

    def __init__(self, fpath_pool, save_pref, annot_db='ant_db_name.db', max_na_pergene=1500,
                 alpha=5e-2, ratio=0.10, black_list=None):
        self.save_pref = save_pref
        self.alpha = alpha
        self.ratio = ratio
        self.black_list = black_list

        self.genome_annot: gut.FeatureDB
        self.dtfm: pd.DataFrame
        self.pvm_raw: pd.DataFrame
        self.pvm_bycoord: pd.DataFrame

        self._pr_load_gtf(annot_db)
        self._pr_load_summary(fpath_pool)
        self._pr_pval_matrix(max_na_pergene)

    def _pr_load_gtf(self, annot_db):
        if not os.path.exists(annot_db):
            gut.create_db(annot_db.replace('.db', '.gtf'), self.genome_annot,
                          disable_infer_transcripts=True, disable_infer_genes=True)
        self.genome_annot = gut.FeatureDB(annot_db)

    @staticmethod
    def _pr_check_read_count(ase_rec: pd.Series, gross_count=5, allelic_count=3):
        allele_count_str = ase_rec['allele_counts']
        allele_count_lst = allele_count_str.split(':')[1].split(';')

        a1_sum, a2_sum = 0, 0
        for count_pair in allele_count_lst:
            a1_count, a2_count = count_pair.split('|')
            a1_sum += int(a1_count)
            a2_sum += int(a2_count)

        return (a1_sum>=allelic_count and a2_sum>=allelic_count) and ((a1_sum+a2_sum)>=gross_count)

    def _pr_load_summary(self, fpath_pool, sep='\t'):
        self.dtfm = pd.concat([pd.read_csv(fp, sep=sep, index_col=False) for fp in fpath_pool],
                              ignore_index=True)

        self.dtfm = self.dtfm.loc[self.dtfm.apply(self._pr_check_read_count, axis=1), :]

    def _pr_fetch_gene_coord(self, gene_id):
        gene = self.genome_annot[gene_id] # Feature attributes, accessed by gene.attributes['gene_name']
        return gene_id, gene.chrom, gene.start, gene.end

    def _pr_pval_matrix(self, max_na_pergene, pval_col='bn_pval'):
        # Get P values.
        self.pvm_raw = (self.dtfm
                        .loc[:, ['gene_id', 'sample_id', pval_col]]
                        .pivot(index='gene_id', columns='sample_id')
                        .drop(self.black_list, axis=0, errors='ignore')
                        .droplevel(0, axis=1))

        self.pvm_raw.columns.name = None
        self.pvm_raw.index.name = None

        # Remove genes without herterozygous locus in `max_na_pergene` individuals at maximum
        self.pvm_bycoord = (self.pvm_raw
                            .loc[self.pvm_raw.isna().sum(axis=1) < max_na_pergene, :]
                            .fillna(1))

        _coord = sorted([self._pr_fetch_gene_coord(gene_id) for gene_id in self.pvm_bycoord.index],
                            key=lambda x: (int(x[1]), int(x[2])))

        _, n_samples = self.pvm_raw.shape
        self.pvm_bycoord = (self.pvm_bycoord
                            .loc[[c[0] for c in _coord], :]
                            .apply(fdr_bh, axis=1, result_type='broadcast'))
        self.pvm_bycoord['ratio'] = self.pvm_bycoord.le(self.alpha).sum(axis=1) / n_samples
        self.pvm_bycoord['chrom'] = [c[1] for c in _coord]
        self.pvm_bycoord['pos'] = [c[2] for c in _coord]

    # Draw a heatmap of p-vales per gene across the genome.
    def _pr_draw_p_val_htmp(self):
        # Transform p-values by -log10()
        min_10_exp = sys.float_info.min_10_exp
        pval_matrix = (self.pvm_raw
                       .fillna(1)
                       .apply(lambda x: [-math.log10(e) if e > 0 else -min_10_exp for e in x]))

        # Heatmap
        figsize = [x / 2 if x < 200 else 100 for x in pval_matrix.shape]
        figsize[0], figsize[1] = figsize[1], figsize[0]

        fig, ax = plt.subplots()
        colorbar = ax.pcolor(pval_matrix)
        fig.colorbar(colorbar, ax=ax)

        return fig

    # Show the number of individuals is ASE for the gene.
    def _pr_draw_ase_gene_count_across_genome(self, figheight=9, figwidth=16):
        mtrx = (self.pvm_bycoord
                .loc[:, ['chrom', 'pos', 'ratio']]
                .apply(axis=1, result_type='expand',
                       func=lambda x: {'P': CHRSPAN[x['chrom']][0] + x['pos'], 'R': x['ratio']}))

        fig, axes = plt.subplots()
        axes.set_title('Percentage of ASE genes across the genome')
        axes.stem(mtrx.loc[:, 'P'], mtrx.loc[:, 'R'] * 100, linefmt='black', markerfmt='k.',
                  use_line_collection=True)

        chr_xticks, chr_labels = [], []
        for idx, (chrom, (xmin, xmax)) in enumerate(CHRSPAN.items()):
            chr_xticks.append(sum([xmin, xmax]) / 2)
            chr_labels.append(chrom)

            color = '0.5' if idx % 2 == 0 else '1'
            axes.axvspan(xmin=xmin, xmax=xmax, ymin=0, facecolor=color, alpha=0.2)

        axes.set_xticks(chr_xticks)
        axes.set_xticklabels(chr_labels, rotation=45, rotation_mode='anchor', ha='right')
        axes.set_xlim((-CHRLEN['1'] / 3, CHRSPAN['22'][1] + CHRLEN['1'] / 3))

        axes.set_xlabel('Genome coordination (gene binned by pos., {})'.format(mtrx.shape[0]))
        axes.set_ylabel('Percentage (out of {})'.format(self.pvm_bycoord.shape[1] - 2))
        axes.set_ylim((0, 65))

        fig.set_figheight(figheight)
        fig.set_figwidth(figwidth)
        plt.tight_layout()

        return fig

    def save_pval_matrices(self, report_fmt='csv'):
        '''Generate files to show quantification results.'''
        logging.info('Report')
        save_pref = self.save_pref
        if self.pvm_raw is not None:
            if report_fmt == 'csv':
                sep = ','
            elif report_fmt == 'tsv':
                sep = '\t'
            else:
                logging.warning('Unknown format: {}, using csv as default'.format(report_fmt))
                sep = ','

            self.pvm_raw.to_csv(save_pref + '-pvm_raw.' + report_fmt, sep=sep)
            self.pvm_bycoord.to_csv(save_pref + '-pvm_sorted_filtered_adj.' + report_fmt, sep=sep)


        return self

    def save_figures(self, fig_fmt: str = 'pdf', figheight=9, figwidth=16):
        '''Draw figures to show the result.'''
        logging.info('Visualize')
        save_pref = self.save_pref

        if self.pvm_raw is not None:
            (self._pr_draw_ase_gene_count_across_genome(figheight, figwidth)
             .savefig(save_pref + '-ase_gene_count_across_genome.' + fig_fmt))
            plt.close('all')

            pval_cols = [x for x in self.pvm_bycoord.columns if x not in ['ratio', 'chrom', 'pos']]
            ((self.pvm_bycoord.loc[:, pval_cols] < self.alpha)
             .sum(axis=0)
             .plot(kind='density')
             .set_title('ASE genes per individual (p<{})'.format(self.alpha))
             .get_figure()
             .savefig(save_pref + '-ase_gene_per_individual.' + fig_fmt))
            plt.close('all')

            ((self.pvm_bycoord.loc[:, pval_cols] < self.alpha)
             .sum(axis=1)
             .plot(kind='density')
             .set_title('Individuals per ASE gene (p<{})'.format(self.alpha))
             .get_figure()
             .savefig(save_pref + '-individual_per_ase_gene.' + fig_fmt))
            plt.close('all')

        return self

