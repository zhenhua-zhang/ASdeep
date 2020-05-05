"""ASE report.
"""

import logging
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

import gffutils as gut
import seaborn as sbn

logger = logging.getLogger("ASEFactory")
logger.setLevel(logging.DEBUG)

fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}",
                        datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_stream = logging.StreamHandler(sys.stderr)
cs_stream.setFormatter(fmt)
cs_stream.setLevel(logging.DEBUG)
logger.addHandler(cs_stream)


class ASEReport:
    """A class to process ASE report.
    """

    def __init__(self, file_path_list, save_prefix, **kwargs):
        self.file_path_list = file_path_list
        self.save_prefix = save_prefix

        self.merged_dtfm = None
        self.p_val_matrix = None
        self.p_val_raw_matrix = None

        self.genome_annot = kwargs.get("genome_annot", None)
        self.ant_db_name = None
        self.ant_sql = None

    def _pr_join_dtfm(self, **kwargs):
        """Merge multiple data frame of ASE quantification report.

        Parameters
        ------
        kwargs (optional): Parameters(i.e. sep, usecols) used in pandas.read_csv()
        """
        sep = kwargs.get("sep", "\t")
        use_cols = kwargs.get("usecols", None)
        dtfm_list = [pd.read_csv(file_path, sep=sep, usecols=use_cols) for file_path in self.file_path_list]

        if dtfm_list:
            return pd.concat(dtfm_list, ignore_index=True)

        logger.warning("The dataframe list is empty.")
        return None

    @staticmethod
    def _pr_check_read_count(ase_rec: pd.Series, gross_count=10, gross_count_per_allele=3):
        allele_count_str = ase_rec["allele_counts"]
        allele_count_lst = allele_count_str.split(":")[1].split(";")

        hete_num = len(allele_count_lst)
        a1_sum, a2_sum = 0, 0
        for count_pair in allele_count_lst:
            a1_count, a2_count = count_pair.split("|")
            a1_sum += int(a1_count)
            a2_sum += int(a2_count)

        return (a1_sum >= gross_count_per_allele and a2_sum >= gross_count_per_allele) or ((a1_sum + a2_sum) > gross_count)

    def _pr_p_val_matrix(self, max_na_per_gene=100, adj_func=fdrcorrection):
        """Get P values."""
        if self.merged_dtfm is not None:
            p_val_raw_matrix = self.merged_dtfm \
                    .loc[:, ["gene_id", "sample_id", "p_val"]] \
                    .pivot(index="gene_id", columns="sample_id")
            # Remove genes without herterozygous locus in `max_na_per_gene` individuals at maximum
            col_nona_count = p_val_raw_matrix.notna().sum(axis=1)
            p_val_matrix = p_val_raw_matrix \
                    .loc[col_nona_count < max_na_per_gene, :] \
                    .fillna(1) \
                    .apply(lambda x: adj_func(x)[-1], axis=0)

            return p_val_raw_matrix, p_val_matrix

        return None, None

    def _pr_draw_p_val_htmp(self):
        # Transform p-values by -log10()
        p_val_matrix = self.p_val_matrix.apply(lambda x: [math.log10(e) * -1 if e > 1e-6 else 7 for e in x])

        # Heatmap
        fig_size = [x / 2 if x < 200 else 100 for x in p_val_matrix.shape]
        fig_size[0], fig_size[1] = fig_size[1], fig_size[0]
        logger.debug(p_val_matrix.shape)
        try:
            grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens", row_cluster=True, col_cluster=True)
        except:
            logger.warning("Failed to do cluster for row or columns, try again without cluster")
            grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens")

        return grid

    def _pr_fetch_gene_coord(self, gene_id):
        gene = self.ant_sql[gene_id]  # Feature attributes, accessed by gene.attributes["gene_name"]
        return gene_id, gene.chrom, gene.start, gene.end

    def _pr_draw_p_val_mhtn(self):
        gene_id_list = self.p_val_matrix.index
        gene_coord_list = [self._pr_fetch_gene_coord(gene_id) for gene_id in gene_id_list]

        gene_coord_sorted_list = sorted(gene_coord_list, key=lambda x: (int(x[1]), int(x[2])))
        gene_id_sorted_by_coord_list = [coord[0] for coord in gene_coord_sorted_list]
        p_val_sorted_by_coord_matrix = self.p_val_matrix \
                .loc[gene_id_sorted_by_coord_list, :] \
                .apply(lambda x: [math.log10(e) * -1 if e > 1e-6 else 7 for e in x])

        fig_size = [x / 2 if x < 200 else 100 for x in p_val_sorted_by_coord_matrix.shape]
        fig_size[0], fig_size[1] = fig_size[1], fig_size[0]

        axes = p_val_sorted_by_coord_matrix \
                .plot(figsize=fig_size, legend=False, rot=45)
        fig = axes.get_figure()

        return fig

    def _pr_desc_stat(self, adj_func=fdrcorrection):
        """Generate descriptive statistic for the merged ASE report.
        """
        if self.p_val_raw_matrix is not None:
            p_val_raw_matrix = self.p_val_raw_matrix \
                    .fillna(1) \
                    .apply(lambda x: adj_func(x)[-1], axis=0)

            ase_gene_per_individual = p_val_raw_matrix[p_val_raw_matrix < 0.05].count(axis=0)
            individual_per_ase_gene = p_val_raw_matrix[p_val_raw_matrix < 0.05].count(axis=1)

            return ase_gene_per_individual, individual_per_ase_gene

        return None, None

    def init(self, **kwargs):
        """Initialize the instance.
        """
        # Create a GFF/GTF database to fetch coordnation of given gene.
        ant_db_name = kwargs.get("ant_db_name", None)
        max_na_per_gene = kwargs.get("max_na_per_gene", 100)

        if ant_db_name is None:
            ant_db_name = os.path.splitext(self.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gut.create_db(self.genome_annot, ant_db_name, disable_infer_transcripts=True, disable_infer_genes=True)

        self.ant_sql = gut.FeatureDB(ant_db_name)
        self.merged_dtfm = self._pr_join_dtfm()
        self.p_val_raw_matrix, self.p_val_matrix = self._pr_p_val_matrix(max_na_per_gene=max_na_per_gene)

        return self

    def report(self, report_fmt="csv"):
        """Generate files to show quantification results.
        """
        save_prefix = self.save_prefix
        if self.p_val_matrix is not None:
            if report_fmt == "csv":
                sep = ","
            elif report_fmt == "tsv":
                sep = "\t"
            else:
                logger.warning("Unknown type of format: {}, using csv as default".format(report_fmt))
                sep = ","

            p_value_matrix_output_name = ".".join([save_prefix, "p_value_matrix", report_fmt])
            self.p_val_matrix.to_csv(p_value_matrix_output_name, sep=sep)

            p_value_raw_matrix_output_name = ".".join([save_prefix, "p_value_raw_matrix", report_fmt])
            self.p_val_raw_matrix.to_csv(p_value_raw_matrix_output_name, sep=sep)

        return self

    def visualize(self, fig_fmt="pdf"):
        """Draw figures to show the result.
        """
        save_prefix = self.save_prefix
        if self.p_val_matrix is not None:
            try:
                htmp_output = ".".join([save_prefix, "heatmap", fig_fmt])
                grid = self._pr_draw_p_val_htmp()
                grid.fig.savefig(htmp_output)
                plt.close("all")
            except Exception as exp:
                logger.error(exp)

            mhtn_output = ".".join([save_prefix, "manhatten", fig_fmt])
            fig = self._pr_draw_p_val_mhtn()
            fig.savefig(mhtn_output)
            plt.close("all")

        if self.p_val_raw_matrix is not None:
            ase_gene_per_individual, individual_per_ase_gene = self._pr_desc_stat()

            ase_gene_per_individual_hist_fig = ".".join([save_prefix, "ase_gene_per_individual", fig_fmt])
            axes = ase_gene_per_individual.plot(kind="hist")
            fig = axes.get_figure()
            fig.savefig(ase_gene_per_individual_hist_fig)
            plt.close("all")

            individual_per_ase_gene_hist_fit = ".".join([save_prefix, "individual_per_ase_gene", fig_fmt])
            axes = individual_per_ase_gene.plot(kind="hist")
            fig = axes.get_figure()
            fig.savefig(individual_per_ase_gene_hist_fit)
            plt.close("all")

        return self
