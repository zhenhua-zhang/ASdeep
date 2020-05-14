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

logger = logging.getLogger("ASEReport")
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
        self.p_val_matrix_raw = None
        self.p_val_matrix_sorted_by_coord = None

        self.genome_annot = kwargs.get("genome_annot", None)
        self.ant_db_name = None
        self.ant_sql = None

    @staticmethod
    def _pr_close_fig():
        plt.cla()
        plt.clf()
        plt.close("all")

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

    def _pr_fetch_gene_coord(self, gene_id):
        gene = self.ant_sql[gene_id]  # Feature attributes, accessed by gene.attributes["gene_name"]
        return gene_id, gene.chrom, gene.start, gene.end

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

            gene_coord_list = [self._pr_fetch_gene_coord(gene_id) for gene_id in p_val_matrix.index]
            gene_coord_list = sorted(gene_coord_list, key=lambda x: (int(x[1]), int(x[2])))

            gene_id_sorted_by_coord_list = [coord[0] for coord in gene_coord_list]
            gene_chrom = [coord[1] for coord in gene_coord_list]

            p_val_matrix_sorted_by_coord = p_val_matrix.loc[gene_id_sorted_by_coord_list, :]
            p_val_matrix_sorted_by_coord["chrom"] = gene_chrom

            return p_val_raw_matrix, p_val_matrix, p_val_matrix_sorted_by_coord

        return None, None, None

    def _pr_draw_p_val_htmp(self):
        # Transform p-values by -log10()
        p_val_matrix = self.p_val_matrix.apply(lambda x: [math.log10(e) * -1 if e > 1e-6 else 7 for e in x])

        # Heatmap
        fig_size = [x / 2 if x < 200 else 100 for x in p_val_matrix.shape]
        fig_size[0], fig_size[1] = fig_size[1], fig_size[0]
        try:
            grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens", row_cluster=True, col_cluster=True)
        except:
            logger.warning("Failed to do cluster for row or columns, try again without cluster")
            grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens")

        return grid

    def _pr_draw_ase_gene_pp_across_genome(self):
        mtrx = self.p_val_matrix_sorted_by_coord

        chrom_width_list = mtrx["chrom"].value_counts()
        chrom_width_list.index = [int(x) for x in chrom_width_list.index]
        chrom_width_list= chrom_width_list.sort_index()

        chrom_vspan = [
            [0, xmax] if idx == 0 else [sum(chrom_width_list[:idx]), sum(chrom_width_list[:idx+1])]
            for idx, xmax in enumerate(chrom_width_list)
        ]

        chrom_xticks = [
            xmax / 2 if idx == 0 else (sum(chrom_width_list[:idx]) + sum(chrom_width_list[:idx+1])) / 2
            for idx, xmax in enumerate(chrom_width_list)
        ]

        mtrx = mtrx.loc[:, [x for x in mtrx.columns if x not in [("chrom", "")]]]
        ase_gene_pp = mtrx[mtrx<0.05].count(axis=1) / (mtrx.shape[1] - 1) * 100
        fig, axes = plt.subplots()
        axes.plot(ase_gene_pp.index, ase_gene_pp, linewidth=0.5)
        
        for idx, (xmin, xmax) in enumerate(chrom_vspan):
            color = "g" if idx % 2 == 0 else "r"
            axes.axvspan(xmin=xmin, xmax=xmax, ymin=0, facecolor=color, alpha=0.1)

        axes.set_xticks(chrom_xticks)
        axes.set_xticklabels(["chr{}".format(x) for x in chrom_width_list.index], rotation=45, rotation_mode="anchor", ha="right")

        return fig

    def _pr_desc_stat(self, adj_func=fdrcorrection):
        """Generate descriptive statistic for the merged ASE report.
        """
        if self.p_val_matrix_raw is not None:
            p_val_raw_matrix = self.p_val_matrix_raw \
                    .fillna(1) \
                    .apply(lambda x: adj_func(x)[-1], axis=0)

            ase_gene_per_individual = p_val_raw_matrix[p_val_raw_matrix < 0.05].count(axis=0)
            individual_per_ase_gene = p_val_raw_matrix[p_val_raw_matrix < 0.05].count(axis=1)

            return ase_gene_per_individual, individual_per_ase_gene

        return None, None

    def init(self, **kwargs):
        """Initialize the instance.
        """
        logger.info("Init")
        # Create a GFF/GTF database to fetch coordnation of given gene.
        ant_db_name = kwargs.get("ant_db_name", None)
        max_na_per_gene = kwargs.get("max_na_per_gene", 100)

        if ant_db_name is None:
            ant_db_name = os.path.splitext(self.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gut.create_db(self.genome_annot, ant_db_name, disable_infer_transcripts=True, disable_infer_genes=True)

        self.ant_sql = gut.FeatureDB(ant_db_name)
        self.merged_dtfm = self._pr_join_dtfm()
        self.p_val_matrix_raw, self.p_val_matrix, self.p_val_matrix_sorted_by_coord = self._pr_p_val_matrix(max_na_per_gene=max_na_per_gene)

        return self

    def report(self, report_fmt="csv"):
        """Generate files to show quantification results.
        """
        logger.info("Report")
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

            p_val_matrix_sorted_by_coord_output_name = ".".join([save_prefix, "p_value_matrix_sorted_by_coord", report_fmt])
            self.p_val_matrix_sorted_by_coord.to_csv(p_val_matrix_sorted_by_coord_output_name, sep=sep)

            p_value_raw_matrix_output_name = ".".join([save_prefix, "p_value_raw_matrix", report_fmt])
            self.p_val_matrix_raw.to_csv(p_value_raw_matrix_output_name, sep=sep)

        return self

    def visualize(self, fig_fmt="pdf"):
        """Draw figures to show the result.
        """
        logger.info("Visualize")
        save_prefix = self.save_prefix

        if self.p_val_matrix is not None:
            try:
                htmp_output = ".".join([save_prefix, "heatmap", fig_fmt])
                grid = self._pr_draw_p_val_htmp()
                grid.fig.savefig(htmp_output)
                plt.close("all")
            except Exception as exp:
                logger.error(exp)

            fig = self._pr_draw_ase_gene_pp_across_genome()
            fig.set_figheight(9)
            fig.set_figwidth(16)
            fig.savefig(".".join([save_prefix, "ase_gene_pp_across_genome", fig_fmt]))
            plt.cla()
            plt.clf()
            plt.close("all")

        if self.p_val_matrix_raw is not None:
            ase_gene_per_individual, individual_per_ase_gene = self._pr_desc_stat()

            axes = ase_gene_per_individual.plot(kind="hist")
            fig = axes.get_figure()
            fig.savefig(".".join([save_prefix, "ase_gene_per_individual", fig_fmt]))
            plt.cla()
            plt.clf()
            plt.close("all")

            individual_per_ase_gene_hist_fit = ".".join([save_prefix, "individual_per_ase_gene", fig_fmt])
            axes = individual_per_ase_gene.plot(kind="hist")
            fig = axes.get_figure()
            fig.savefig(individual_per_ase_gene_hist_fit)
            plt.cla()
            plt.clf()
            plt.close("all")

        return self
