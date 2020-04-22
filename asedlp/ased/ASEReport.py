"""ASE report.
"""

import logging
import math
import os
import sys

import matplotlib.pyplot as plt
import pandas as pds
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
        dtfm_list = [pds.read_csv(file_path, sep=sep, usecols=use_cols) for file_path in self.file_path_list]

        if dtfm_list:
            return pds.concat(dtfm_list, ignore_index=True)

        logger.warning("The dataframe list is empty.")
        return None

    def _pr_p_val_matrix(self, max_na=150, adj_func=fdrcorrection):
        """Get P values."""
        p_val_matrix = self.merged_dtfm \
                .loc[:, ["gene_id", "sample_id", "p_val"]] \
                .pivot(index="gene_id", columns="sample_id")
        # Remove genes without herterozygous locus in `max_na` individuals at maximum
        col_nona_count = p_val_matrix.notna().sum(axis=1)
        p_val_matrix = p_val_matrix \
                .loc[col_nona_count < max_na, :] \
                .fillna(1) \
                .apply(lambda x: adj_func(x)[-1], axis=0)

        return p_val_matrix

    def _pr_draw_p_val_htmp(self, max_na=150):
        # Transform p-values by -log10()
        merged_dtfm = self.p_val_matrix.apply(lambda x: [math.log10(e) * -1 if e > 1e-6 else 7 for e in x])

        # Heatmap
        fig_size = [x / 2 if x < 200 else 100 for x in merged_dtfm.shape]
        fig_size[0], fig_size[1] = fig_size[1], fig_size[0]
        grid = sbn.clustermap(merged_dtfm, figsize=fig_size, cmap="Greens", row_cluster=True, col_cluster=True)
        return grid

    def _pr_fetch_gene_coord(self, gene_id):
        gene = self.ant_sql[gene_id]
        return gene_id, gene.name, gene.chrom, gene.start, gene.end

    def _pr_draw_p_val_mhtn(self, max_na=150):
        gene_id_list = self.p_val_matrix \
                .loc[:, "gene_id"] \
                .drop_duplications()

        gene_coord_list = [self._pr_fetch_gene_coord(gene_id) for gene_id in gene_id_list]

        # XXX: The chromosome should be numerical
        gene_coord_list = sorted(gene_coord_list, key=lambda x: (int(x[2]), int(x[3])))
        gene_id_name_list = [[coord[0], coord[1]] for coord in gene_coord_list]
        fig, ax = plt.subplots()
        ax.plot()

        return self

    def _pr_desc_stat(self):
        """Generate descriptive statistic for the merged ASE report.
        """
        return 0

    def init(self, ant_db_name=None):
        """Initialize the instance.
        """
        # Create a GFF/GTF database to fetch coordnation of given gene.
        if ant_db_name is None:
            ant_db_name = os.path.splitext(self.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gut.create_db(self.genome_annot, ant_db_name, disable_infer_transcripts=True, disable_infer_genes=True)

        self.ant_sql = gut.FeatureDB(ant_db_name)
        self.merged_dtfm = self._pr_join_dtfm()
        logger.debug(self.merged_dtfm.columns)
        self.p_val_matrix = self._pr_p_val_matrix()
        logger.debug(self.p_val_matrix.columns)

        return self

    def report(self, output_pre="output", report_fmt="csv"):
        """Generate files to show quantification results.
        """
        if report_fmt == "csv":
            sep = ","
        elif report_fmt == "tsv":
            sep = "\t"
        else:
            logger.warning("Unknown type of format: {}, using csv as default".format(report_fmt))
            sep = ","

        p_value_matrix_output_name = ".".join([output_pre, "p_value_matrix", report_fmt])
        self.p_val_matrix.to_csv(p_value_matrix_output_name, sep=sep)

        return self

    def visualize(self, output_pre="output", fig_fmt="pdf"):
        """Draw figures to show the result.
        """
        htmp_output = ".".join([output_pre, "heatmap", "pdf"])
        grid = self._pr_draw_p_val_htmp()
        grid.fig.savefig(htmp_output)

        return self
